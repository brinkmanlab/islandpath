#!/usr/bin/env perl
=pod

    IslandPath-DIMOB

    Script to run IslandPath-DIMOB genomic island prediction on a given genome

    Usage:
    ./Dimob.pl <genome.gbk> <outputfile.txt>
    Example:
    ./Dimob.pl example/NC_003210.gbk NC_003210_GIs.txt

    Claire Bertelli
    Brinkman Laboratory
    Simon Fraser University
 	Email: claire.bertelli@sfu.ca

    Last maintained: October 24th, 2016

=cut


use strict;
use warnings FATAL => 'all';
use Getopt::Long;

use lib './lib';
use GenomeUtils;
use Dimob;

use Data::Dumper;
use File::Copy;
use File::Basename;
use File::Spec;
use File::Path;

MAIN: {

    # config files
    my $cfname = 'Dimob.config';
    my $logger;
#    my $logger_cfg = 'logger.conf';

    # usage help
    my $usage = "Usage:\n./Dimob.pl <genome.gbk> <outputfile.txt>\nExample:\n./Dimob.pl example/NC_003210.gbk NC_003210_GIs.txt\n";

#   my $arguments = GetOptions()
    my ($inputfile, $outputfile) = @ARGV;

    # Check that input file and output file are specified or die and print help message
    unless(defined($inputfile) && defined($outputfile)){
        print $usage;
        exit;
    }

    # Check that input file is readable
    unless( -f $inputfile && -r $inputfile ) {
        print "Error: $inputfile is not readable";
        exit;
    }

    # if it does not exist already, create a tmp directory to store intermediate results, copy the input file to the tmp
    #$logger->info("Creating temp directory with needed files");
    my $tmp_path = "tmp_dimob";
    if (! -d $tmp_path)
    {
        my $dirs = eval { mkpath($tmp_path) };
        die "Failed to create $tmp_path: $@\n" unless $dirs;
    }
    copy($inputfile,$tmp_path) or die "Failed to copy $inputfile: $!\n";
    my($filename, $dirs, $suffix) = fileparse($inputfile, qr/\.[^.]*/);
    $inputfile = File::Spec->catfile($tmp_path,$filename);

    # create a dimob object
    my $dimob_obj = Dimob->new({cfg_file => $cfname, workdir => './tmp_dimob',
            MIN_GI_SIZE => 2000}
    );
    my $cfg = Dimob::Config->config;


    # Check that the logger exists and initializes it
    if($cfg->{logger_conf} && ( -r $cfg->{logger_conf})) {
        Log::Log4perl::init($cfg->{logger_conf});
        $logger = Log::Log4perl->get_logger;
        $logger->debug("Logging initialized");
    }



    ######
    # From an embl or genbank file regenerate a ptt, ffn, and faa file needed by dimob.pm

    # create a genome object from package GenomeUtils
    my $genome_obj = GenomeUtils->new(
                                            {workdir => './tmp_dimob'}
                                            );

    # check the gbk/embl file format
    $genome_obj->read_and_check($genome_obj, $inputfile);

    # read the gbk/embl file and convert it to all files needed
    my $genome_name = $genome_obj->{'base_filename'};

    $genome_obj->read_and_convert($inputfile, $genome_name);

    ######
    # Runs IslandPath-DIMOB on the genome files

    $logger->info("Running IslandPath-DIMOB");
    my @islands = $dimob_obj->run_dimob($inputfile);

    $logger->info("Printing results");

    my $i = 1;
    open my $fhgd, '>', $outputfile or die "Cannot open output.txt: $!";
    foreach my $island (@islands) {
        my $start = $island->[0];
        my $end = $island->[1];
        print $fhgd "GI_$i\t$start\t$end\n";
        $i++;
    }
    close $fhgd;

    $logger->info("Removing tmp files");
    unless(unlink glob "$inputfile.*") {
        $logger->error("Can't remove $inputfile: $!");
    }
    unless(rmdir $tmp_path) {
        $logger->error("Can't remove $inputfile: $!");
    }

}
