#!/usr/bin/env perl
=head1 NAME

    IslandPath-DIMOB

=head1 DESCRIPTION

    Script to run IslandPath-DIMOB genomic island prediction on a given genome

=head1 SYNOPSIS

    ./Dimob.pl <genome.gbk> <outputfile.txt>
    Example:
    ./Dimob.pl example/NC_003210.gbk NC_003210_GIs.txt

=head1 AUTHORS

	Claire Bertelli
	Email: claire.bertelli@sfu.ca
    Brinkman Laboratory
    Simon Fraser University

=head1 LAST MAINTAINED

    January 16th, 2019

=cut


use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Data::Dumper;
use File::Copy;
use File::Basename;
use File::Spec;
use File::Path;
use Log::Log4perl qw(get_logger :nowarn);
use File::Temp qw/ :mktemp /;
use Cwd;

# use local Dimob libraries
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use GenomeUtils;
use Dimob;


MAIN: {

    # config files
    my $cwd = getcwd;
    my $cfname = "$RealBin/Dimob.config";
    my $logger;
    #my $logger_conf = "$RealBin/logger.conf";

    # usage help
    my $usage = "Usage:\n./Dimob.pl <genome.gbk> <outputfile.txt>\nExample:\n./Dimob.pl example/NC_003210.gbk NC_003210_GIs.txt\n";

    my ($inputfile, $outputfile) = @ARGV;

    # Check that input file and output file are specified or die and print help message
    unless(defined($inputfile) && defined($outputfile)){
        print $usage;
        exit;
    }

    # Transform relative path to absolute path and 
    # check that input file is readable
    $inputfile = File::Spec -> rel2abs($inputfile);
    unless( -f $inputfile && -r $inputfile ) {
        print "Error: $inputfile is not readable";
        exit;
    }

    # Create a dimob object
    my $dimob_obj = Dimob->new(
        {cfg_file => $cfname,
            bindir => $RealBin,
            workdir => $cwd,
            MIN_GI_SIZE => 2000,
            extended_ids => 1
        }
    );
    
    # Recover the config from file, initialized during creation dimob_obj
    my $cfg = Dimob::Config->config;
    $cfg->{logger_conf} = "$RealBin/" . $cfg->{logger_conf};
    $cfg->{hmmer_db} = "$RealBin/" . $cfg->{hmmer_db};

    # Check that the logger exists and initializes it
    print $cfg->{logger_conf};
    if($cfg->{logger_conf} && ( -r $cfg->{logger_conf})) {
        Log::Log4perl::init($cfg->{logger_conf});
        $logger = Log::Log4perl->get_logger;
        $logger->debug("Logging initialized");
    }


    # Create a tmp directory to store intermediate results, copy the input file to the tmp
    $logger->info("Creating temp directory with needed files");
    my($filename, $dirs, $suffix) = fileparse($inputfile, qr/\.[^.]+$/);

    my $tmp_path = mkdtemp($dirs . "dimob_tmpXXXXXXXXXX");
    if (! -d $tmp_path)
    {
        my $dirs = eval { mkpath($tmp_path) };
        die "Failed to create $tmp_path: $@\n" unless $dirs;
    }
    copy($inputfile,$tmp_path) or die "Failed to copy $inputfile: $!\n";
    $inputfile = File::Spec->catfile($tmp_path,$filename);

    # update workdir in genome_obj with the temporary directory
    $dimob_obj -> {workdir} = $tmp_path;

    ######
    # From an embl or genbank file regenerate a ptt, ffn, and faa file needed by dimob.pm

    # create a genome object from package GenomeUtils
    my $genome_obj = GenomeUtils->new();

    $logger->info("This is the $inputfile");
    # check the gbk/embl file format
    $genome_obj->read_and_check($genome_obj, $inputfile . $suffix);

    # read the gbk/embl file and convert it to all files needed
    my $genome_name = $genome_obj->{'base_filename'};

    $genome_obj->read_and_convert($inputfile . $suffix, $genome_name);

    ######
    # Runs IslandPath-DIMOB on the genome files

    $logger->info("Running IslandPath-DIMOB");
    my @islands = $dimob_obj->run_dimob($inputfile);

    $logger->info("Printing results");

    my $i = 1;
    open my $fhgd, '>', $outputfile or die "Cannot open output.txt: $!";

    if ($outputfile =~ /\.txt$/) {
        #legacy output
        $logger->info("Warning: txt output is now depreciated. Support has been added to output GFF3 formatted documents. Use (any) other extension to enable GFF output. See: https://github.com/brinkmanlab/islandpath/issues/7");
        foreach my $island (@islands) {
            my $start = $island->[1];
            my $end = $island->[2];
            print $fhgd "GI_$i\t$start\t$end\n";
            $i++;
        }
    } else {
        #GFF output
        print $fhgd "##gff-version 3\n";
        foreach my $island (@islands) {
            my $label = $island->[0];
            my $start = $island->[1];
            my $end = $island->[2];
            my $strand = $island->[3];
            #TODO use proper chromosome sequence id
            print $fhgd "$label\tislandpath\tgenomic_island\t$start\t$end\t.\t$strand\t.\tID=$label\_gi$i\n";
            $i++;
        }
    }

    close $fhgd;

    $logger->info("Removing tmp files");
    unless(unlink glob "$inputfile.*") {
        $logger->error("Can't remove $inputfile: $!");
    }
    unless(rmdir $tmp_path) {
        $logger->error("Can't remove $tmp_path: $!");
    }
}
