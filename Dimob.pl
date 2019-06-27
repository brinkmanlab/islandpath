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

    December 16th, 2016

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
    my $usage = "Usage:\n./Dimob.pl <genome.gbk> <outputfile.txt> [cutoff_dinuc_bias == int] [min_length == int] \nExample:\n./Dimob.pl example/NC_003210.gbk NC_003210_GIs.txt 8 8000\n";

    my ($inputfile, $outputfile, $cutoff_dinuc_bias, $min_length) = @ARGV;

    # Check that input file and output file are specified or die and print help message
    unless(defined($inputfile) && defined($outputfile)){
        print $usage;
        exit;
    }
    
    ## min_length
	if (!$min_length) { $min_length = 8000; }
    
    # Create a dimob object
    my $dimob_obj = Dimob->new(
        {
        cfg_file => $cfname,
        bindir => $RealBin,
        workdir => $cwd,
        MIN_GI_SIZE => $min_length,
        extended_ids => 1
        }
    );

	# Recover the config from file, initialized during creation dimob_obj
    my $cfg = Dimob::Config->config;
    $cfg->{logger_conf} = $RealBin."/".$cfg->{logger_conf};
    $cfg->{hmmer_db} = "$RealBin/".$cfg->{hmmer_db};

    # Check that the logger exists and initializes it
    #print $cfg->{logger_conf};
    if($cfg->{logger_conf} && ( -r $cfg->{logger_conf})) {
        Log::Log4perl::init($cfg->{logger_conf});
        $logger = Log::Log4perl->get_logger;
        #$logger->debug("Logging initialized");
    }

    $logger->debug("IslandPath-DIMOB initialized");
	
	## min_length
	if (!$min_length) { 
		$min_length = 8000;
	    $logger->debug("Use default min_length: 8000 bp");
	} else {
	    $logger->debug("Use min_length: $min_length");
	}

    # Transform relative path to absolute path and 
    # check that input file is readable
    $inputfile = File::Spec -> rel2abs($inputfile);
    unless( -f $inputfile && -r $inputfile ) {
        print "Error: $inputfile is not readable";
        exit;
    }

    ## check if $cutoff_genes provided
	if (!$cutoff_dinuc_bias) {
		$cutoff_dinuc_bias = 8;
        $logger->debug("Use cutoff_dinuc_bias default: 8");
	} else {
        $logger->debug("Use cutoff_dinuc_bias provided: ".$cutoff_dinuc_bias);
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
    my @islands = $dimob_obj->run_dimob($inputfile, $outputfile, $cutoff_dinuc_bias);

    $logger->info("Printing results");

    my $i = 1;
    open my $fhgd, '>', $outputfile or die "Cannot open $outputfile: $!";
	my $discard_file = $outputfile."_discard.txt";
	open DISCARD, '>', $discard_file or die "Cannot open $discard_file: $!";
	print DISCARD "seq\tstart\tend\tlength\n";
    foreach my $island (@islands) {
        my $start = $island->[0];
        my $end = $island->[1];
        my $seq = $island->[2];
        
        ## discard if smaller than min length set
        my $diff = $end - $start;        
        if ($diff < $min_length) {
        	print DISCARD "$seq\t$start\t$end\t$diff\n";
        } else {
	        print $fhgd "GI_$i\t$seq\t$start\t$end\n";
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

## 
## New implementations
## 
## Add multicontig function
## Fix smartmacth experimental warning message
## Fix dinuc bias loop iteration bug
## Use cutoff_dinuc_bias as a variable
## Output csv information for dinucleotide bias.
## Provide additional information in output GI information
## 
