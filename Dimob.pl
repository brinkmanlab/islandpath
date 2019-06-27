#!/usr/bin/env perl
=head1 NAME

    IslandPath-DIMOB

=head1 DESCRIPTION

    Script to run IslandPath-DIMOB genomic island prediction on a given genome

=head1 SYNOPSIS

	perl Dimob.pl <genome.gbk> <output_name> [cutoff_dinuc_bias] [min_length]

	Default values:
		cutoff_dinuc_bias = 8
		min_length = 8000

	Example:
		perl Dimob.pl example/NC_003210.gbk NC_003210_GIs
		perl Dimob.pl example/NC_003210.gbk NC_003210_GIs 6 10000
		perl Dimob.pl example/NC_000913.embl NC_000913_GIs 6 10000

=head1 AUTHORS

	Claire Bertelli [Original author]
	Email: claire.bertelli@sfu.ca
    Brinkman Laboratory
    Simon Fraser University
    
   	Jose F. Sanchez-Herrero [Developer of this fork]
	Email: jsanchez@igtp.cat
	Bioinformatics Facility Unit, Institut German Trias i Pujol (IGTP) 
	Badalona, Barcelona, Spain	


=head1 LAST MAINTAINED

    June 27th, 2019

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

## 
## New implementations
## 
## Add multicontig function
## Fix smartmacth experimental warning message
## Fix dinuc bias loop iteration bug
## Use GI min_length as a variable
## Use cutoff_genes_dinuc as a variable
## merge with gff branch 
## Output csv information for dinucleotide bias.
## Provide additional information in output GI information

MAIN: {

    # config files
    my $cwd = getcwd;
    my $cfname = "$RealBin/Dimob.config";
    my $logger;
    #my $logger_conf = "$RealBin/logger.conf";

    # usage help
    my $usage = "Usage:\nperl Dimob.pl <genome.gbk> <output_name> [cutoff_dinuc_bias] [min_length]\n";
    $usage .= "\nDefault values:\n\tcutoff_dinuc_bias = 8\n\tmin_length = 8000\n\n";
    $usage .= "Example:\n\tperl Dimob.pl example/NC_003210.gbk NC_003210_GIs\n";
    $usage .= "\tperl Dimob.pl example/NC_003210.gbk NC_003210_GIs 6 10000\n";
    $usage .= "\tperl Dimob.pl example/NC_000913.embl NC_000913_GIs 6 10000\n\n";

    my ($inputfile, $output_name, $cutoff_dinuc_bias, $min_length) = @ARGV;

    # Check that input file and output file are specified or die and print help message
    unless(defined($inputfile) && defined($output_name)){
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
    my @islands = $dimob_obj->run_dimob($inputfile, $output_name, $cutoff_dinuc_bias);

    $logger->info("Printing results");

	my $outputfile = $output_name.".txt";
    ## txt output
    open OUT_TXT, '>', $outputfile or die "Cannot open $outputfile: $!";
	print OUT_TXT "##GI\tseq\tstart\tend\tstrand\n";
		
	## gff output
	my $gff_file = $output_name.".gff3";
	open GFF, '>', $gff_file or die "Cannot open $gff_file: $!";
	print GFF "##gff-version 3\n";

	## discarded regions
	my $discard_file = $output_name."_discard.txt";
	open DISCARD, '>', $discard_file or die "Cannot open $discard_file: $!";
	print DISCARD "##seq\tstart\tend\tlength\tstrand\n";

	## loop through islands and print
    my $i = 1;
    foreach my $island (@islands) {
        my $seq = $island->[0];
        my $start = $island->[1];
        my $end = $island->[2];
        my $strand = $island->[3];
        
        ## discard if smaller than min length set
        my $diff = $end - $start;        
        if ($diff < $min_length) {
        	print DISCARD "$seq\t$start\t$end\t$diff\t$strand\n";
        } else {

		    #$logger->info("Warning: txt output is now depreciated. Support has been added to output GFF3 formatted documents. Use (any) other extension to enable GFF output. See: https://github.com/brinkmanlab/islandpath/issues/7");
	        print OUT_TXT "GI_$i\t$seq\t$start\t$end\t$strand\n";
	        print GFF "$seq\tislandpath\tgenomic_island\t$start\t$end\t.\t$strand\t.\tID=$seq\_gi$i\n";
	        $i++;
        }
    }
	## close filehandles
    close (GFF); close (OUT_TXT); close (DISCARD);
	
	## finish
    $logger->info("Removing tmp files");
 	unless(unlink glob "$inputfile.*") {
        $logger->error("Can't remove $inputfile: $!");
    }
    unless(rmdir $tmp_path) {
        $logger->error("Can't remove $tmp_path: $!");
    }
}

