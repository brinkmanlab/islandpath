#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
use Getopt::Long;

use lib './lib';
use GenomeUtils;
use Dimob;


my $usage = "Usage:\n./Dimob.pl <genome.gbk> <outputfile.txt>\nExample:\n./Dimob.pl example/NC_003210.gbk NC_003210_GIs.txt\n";

my ($inputfile, $outputfile) = @ARGV;

unless(defined($inputfile) && defined($outputfile)){
    print $usage;
    exit;
}

# create a genome object from package GenomeUtils
my $genome_obj = GenomeUtils->new(
                            			{workdir => './tmp'}
                            			);

# check the gbk/embl file format
$genome_obj->read_and_check($genome_obj, $inputfile);

# read the gbk/embl file and convert it to all files needed
my $genome_name = $inputfile;
$genome_obj->read_and_convert($inputfile, $genome_name);

# create a dimob object
my $dimob_obj = Dimob->new({workdir => './tmp',
                         MIN_GI_SIZE => 8000}
                         );

# run islandPath-Dimob
$inputfile =~ s/\/\//\//g;
my ($file, $extension ) = $inputfile =~ /(.+)\.(\w+)/;
$dimob_obj->run_dimob($file);

