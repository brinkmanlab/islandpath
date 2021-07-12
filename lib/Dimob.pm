=head1 NAME

    Dimob

=head1 DESCRIPTION

    Object to run Dimob against a given genome

=head1 SYNOPSIS

    use Dimob;

    $dimob_obj = Dimob->new({workdir => '/tmp/workdir',
                                           MIN_GI_SIZE => 8000});

    $dimob_obj->run_dimob($rep_accnum);
    
=head1 AUTHORS

	Claire Bertelli
	Email: claire.bertelli@sfu.ca
    and
    Matthew Laird
    Email: lairdm@sfu.ca
    Brinkman Laboratory
    Simon Fraser University

=head1 LAST MAINTAINED

    Oct 24, 2016

=cut

package Dimob;

use strict;
use warnings;
use Moose;
use Log::Log4perl qw(get_logger :nowarn);
use File::Temp qw/ :mktemp /;
use Data::Dumper;
use File::Spec;


use Dimob::genomicislands;
use Dimob::mobgene;
use Dimob::Config;
use GenomeUtils;


my $cfg; my $cfg_file; my $logger;

my $module_name = 'Dimob';

sub BUILD {
    my $self = shift;
    my $args = shift;

    # Initialize the configuration file
    # CHECKING INITIALIZATION
    Dimob::Config->initialize({cfg_file => $args->{cfg_file}});
    $cfg = Dimob::Config->config;
    $cfg -> {rootdir} = $args -> {bindir};
    $cfg_file = File::Spec->rel2abs(Dimob::Config->config_file);

    # Initialize the logfile with the current date
    $logger = Log::Log4perl->get_logger;
    $self->log_rotate();
    $logger->trace("IslandPath-DIMOB initialized");

    die "Error, work dir not specified:  $args->{workdir}"
			unless( -d $args->{workdir} );
    $self->{workdir} = $args->{workdir};

    $self->{MIN_GI_SIZE} = $args->{MIN_GI_SIZE};

    # Do we need to use extended ids because
    # there could be duplicate gis. All the files 
    # we generate have coordinates in the faa header
    if($args->{extended_ids}) {
			$logger->trace("Using extended ids");
			$self->{extended_ids} = 1;
    }
    
}

# The generic run to be called from the scheduler
# magically do everything.

sub run {
    my $self = shift;
    my $filename = shift;
    my $callback = shift;

    my @islands = $self->run_dimob($filename);

    print Dumper \@islands;

### TO REMOVE
    if(@islands) {
			# If we get a undef set it doesn't mean failure, just
			# nothing found.  Write the results to the callback
			# if we have any
			if($callback) {
					$callback->record_islands($module_name, @islands);
			}
    }
### TO REMOVE


    # We just return 1 because any failure for this module
    # would be in the form of an exception thrown.
    return 1;

}

sub run_dimob {
    my $self = shift;
    my $filename = shift;
    my @tmpfiles;

    # We're given the filename, look up the files
    unless($filename) {
    	$logger->error("Error, can't find genome $filename");
    	return ();
    }    
#########################Need to look if file exists and size different from 0
    my $format_str = $cfg->{expected_exts};
    $logger->trace("Genome $filename, found formats: $format_str");

    # To make life easier, break out the formats available
    my $formats;
    foreach (split /\s+/, $format_str) { $_ =~ s/^\.//; $formats->{$_} = 1; }

    # Ensure we have the needed files
    unless($formats->{ffn}) {
    	$logger->error("Error, we don't have the needed ffn file...");
    	return ();
    }
    unless($formats->{faa}) {
    	$logger->error("Error, we don't have the needed faa file...");
    	return ();
    }
    unless($formats->{ptt}) {
    	$logger->error("Error, we don't have the needed ptt file...");
    	return ();
    }

    $logger->debug("Looking for mobility genes using hmmer");
    # We need a temporary file to hold the hmmer output
    my $hmmer_outfile = $self->_make_tempfile();
    push @tmpfiles, $hmmer_outfile;

    # Now the command and database to use....
    my $cmd = $cfg->{hmmer_cmd};
    my $hmmer_db = $cfg->{hmmer_db};
    $cmd .= " $hmmer_db $filename.faa >$hmmer_outfile";
    $logger->debug("Running hmmer command $cmd");
    my $rv = system($cmd);
		
    #	or $logger->logdie("Error running hmmer: $!");
		# if($rv != 0) {
		#		$logger->logdie("Error running hmmer, rv: $rv");
		# }

    unless( -s $hmmer_outfile ) {
    	$logger->logdie("Error, hmmer output seems to be empty");
    }

    my $mob_list;

    $logger->debug("Parsing hmmer results with Mobgene");
    my $mod_args = {};
    if($self->{extended_ids}) {
    	$mod_args->{extended_ids} = $self->{extended_ids};
    }

    my $mobgene_obj = Dimob::Mobgene->new($mod_args);
    my $mobgenes = $mobgene_obj->parse_hmmer( $hmmer_outfile, $cfg->{hmmer_evalue} );
    foreach(keys %$mobgenes){
    	$mob_list->{$_}=1;   
    }

    #get a list of mobility genes from ptt file based on keyword match
    $logger->debug("Retrieving mobility genes from ptt file");
    my $mobgene_ptt = $mobgene_obj->parse_ptt("$filename.ptt");

    foreach(keys %$mobgene_ptt){
    	$mob_list->{$_}=1;   
    }

    #calculate the dinuc bias for each gene cluster of 6 genes
    #input is a fasta file of ORF nucleotide sequences
    $logger->debug("Calculating dinucleotide bias");
    my $dinuc_results = cal_dinuc("$filename.ffn");
    my @dinuc_values;
    foreach my $val (@$dinuc_results) {
    	push @dinuc_values, $val->{'DINUC_bias'};
    }

    #calculate the mean and std deviation of the dinuc values
    my $median = cal_median( \@dinuc_values );
    my $sd   = cal_stddev( \@dinuc_values );

    #generate a list of dinuc islands with ffn fasta file def line as the hash key
    my $gi_orfs = dinuc_islands( $dinuc_results, $median, $sd, 8 );

    #convert the def line to gi numbers (the data structure is maintained)
    my $extended = $self->{extended_ids} ? 1 : undef;
    my $dinuc_islands = defline2gi( $gi_orfs, "$filename.ptt", $extended );

    #check the dinuc islands against the mobility gene list
    #any dinuc islands containing >=1 mobility gene are classified as
    #dimob islands
    $logger->debug("Looking for regions with dinuc bias and mobility genes");
    my $dimob_islands = dimob_islands( $dinuc_islands, $mob_list );

    my @gis;
    foreach (@$dimob_islands) {

    	#get the pids from the  for just the start and end genes
    	unless($_->[0]{start} && $_->[-1]{end}) {
    		$logger->warn("Warning, GI is missing either start or end: ($_->[0]{start}, $_->[-1]{end})");
    		next;
    	}

			push (@gis, [ $_->[0]{seq_id}, $_->[0]{start}, $_->[-1]{end}, $_->[0]{strand}]);
			#my $start = $_->[0]{start};
			#my $end = $_->[-1]{end};		 
			#print "$start\t$end\n";    
    }
		
    # And cleanup after ourself
    $logger->debug("Cleaning temporary files");

    if($cfg->{clean_tmpfiles}) {
        $logger->trace("Cleaning up temp files for Dimob");
        $self->_remove_tmpfiles(@tmpfiles);
    }

    return @gis;
}

# Make a temp file in our work directory and return the name

sub _make_tempfile {
	my $self = shift;

	# Let's put the file in our workdir
	my $tmp_file = mktemp($self->{workdir} . "/blasttmpXXXXXXXXXX");
	
	# And touch it to make sure it gets made
	`touch $tmp_file`;

	return $tmp_file;
}

sub _remove_tmpfiles {
	my $self = shift;
	my @tmpfiles = @_;

	foreach my $file (@tmpfiles) {
		unless(unlink $file) {
		$logger->error("Can't unlink file $file: $!");
		}
	}
}

1;

sub log_rotate {
    my $self = shift;

    # Build the logfile name we want to use
    my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
    my $datestr = "$year$mon$mday-$hour$min$sec";
    my $logfile = File::Spec->catpath(undef, $cfg->{logdir}, "dimob.$datestr.log");

    my $app = Log::Log4perl->appender_by_name("errorlog");

    unless($app) {
        $logger->warn("Logging doesn't seem to be defined, you might not even see this");
        return;
    }

    if($app->filename eq $logfile) {
        $logger->info("Asked to rotate the logfile, nothing to do, keeping " . $app->filename);
    } else {
        $logger->info("Rotating log file, current file: " . $app->filename . ", switching to: " . $logfile);
        $app->file_switch($logfile);
        $logger->info("Initializing logfile, rotated.");
    }
}
