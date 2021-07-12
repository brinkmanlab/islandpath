=head1 DESCRIPTION

    A repackaging of Will's original Dimob software
    loosely associated subroutines used by Dimob.pm for generating dinuc and dimob islands

=head1 SYNOPSIS

    use Dimob::genomicislands;

	see Dimob.pm

=head1 AUTHOR

    Claire Bertelli
    Email: claire.bertelli@sfu.ca
    and
    Matthew Laird
    Email: lairdm@sfu.ca
    Brinkman Laboratory
    Simon Fraser University

=head1 LAST MAINTAINED

    Dec 17, 2016

=cut


package Dimob::genomicislands;

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::SeqWords2;
use File::Basename;
use Statistics::Descriptive;
use Getopt::Long;
use Carp;

use Dimob::tabdelimitedfiles;

our ( @ISA, @EXPORT, @EXPORT_OK );
use Exporter;
@ISA       = qw(Exporter);
@EXPORT    = qw(cal_dinuc cal_mean cal_median cal_stddev dinuc_islands dimob_islands defline2gi);
@EXPORT_OK = qw(cal_dinuc cal_mean cal_median cal_stddev dinuc_islands dimob islands defline2gi);

#use Data::Dumper; ##enable only when trouble shooting

sub cal_dinuc {

#input a fasta or genbank file containing nucleotide sequences of genes
#output an array of hash containing ORF_label and DINUC_bias of each 6-gene cluster
#an array is needed to keep the gene order
#optionally, one can supply and output filename to have the results written to a file

	my $input_fasta = shift @_;
	my $output_file = shift;
	my $fasta_name  = basename($input_fasta);
	##if an output file name is provided, write the results to a file
	##otherwise, just output the results as a hash
	if ($output_file) {
		open( OUTFILE, ">$output_file" )
		  or croak "Can not open output file $output_file.\n";
	}
	my $genome = Bio::SeqIO->new(
		'-file'   => $input_fasta,
		'-format' => 'Fasta'
	  )
	  or croak "no $input_fasta\n";
	my $seqobj;
	my $seqleng;     # sequence length
	my $seqshift;    # seqobj with 1 nucleotide shift (N2 to Nn)
	my $seqtrunc
	  ; # seqobj with last nucleotide truncated; used only for seq with odd number nucleotides
	my $seqrev;         # reverse complement of the seqobj
	my $seqrevshift;    # seqrev with 1 nuc shift
	my $seqrevtrunc;
	my $seq_word;       #seqword object
	my $seq_wordshift;
	my $seq_wordrev;
	my $seq_wordrevshift;
	my $monohash_ref;
	my $monorevhash_ref;
	my $monoshifthash_ref;
	my $monorevshifthash_ref;
	my $dihash_ref;
	my $dishifthash_ref;
	my $direvhash_ref;
	my $direvshifthash_ref;
	my %genomemono;     # hash of cumulative monomer counts of all ORFs
	my %genomedi;       # hash of cumulative dinuc counts of all ORFs
	my @allorfsdi
	  ;    # array of hash of dinuc counts of each of the ORFs in a genome
	my @allorfsmono
	  ;    #array of hash of monomer counts of each of the ORFs in a genome

	## since SeqWords only returns keys of dinucleotides that are present
	## in the sequence, it is necessary to fill the hash with all
	## dinucleotide keys, including the ones not present and give them
	## a value of 0, the whole set of keys is established below as an array.

	my @dinuc_keys = qw/ AA AC AG AT CC CA CG CT GG GA GC GT TT TA TC TG /;
	my $dinuc_key;
	my @ORF_ids;    #array of IDs
	while ( $seqobj = $genome->next_seq() ) {
		push @ORF_ids, $seqobj->id;
		$seqleng = $seqobj->length();
		##if even number of dinucleotides in the ORF
		if ( $seqleng % 2 == 0 ) {
			$seqshift = Bio::PrimarySeq->new(
				-seq      => ( $seqobj->subseq( 2, $seqleng - 1 ) ),
				-alphabet => 'dna',
				-id       => ( $seqobj->display_id )
			);

			#print SEQOUTFILE "The shifted sequence is ".$seqshift->seq()."\n";
			$seq_word      = Bio::Tools::SeqWords2->new( -seq => $seqobj );
			$seq_wordshift = Bio::Tools::SeqWords2->new( -seq => $seqshift );
			$seqrev        = $seqobj->revcom;

#print SEQOUTFILE "The revcom of ".$seqobj->display_id()."is ".$seqrev->seq()."\n";
			$seqrevshift = Bio::PrimarySeq->new(
				-seq      => ( $seqrev->subseq( 2, $seqleng - 1 ) ),
				-alphabet => 'dna',
				-id       => ( $seqrev->display_id )
			);

	  #print SEQOUTFILE "The shifted rev sequence is ".$seqrevshift->seq()."\n";
			$seq_wordrev = Bio::Tools::SeqWords2->new( -seq => $seqrev );
			$seq_wordrevshift =
			  Bio::Tools::SeqWords2->new( -seq => $seqrevshift );
			$monohash_ref         = $seq_word->count_words('1');
			$monorevhash_ref      = $seq_wordrev->count_words('1');
			$monoshifthash_ref    = $seq_wordshift->count_words('1');
			$monorevshifthash_ref = $seq_wordrevshift->count_words('1');
			$dihash_ref           = $seq_word->count_words('2');
			$dishifthash_ref      = $seq_wordshift->count_words('2');
			$direvhash_ref        = $seq_wordrev->count_words('2');
			$direvshifthash_ref   = $seq_wordrevshift->count_words('2');
		}
		else {
			$seqtrunc = Bio::PrimarySeq->new(
				-seq      => ( $seqobj->subseq( 1, $seqleng - 1 ) ),
				-alphabet => 'dna',
				-id       => ( $seqobj->display_id )
			);
			$seqshift = Bio::PrimarySeq->new(
				-seq      => ( $seqobj->subseq( 2, $seqleng ) ),
				-alphabet => 'dna',
				-id       => ( $seqobj->display_id )
			);

			#print SEQOUTFILE "The shifted sequence is ".$seqshift->seq()."\n";
			$seq_word      = Bio::Tools::SeqWords2->new( -seq => $seqtrunc );
			$seq_wordshift = Bio::Tools::SeqWords2->new( -seq => $seqshift );
			$seqrev        = $seqobj->revcom;

#print SEQOUTFILE "The revcom of ".$seqobj->display_id()."is ".$seqrev->seq()."\n";
			$seqrevtrunc = Bio::PrimarySeq->new(
				-seq      => ( $seqrev->subseq( 1, $seqleng - 1 ) ),
				-alphabet => 'dna',
				-id       => ( $seqrev->display_id )
			);
			$seqrevshift = Bio::PrimarySeq->new(
				-seq      => ( $seqrev->subseq( 2, $seqleng ) ),
				-alphabet => 'dna',
				-id       => ( $seqrev->display_id )
			);

	  #print SEQOUTFILE "The shifted rev sequence is ".$seqrevshift->seq()."\n";
			$seq_wordrev = Bio::Tools::SeqWords2->new( -seq => $seqrevtrunc );
			$seq_wordrevshift =
			  Bio::Tools::SeqWords2->new( -seq => $seqrevshift );
			$monohash_ref         = $seq_word->count_words('1');
			$monorevhash_ref      = $seq_wordrev->count_words('1');
			$monoshifthash_ref    = $seq_wordshift->count_words('1');
			$monorevshifthash_ref = $seq_wordrevshift->count_words('1');
			$dihash_ref           = $seq_word->count_words('2');
			$dishifthash_ref      = $seq_wordshift->count_words('2');
			$direvhash_ref        = $seq_wordrev->count_words('2');
			$direvshifthash_ref   = $seq_wordrevshift->count_words('2');
		}
		my %monohash         = %$monohash_ref;
		my %monorevhash      = %$monorevhash_ref;
		my %monoshifthash    = %$monoshifthash_ref;
		my %monorevshifthash = %$monorevshifthash_ref;
		my %dihash           = %$dihash_ref;
		my %dishift          = %$dishifthash_ref;
		my %direvhash        = %$direvhash_ref;
		my %direvshift       = %$direvshifthash_ref;
		my $mono;
		my $di;

		# combine mononculeotide counts from all 4 strands
		foreach $mono ( keys %monorevhash ) {
			$monohash{$mono} += $monorevhash{$mono};
		}
		foreach $mono ( keys %monoshifthash ) {
			$monohash{$mono} += $monoshifthash{$mono};
		}
		foreach $mono ( keys %monorevshifthash ) {
			$monohash{$mono} += $monorevshifthash{$mono};
		}

		# combine dinculeotide counts from all 4 strands
		foreach $di ( keys %direvhash ) {
			$dihash{$di} += $direvhash{$di};
		}
		foreach $di ( keys %dishift ) {
			$dihash{$di} += $dishift{$di};
		}
		foreach $di ( keys %direvshift ) {
			$dihash{$di} += $direvshift{$di};
		}

		# make sure the hash has a full set of dinucleotide keys
		foreach $dinuc_key (@dinuc_keys) {
			if ( exists $dihash{$dinuc_key} ) {
			}
			else {
				$dihash{$dinuc_key} = 0;
			}
		}

		#		  #print out the mononuc counts of each ORF - works
		#		  foreach my $key(sort keys %monohash)
		#		  {
		#		   	  print OUTFILE "$key\t $monohash{$key}\n";
		#		  }
		#		  #print out the dinuc counts of each ORF - works
		#		  foreach my $key(sort keys %dihash){
		#		  	  print OUTFILE "$key\t $dihash{$key}\n";
		#		  }

		# push onto the arrays current ORF's di and mono nucleotide profile
		push @allorfsdi,   \%dihash;
		push @allorfsmono, \%monohash;

		# add mononucleotide counts of single orf to overall genome profile
		while ( ( my $key, my $val ) = each %monohash ) {
			$genomemono{$key} += $val;
		}

		# add dinucleotide counts of single orf to overall genome profile
		while ( ( my $key, my $val ) = each %dihash ) {
			$genomedi{$key} += $val;
		}
	}

	#print Dumper(@allorfsmono);
	#print Dumper(@allorfsdi);
	#print Dumper(%genomemono);
	#print Dumper(%genomedi);

	#calculate the overall genome dinucleotide profile
	my $base1;
	my $base2;
	my $totalnuc;
	my $totaldinuc;
	my %genome_profile;    #pooled dinucleotide relative abundance of all ORFs
	                       #count up total number of mono and di-nucleotides
	while ( ( my $key, my $val ) = each %genomemono ) {
		$totalnuc += $val;
	}
	while ( ( my $key, my $val ) = each %genomedi ) {
		$totaldinuc += $val;
	}

	#print "$totalnuc\n";
	#print "$totaldinuc\n";

	foreach $base1 ( keys %genomemono ) {
		foreach $base2 ( keys %genomemono ) {
			my $dinuc = $base1 . $base2;

			#print "$dinuc\n";
			$genome_profile{$dinuc} =
			  ( $genomedi{$dinuc} / $totaldinuc ) /
			  ( ( $genomemono{$base1} / $totalnuc ) *
				  ( $genomemono{$base2} / $totalnuc ) );
		}
	}

	# print out dinucleotide genome signature - works
	# print Dumper(%genome_profile);

	#calculate the orf dinucleotide profile
	#take 6 ORFs at a time and stop the process at n-6 position
	my $i = @allorfsdi;
	my @allorfdi_prof;
	for ( my $index = 0 ; $index <= $i - 6 ; $index++ ) {
		my %orfdi_prof;
		my %orf_di;
		my %orf_mono;
		for ( my $k = 0 ; $k <= 5 ; $k++ ) {
			my %orf_di2   = %{ $allorfsdi[ $index + $k ] };
			my %orf_mono2 = %{ $allorfsmono[ $index + $k ] };
			while ( ( my $key, my $val ) = each %orf_di2 ) {
				$orf_di{$key} += $val;
			}
			while ( ( my $key, my $val ) = each %orf_mono2 ) {
				$orf_mono{$key} += $val;
			}
		}

		#count up total number of mono and di-nucleotides
		my $orfnuc;
		my $orfdinuc;
		while ( ( my $key, my $val ) = each %orf_mono ) {
			$orfnuc += $val;
		}
		while ( ( my $key, my $val ) = each %orf_di ) {
			$orfdinuc += $val;
		}

		my $nuc1;
		my $nuc2;
		foreach $nuc1 ( keys %orf_mono ) {
			foreach $nuc2 ( keys %orf_mono ) {
				my $dinuc = $nuc1 . $nuc2;

				#print OUTFILE $dinuc; #used to check that concatination works
				$orfdi_prof{$dinuc} =
				  ( $orf_di{$dinuc} / $orfdinuc ) /
				  ( ( $orf_mono{$nuc1} / $orfnuc ) *
					  ( $orf_mono{$nuc2} / $orfnuc ) );
			}
		}
		push @allorfdi_prof, \%orfdi_prof;
	}

	#print out the dinucleotide signature of each ORF - works
	#print Dumper(@allorfdi_prof);

	#calculate the dinucleotide bias of each orf
	my $j = @allorfdi_prof;
	my @biases;
	for ( my $index = 0 ; $index <= $j - 1 ; $index++ ) {
		my %orf_profile = %{ $allorfdi_prof[$index] };
		my $dinuc;
		my $bias = 0;
		foreach $dinuc ( keys %orf_profile ) {
			$bias += abs( $orf_profile{$dinuc} - $genome_profile{$dinuc} );
		}
		$biases[$index] = ( $bias / 16 );
	}

	if ($output_file) {
		print OUTFILE "ORFs_dinucleotide_analysis_for $fasta_name\n";
		my $bias2;
		my $count = 0;
		foreach $bias2 (@biases) {
			$count++;
			$bias2 = $bias2 * 1000;
			printf OUTFILE
			  "ORF%5d-ORF%5d($ORF_ids[$count]-$ORF_ids[$count+5])=%5.2f\n",
			  $count, ( $count + 5 ), $bias2;
		}
		close OUTFILE;
		return;
	}
	else {
		my $bias2;
		my $count = 0;
		my $range = 5;
		my @results;
		my $key_string;
		foreach $bias2 (@biases) {
			$count++;
			$bias2      = $bias2 * 1000;
			$key_string =
			    "ORF" . $count . "("
			  . $ORF_ids[$count] . ")-ORF"
			  . ( $count + $range ) . "("
			  . $ORF_ids[ $count + $range ] . ")";
			push @results, { ORF_label => $key_string, DINUC_bias => $bias2 };
		}
		return \@results;
	}
}

sub dinuc_islands {

#Determine which ORFs are in Dinuc Biased region and return their IDs and bias values
#input a array of a hash containing ORF names (IDs) and their dinuc bias values (from cal_dinuc)
#also input the mean and the standard deviation values (from cal_mean and cal_stddev)

	my $ORFs_dinuc_array = shift @_;
	my $mean             = shift @_;
	my $sd               = shift @_;
	my $cutoff           = shift @_ || 8;    #default is 8
	my @dinuc;    #index of array elements with dinuc bias
	push @dinuc, 0;

	#in this round, keep track all ORFs that have dinuc bias
	#regardless the cutoff, we'll eliminate the fragments smaller
	#than the cutoff later.
	my $i =
	  0;    #for keeping track the index of array elements that have dinuc bias
	foreach my $ORF_dinuc (@$ORFs_dinuc_array) {
		my $orfdbias = $ORF_dinuc->{'DINUC_bias'};
		if ( $orfdbias > $mean + ( $sd * 2 ) ) {
			my @temp;
			my $temp;
			push @temp, $i;
			push @temp, $i + 1;
			push @temp, $i + 2;
			push @temp, $i + 3;
			push @temp, $i + 4;
			push @temp, $i + 5;

			#print STDERR @temp;
			foreach $temp (@temp) {
				if ( $temp > $dinuc[-1] ) {
					push @dinuc, $temp;
				}
			}
		}
		elsif ( $orfdbias > $mean + $sd ) {
			my @temp;
			my $temp;
			push @temp, $i;
			push @temp, $i + 1;
			push @temp, $i + 2;
			foreach $temp (@temp) {
				if ( $temp > $dinuc[-1] ) {
					push @dinuc, $temp;
				}
			}
		}
		$i++;
	}
	shift @dinuc;    #remove the initial zero at the beginning of the array

	#here, we will only keep clusters of dinuc biased ORFs greater than the size cutoff
	my @temp;
	my @dinuc_abovecut;
	my $count = 0;
	my @islands;

	my $island_index = 0;
	foreach my $element (@dinuc) {
		if ($element != $dinuc[-1]){
			# test if numbers in the array are consecutive
			if ( ( $element + 1 ) == $dinuc[ $count + 1 ] ) {
				push @temp, $element;
			}
			else {
				if ( scalar(@temp) >= $cutoff - 1 ) {
					push @temp, $element;
					push( @dinuc_abovecut, @temp );
					@temp = ();
				}
				else {
					@temp = ();
				}
			}
			$count++;
		}
		else {
			if ( scalar(@temp) >= $cutoff - 1 ) {
				push @temp, $element;
				push( @dinuc_abovecut, @temp );
				@temp = ();
			}
			else {
				@temp = ();
			}
		}
	}

	# We will merge regions of dinuc_abovecut that are less than 5 genes apart
	# We do not export the last gene of each dinuc_abovecut as adding the last gene significantly
	# diminishes precision more than it increases recall.
	for ( my $i = 0 ; $i < (scalar(@dinuc_abovecut)-1) ; $i++ ) {
		if ( ( $dinuc_abovecut[$i] + 1) == ( $dinuc_abovecut[ $i + 1 ] ) ) {
			push @{ $islands[$island_index] }, $ORFs_dinuc_array->[ $dinuc_abovecut[$i]];
		}
		elsif ( ( $dinuc_abovecut[$i] + 6) >= ( $dinuc_abovecut[ $i + 1 ] ) ) {
			my $dif = $dinuc_abovecut[ $i + 1 ]-$dinuc_abovecut[$i]-1;
			for (my $j = 0; $j <= $dif; $j++ ) {
				push @{ $islands[$island_index] }, $ORFs_dinuc_array->[ $dinuc_abovecut[$i]+$j];
			}
		}
		#elsif (( $dinuc_abovecut[$i] ) == ( $dinuc_abovecut[ $i - 1 ]+1 ) ) {
			#push @{ $islands[$island_index] }, $ORFs_dinuc_array->[ $dinuc_abovecut[$i] ];
    	#}
		else {
			$island_index++;
		}
	}
	# We do not export the last gene of the last dinuc island for the same reasons as mentioned above
	# push @{ $islands[$island_index] }, $ORFs_dinuc_array->[ $dinuc_abovecut[-1] ];

	return \@islands;
}



sub dimob_islands {

	#take a list of dinuc islands (from sub defline2gi) and
	#a list of mobility genes (from mobgene.pm's output
	# - basically a hash structure with gi numbers as the keys)

	my $dinuc_island_orfs = shift;
	my $mobgenes          = shift;
	#	print Dumper $dinuc_island_orfs;
	my @dimob_island_orfs;
	
	foreach my $island (@$dinuc_island_orfs){
		my $is_dimob= 0; #false
		foreach my $orf (@$island){
			my $orf_ginum = $orf->{'ORF_label'};
			if (exists $mobgenes->{$orf_ginum}){
				$is_dimob = 1; #true
			}
		}
		if ($is_dimob){
			push @dimob_island_orfs, $island;
		}
	}
	return \@dimob_island_orfs;
}

sub defline2gi {

   #take the outpupt of dinuc_islands and convert the ffn def line to gi numbers
   #then return the same data structure with ORF_
   #req input: a ptt file for annotation and output from dinuc_islands
	my $dinucislands = shift @_;
	my $pttfilename  = shift @_;
	my $extended_ids = @_ ? shift : 0;
	my $header_line  = 3;
	my @result_islands;
	my ( $header_arrayref, $pttfh ) =
	  extract_headerandbodyfh( $pttfilename, $header_line );
	my $ptt_table_hashref = table2hash_rowfirst( $header_arrayref, $pttfh, 1 );
	foreach my $island (@$dinucislands) {
		my @result_orfs;
		ORF: foreach my $orf_index (@$island) {
			my $def_label = $orf_index->{'ORF_label'};
			next ORF unless($def_label);
			my ( $orf1, $orf2 ) = split '-ORF', $def_label;
			my $orf_start;
			my $orf_end;
			my $pid;
			my $seq_id;
			my $strand = '.';
#			print "$orf1\n";

			if ( $orf1 =~ /.*\|([^:]*):(\d+)\.\.(\d+)\)/ ) {
				$seq_id = $1;
				$strand = '+';
				$orf_start = $2;
				$orf_end   = $3;
				my $coordinate = "$orf_start..$orf_end";
			 	$pid = $ptt_table_hashref->{$coordinate}->{'PID'};
				if($extended_ids) {
#				    print "Using extended\n";
				    $pid .= "_$orf_start" . '..' . $orf_end;
				}
				unless(defined($pid)){
				    #warn "Could not find pid";
				}
			}elsif ( $orf1 =~ /.*\|([^:]*):c(\d+)\.\.(\d+)\)/ ) {
				$seq_id = $1;
				$strand = '-';
				$orf_start = $2;
				$orf_end   = $3;
				my $coordinate = "$orf_start..$orf_end";
				$pid = $ptt_table_hashref->{$coordinate}->{'PID'};
				if($extended_ids) {
#				    print "Using extended (comp)\n";
				    $pid .= "_$orf_start" . '..' . $orf_end;
				}
				unless(defined($pid)){
				    #warn "Could not find pid";
				}
			    }
			#Morgan Hack: sometimes we don't need look up pid by coordinates
			elsif($orf1 =~ /\((\d+)\)/ ){
                $pid=$1;	
			} else {

			}

			#print "$orf_start and $orf_end\n";
		
			$orf_index->{'ORF_label'} = $pid;
			$orf_index->{'seq_id'} = $seq_id;
			$orf_index->{'strand'} = $strand;
			$orf_index->{'start'}=$orf_start;
			$orf_index->{'end'}=$orf_end;
			push @result_orfs, $orf_index;
		}
		push @result_islands, [@result_orfs];
	}    
	return \@result_islands;
}

sub cal_mean {
	my $input_array = shift;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data($input_array);
	my $mean = $stat->mean();
	return $mean;
}

sub cal_stddev {
	my $input_array = shift;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data($input_array);
	my $stddev = $stat->standard_deviation();
	return $stddev;
}

sub cal_median {
	my $input_array = shift;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data($input_array);
	my $median = $stat->median();
	return $median;
}
