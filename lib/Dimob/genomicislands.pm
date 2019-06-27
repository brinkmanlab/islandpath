=head1 DESCRIPTION

    A repackaging of Will's original Dimob software
    loosely associated subroutines used by Dimob.pm for generating dinuc and dimob islands

=head1 SYNOPSIS

    use Dimob::genomicislands;

	see Dimob.pm

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


package Dimob::genomicislands;

## Added: 25 June 2019 by J.F.Sanchez-Herrero
## Add multicontig. Calculate dinuc bias within contig

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::SeqWords2;
use File::Basename;
use Statistics::Descriptive;
use Getopt::Long;
use Carp;
use Data::Dumper;

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
		open( OUTFILE, ">$output_file" ) or croak "Can not open output file $output_file.\n";
	}
	my $genome = Bio::SeqIO->new(
		'-file'   => $input_fasta,
		'-format' => 'Fasta'
	  ) or croak "no $input_fasta\n";
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
	my @allorfsdi;    # array of hash of dinuc counts of each of the ORFs in a genome
	my @allorfsmono;    #array of hash of monomer counts of each of the ORFs in a genome

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
			$seq_wordrevshift = Bio::Tools::SeqWords2->new( -seq => $seqrevshift );
			$monohash_ref         = $seq_word->count_words('1');
			$monorevhash_ref      = $seq_wordrev->count_words('1');
			$monoshifthash_ref    = $seq_wordshift->count_words('1');
			$monorevshifthash_ref = $seq_wordrevshift->count_words('1');
			$dihash_ref           = $seq_word->count_words('2');
			$dishifthash_ref      = $seq_wordshift->count_words('2');
			$direvhash_ref        = $seq_wordrev->count_words('2');
			$direvshifthash_ref   = $seq_wordrevshift->count_words('2');
		} else {
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
			$seq_wordrevshift = Bio::Tools::SeqWords2->new( -seq => $seqrevshift );
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
			if (!exists $dihash{$dinuc_key} ) {
				$dihash{$dinuc_key} = 0;
		}}

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

	##print "ORFs dinucleotide analysis for $fasta_name\n";
	my $bias2;
	my $count = 0;
	my @results;
	print OUTFILE "##id1,id2,ORF1,contig_ORF1,start_1,ORF2,contig_ORF2,start_2,bias\n";
	foreach $bias2 (@biases) {
		## $count++; ## if here, ORF_ids[0] is missing
		my $tmp_bias2 = $bias2 * 1000;
		my $tmp=$count+6;
		my $orf1 = $ORF_ids[$count];
		my $orf2 = $ORF_ids[$count+5];
		my $seq1; my $seq2; my $start_1; my $start_2;

		## get contig id
		if ( $orf1 =~ /.*\|(.*)\:c{0,1}(\d+)\.\.(\d+).*/ ) { $seq1 = $1; $start_1 = $2;}
		if ( $orf2 =~ /.*\|(.*)\:c{0,1}(\d+)\.\.(\d+).*/ ) { $seq2 = $1; $start_2 = $2; }
		
		## skip if different contig ids
		if ($seq1 ne $seq2) {  
			$count++; ## fix missing $ORF_ids[0]
			next; 
		} 
		my $key_string = $count.",".$tmp.",".$orf1.",".$seq1.",".$start_1.",".$orf2.",".$seq2.",".$start_2.",".$tmp_bias2; ## return csv instead
		print OUTFILE $key_string."\n";
		push (@results, { ORF_label => $key_string, DINUC_bias => $tmp_bias2 });
		$count++; ## fix missing $ORF_ids[0]
	}
	close OUTFILE;
	return \@results;
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

	my $i = 0;    #for keeping track the index of array elements that have dinuc bias
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
		} elsif ( $orfdbias > $mean + $sd ) {
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
	my $previous_contig = "";
	my $island_index = 0;
	for ( my $i = 0 ; $i < (scalar(@dinuc_abovecut)-1) ; $i++ ) {
		
		## control if it is ok
		my $def_label;
		if (exists $ORFs_dinuc_array->[$dinuc_abovecut[$i]]->{'ORF_label'}) {
			$def_label = $ORFs_dinuc_array->[$dinuc_abovecut[$i]]->{'ORF_label'};
		} else { 
			next;
		}

		## control for different contig sequences
		my ($id1, $id2, $orf1, $contig_ORF1, $orf2, $contig_ORF2, $bias) = split ',', $def_label;
		if (!$previous_contig) { $previous_contig = $contig_ORF1; }
		if ($previous_contig eq $contig_ORF1) {
		} else { 
			$previous_contig = $contig_ORF1;
			$island_index++;
		}
		
		## fill information
		if ( ( $dinuc_abovecut[$i] + 1) == ( $dinuc_abovecut[ $i + 1 ] ) ) {
			push @{ $islands[$island_index] }, $ORFs_dinuc_array->[ $dinuc_abovecut[$i] ];
		
		# We will merge regions of dinuc_abovecut that are less than 5 genes apart
		# and within the same contig region
		} elsif ( ( $dinuc_abovecut[$i] + 6) >= ( $dinuc_abovecut[ $i + 1 ] ) ) {

			my $dif = $dinuc_abovecut[ $i + 1 ]-$dinuc_abovecut[$i]-1;
			for (my $j = 0; $j <= $dif; $j++ ) {
		
				my $def_label_here = $ORFs_dinuc_array->[ $dinuc_abovecut[$i]+$j]->{'ORF_label'};
				my ($id1, $id2, $orf1, $contig_ORF1, $orf2, $contig_ORF2, $bias) = split ',', $def_label_here;
				if ($previous_contig ne $contig_ORF1) {
					#print "Contig jump! $previous_contig & $contig_ORF1\n";
				} else {
					## only add info if in the same contig
					push @{ $islands[$island_index] }, $ORFs_dinuc_array->[ $dinuc_abovecut[$i]+$j];
				}
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
	my @dimob_island_orfs;
	
	###############
	## dinuc_island_orfs
	###############
	## 'ORF_label' => '16804840_2888639..2889001',
    ## 'DINUC_bias' => '133.287399842604',
    ## 'seq' => 'NC_003210',
    ## 'start' => '2888639',
    ## 'end' => '2889001'

	###############
	## mob list
	###############
	# 'ref|NP_463859.1|gi|16802374|seq:NC_003210:356567..356872' => 1,
	# 'ref|YP_008475635.1|gi|533240825|seq:NC_003210:c1132787..1132987' => 1,
	# '16804339_NC_003210' => 1,
	
	foreach my $island (@$dinuc_island_orfs){
		my $is_dimob= 0; #false
		foreach my $orf (@$island){
			my $orf_ginum = $orf->{'ORF_label'};
			#print "ORF ginum: ".$orf_ginum."\n";
			if (exists $mobgenes->{$orf_ginum}){
				$is_dimob = 1; #true
			} else {
				my $orf_name = $orf->{'orf1'};
				#print "ORF name: ".$orf_name."\n";
				if (exists $mobgenes->{$orf_name}){
					$is_dimob = 1; #true
				} else {
					#print "Dont know...\n";
				}
			}
		}
	    #any dinuc islands containing >=1 mobility gene are classified as dimob islands
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
	my $header_line  = 1; ## header is in line 1
	my @result_islands;
	my ( $header_arrayref, $pttfh ) = extract_headerandbodyfh( $pttfilename, $header_line );
	my $ptt_table_hashref = table2hash_rowfirst( $header_arrayref, $pttfh, 2); ## in column1 is the sequence
	
	foreach my $island (@$dinucislands) {
		my @result_orfs;
		ORF: foreach my $orf_index (@$island) {
			my $def_label = $orf_index->{'ORF_label'};
			next ORF unless($def_label);
			##my ( $orf1, $orf2 ) = split '-ORF', $def_label;
			my ($id1, $id2, $orf1, $contig_ORF1, $start_1, $orf2, $contig_ORF2, $start_2, $bias) = split ',', $def_label;
			my ($orf_start, $orf_end, $pid, $coordinate, $contig, $annotation);

			## just in case
			if ($contig_ORF1 ne $contig_ORF2) { print "Different contigs...\n"; next; }
			if ( $orf1 =~ /.*\|(.*)\:c{0,1}(\d+)\.\.(\d+).*/ ) {
				$contig = $1;
				$orf_start = $2;
				$orf_end   = $3;
				$coordinate = "$orf_start..$orf_end"."_".$contig;
			 	$pid = $ptt_table_hashref->{$coordinate}->{'PID'};
			}

			#Morgan Hack: sometimes we don't need look up pid by coordinates
			elsif($orf1 =~ /\((\d+)\)/ ){
                $pid=$1;	
			} else {
			
			}

			#print "$orf_start and $orf_end\n";
			unless(defined($coordinate)){
			    #warn "Could not find pid";
			}		
			$orf_index->{'ORF_label'} = $coordinate;
			$orf_index->{'start'}=$orf_start;
			$orf_index->{'end'}=$orf_end;
			$orf_index->{'seq'}=$contig;
			$orf_index->{'orf1'}=$orf1;
			$orf_index->{'PID'}= $pid;
			$orf_index->{'strand'}=$ptt_table_hashref->{$coordinate}->{'Strand'};
		 	$annotation = $ptt_table_hashref->{$coordinate}->{'Synonym'}.";".
				 	$ptt_table_hashref->{$coordinate}->{'Gene'}.";".
				 	$ptt_table_hashref->{$coordinate}->{'Product'};

			$orf_index->{'annot'}=$annotation;
			
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
