=head1 NAME

    GenomeUtils

=head1 DESCRIPTION

    Object to load, convert and store custom genomes into
    the internal formats needed by IslandPath

=head1 SYNOPSIS

    use GenomeUtils;

    my $genome_obj = GenomeUtils->new(
                            { workdir => '/tmp/dir/'});

    # $genome_name optional name, will be custom_genome otherwise
    $genome_obj->read_and_convert($filename, $genome_name);

    my $success = $genome_obj->insert_custom_genome();

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

package GenomeUtils;

use strict;
use warnings;
use Moose;
use Log::Log4perl qw(get_logger :nowarn);
use File::Copy;
use File::Basename;
#use Array::Utils qw(:all);
use Data::Dumper;
use File::Temp qw/ :mktemp /;
use File::Spec;
use Dimob::Config;

use Bio::Perl;
use Bio::SeqIO;
use Bio::Seq;

my $cfg; my $logger;

sub BUILD {
    my $self = shift;
    my $args = shift;

    $cfg = Dimob::Config->config;

    $logger = Log::Log4perl->get_logger;

#    die "Error, work dir not specified: $args->{workdir}"
#				unless( -d $args->{workdir} );
#    $self->{workdir} = $args->{workdir};

}


sub read_and_check {
    my $self = shift;
    my $genome_obj = shift;
    my $filename = shift;

    # Separate extension from filename
    $filename =~ s/\/\//\//g;
    # A bit hacky because new NCBI paths have . in the names now
    # Cut the file name off
    my ($volume,$directories,$basefile) =
        File::Spec->splitpath( $filename );
    # Check if it has an extension
    my ( $file, $extension ) = $basefile =~ /(.+)\.(\w+)/;
    # Put the base filename back on to the path
    $file = File::Spec->catpath(undef, $directories, $file);

    unless($extension) {
			$logger->info("Didn't receive file type for $filename");

			# Check if an embl or genbank file exists for this genome
			if(-f $filename . '.gbk' &&
				 -s $filename . '.gbk') {
					$logger->info("We seem to have a genbank file, preferred format");
					$extension = 'gbk';
					$file = $filename;
					$filename .= '.gbk';
			} elsif(-f $filename . '.embl' &&
				-s $filename . '.embl') {
					$logger->info("We seem to have a embl file");
					$extension = 'embl';
					$file = $filename;
					$filename .= '.embl';
			} else {
					$logger->logdie("Can't find file format for $filename, this is very bad");
			}
    }

    $logger->debug("From filename $filename got $file, $extension");

    $self->{base_filename} = $file;
   
    my $in;

    if ( $extension =~ /embl/ ) {
	
			$logger->trace("Reading embl format file");
			$in = Bio::SeqIO->new(
					-file   => $filename,
					-format => 'EMBL'
					);
			$logger->info("The genome sequence in $filename has been read.");
    } elsif ( ($extension =~ /gbk/) || ($extension =~ /gb/) || ($extension =~ /gbf/) || ($extension =~ /gbff/) ) {
			$logger->trace("Reading genbank format file");

			# Special case, our general purpose code likes .gbk...
			if($extension !~ /gbk/) {
					move($filename, "$file.gbk");
					$filename = "$file.gbk";
			}

			$in = Bio::SeqIO->new(
					-file   => $filename,
					-format => 'GENBANK'
					);
			$logger->info("The genome sequence in $filename has been read.");
    } else {
			$logger->logdie("Can't figure out if file is genbank (.gbk) or embl (.embl) [FILEFORMATERROR]");
    }

    # Count the contigs to see if this is
    # an incomplete genome
    my $contigs = 0;

    # Did we find any CDS records?
    my $found_cds = 0;

    my $full_seq_recs;

    SEQ: while ( my $seq = $in->next_seq() ) {
			$contigs += 1;
			$logger->trace("Checking contig " . $seq->accession_number);

			#Only keep those features coding for proteins
			my @cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;

			# We found a cds record in at least one contig
			$found_cds = 1 if(@cds);

			# See if we have a full sequence in the genbank/embl file
			if($seq->seq()) {
					$logger->trace("Found sequence in the genbank/embl file");
					# All good, next.
					next;
			} elsif($full_seq_recs || $self->load_fna($file, \$full_seq_recs)) {
					$logger->trace("Sequence missing, keys available: " . join(',', keys %{$full_seq_recs}));
					$genome_obj->genome_status('MISSINGSEQ');
		#	    print Dumper $full_seq_recs;
					# Do we have sequence information loaded from
					# an fna file?

					# In case the primary accession is not the one used in the fna file...
					foreach my $acc ($seq->accession_number, $seq->get_secondary_accessions) {
				$logger->trace("Looking up seq for $acc");

				# Because the FASTA header could have multiple
				# identifiers, we have to loop through all keys and try
				# to find one with one of our accessions in it... ugh.
				for my $ids (keys %{$full_seq_recs}) {
						if($ids =~ /$acc/) {
					$logger->info("We found identifier $acc in $ids");
					next SEQ;
						}
				}
		#		next SEQ if($full_seq_recs->{$acc});
					}

					$logger->logdie("Error, no sequence for contig [" . $seq->accession_number . '], fna file was found [NOSEQFNA]');
			} else {
					$logger->logdie("Error, no sequence for contig [" . $seq->accession_number . '], fna file was not found [NOSEQNOFNA]');
			}

    }

    unless($found_cds) {
			$logger->logdie("Error, no cds records found for file $filename [NOCDSRECORDS]");
    }

    # Else return the number of contigs found
    $logger->info("We have $contigs contigs");
    return $contigs;

}


# Try to find an fna file for the genome and load
# the records in to a hash

sub load_fna {
    my $self = shift;
    my $basefile = shift;
    my $seq_recs_ref = shift;

    my $fna_file = $basefile . '.fna';

    $logger->debug("Seeing if fna file $fna_file exists, and loading");

    unless(-r $fna_file && -s $fna_file) {
			$logger->warn("fna file $fna_file not found!");
			return 0;
		}

		my $in = Bio::SeqIO->new(
			-file => $fna_file,
			-format => 'FASTA'
			) or $logger->logdie("Error, can't open fna $fna_file using bioperl: $!");

		while(my $seq = $in->next_seq()) {
			return 0 unless($seq->id && $seq->seq());

			my $trimmed_id = $seq->id;
			$trimmed_id =~ s/\.(\d+)$//;

			$logger->trace("Saving sequence for $trimmed_id (removed version numnber $1)");

			$$seq_recs_ref->{$trimmed_id} = $seq;
		#	$$seq_recs_ref->{$trimmed_id} = $seq->seq();
    }

    $logger->trace("Sequence keys available: " . join(',', keys %{$$seq_recs_ref}));

    return 1;
}

# Read in a genbank or Embl file and convert it to
# the other needed formats
#
# Puts the resulting files in the same directory
# as the input genome
#
# Most code recycled from run_custom_islandviewer.pl
# function gbk_or_embl_to_other_formats

sub read_and_convert {
    my $self = shift;
    my $filename = shift;
	my $genome_name = (@_ ? shift : 'custom_genome');

    #separate extension from filename
    $filename =~ s/\/\//\//g;
    my ( $file, $extension ) = $filename =~ /(.+)\.(\w+)/;

	unless($extension) {
		$logger->info("Didn't receive file type for $filename");

		# Check if an embl or genbank file exists for this genome
		if(-f $filename . '.gbk' &&
			-s $filename . '.gbk') {
			$logger->info("We seem to have a genbank file, preferred format");
			$extension = 'gbk';
			$file = $filename;
			$filename .= '.gbk';
		} elsif(-f $filename . '.embl' &&
			-s $filename . '.embl') {
			$logger->info("We seem to have a embl file");
			$extension = 'embl';
			$file = $filename;
			$filename .= '.embl';
		} else {
			$logger->logdie("Can't find file format for $filename, this is very bad");
		}
	}

	$logger->debug("From filename $filename got $file, $extension");

	# We're going to check what files we have, then only
    # generate the ones we need.  Because this code is
    # so nicely compact and to avoid duplication, we're
    # just going to trick the code so if a file exists
    # it gets written to /dev/null
    $self->{base_filename} = $file;
    my $formats = $self->parse_formats($self->find_file_types());

    $logger->trace("Found formats: " . join(',', keys %{$formats}));

    my $in;

    if ( $extension =~ /embl/ ) {

		$in = Bio::SeqIO->new(
				-file   => $filename,
				-format => 'EMBL'
				);
		$logger->info("The genome sequence in $filename has been read.");
		$logger->info("Now generating required file formats (faa, ffn and ptt)");
	} elsif ( ($extension =~ /gbk/) || ($extension =~ /gb/) || ($extension =~ /gbff/)) {
		# Special case, our general purpose code likes .gbk...
		if( ($extension =~ /gb/) || ($extension =~ /gbff/)) {
				move($filename, "$file.gbk");
				$filename = "$file.gbk";
		}

		$in = Bio::SeqIO->new(
				-file   => $filename,
				-format => 'GENBANK'
				);
		$logger->info("The genome sequence in $filename has been read.");
		$logger->info("Now generating required file formats (faa, ffn and ptt)");

	} else {
		$logger->logdie("Can't figure out if file is genbank (.gbk) or embl (.embl)");
	}

	my $outfile = ($formats->{faa} ? '/dev/null' : $file . '.faa');
	my $faa_out = Bio::SeqIO->new(
		-file   => ">" . $outfile,
		-format => 'FASTA'
	);
	$outfile = ($formats->{ffn} ? '/dev/null' : $file . '.ffn');
	my $ffn_out = Bio::SeqIO->new(
		-file   => ">" . $outfile,
		-format => 'FASTA'
	);
	#$outfile = ($formats->{fna} ? '/dev/null' : $file . '.fna');
	#my $fna_out = Bio::SeqIO->new(
	#		-file   => ">" . $outfile,
	#		-format => 'FASTA'
	#		);
	$outfile = ($formats->{ptt} ? '/dev/null' : $file . '.ptt');
	open( my $PTT_OUT, '>', $outfile );

	while ( my $seq = $in->next_seq() ) {
		# open all files

		my $total_length = $seq->length();
		my $seq_id = $seq->object_id();

#		my $total_seq    = $seq->seq();

		#Create fna file
		#$success = $fna_out->write_seq($seq);
		#if ($success == 0) {
		#		$logger->error(".fna file is not generated successfully.");
		#}

		#Only keep those features coding for proteins
		my @cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;

		#Remove any pseudogenes
		my @tmp_cds;
		foreach (@cds) {
				unless ( $_->has_tag('pseudo') ) {
			push( @tmp_cds, $_ );
				}
		}
		@cds = @tmp_cds;

		my $num_proteins = scalar(@cds);

		#Create header for ptt file
		print $PTT_OUT $seq->description, " - 1..", $seq->length, "\n";
		print $PTT_OUT $num_proteins, " proteins\n";
		print $PTT_OUT join(
				"\t", qw(Location Strand Length PID Gene Synonym Code COG
					Product)
				),
				"\n";

		my $count = 0;

		#Step through each protein
		foreach my $feat (@cds) {
			$count++;

			#Get the general features
			my $start  = $feat->start;
			my $end    = $feat->end;
			my $strand = $feat->strand;
			my $length = $feat->length;

			if ($length <= 2) {
				throw Bio::Root::Exception("Something's wrong with one of the protein sequences! CDS info: start=$start end=$end strand=$strand");
			}

			#Get more features associated with gene (not all of these will necessarily exist)
			my ( $product, $protein_accnum, $gene_name, $locus_tag ) = ( '', '', '', '' );
				($product) = $feat->get_tag_values('product')
			if $feat->has_tag('product');
				($protein_accnum) = $feat->get_tag_values('protein_id')
			if $feat->has_tag('protein_id');
				($gene_name) = $feat->get_tag_values('gene')
			if $feat->has_tag('gene');
				($locus_tag) = $feat->get_tag_values('locus_tag')
			if $feat->has_tag('locus_tag');

			my $gi = $count;
			$gi = $1 if tag( $feat, 'db_xref' ) =~ m/\bGI:(\d+)\b/;

			my $ref_accnum = "UN_" . sprintf("%06d", $count) . ".0";
			if($feat->has_tag('protein_id')) {
				$ref_accnum = tag( $feat, 'protein_id' );
			}

			my $strand_expand  = $strand >= 0 ? '+' : '-';
			my $strand_expand2 = $strand >= 0 ? ''  : 'c';
			my $desc = "$seq_id\:$strand_expand2" . "$start..$end";

			$desc = "ref\|$ref_accnum\|gi\|$gi\|" . $desc;

			#Create the ffn seq
			my $ffn_seq = $seq->trunc( $start, $end );
			if ( $strand == -1 ) {
				$ffn_seq = $ffn_seq->revcom;
			}
			$ffn_seq->id($desc);
			$ffn_seq->desc($product);
			my $success = 0;
			$success = $ffn_out->write_seq($ffn_seq);
			if ($success == 0) {
				$logger->error(".ffn file is not generated successfully.");
			}

			#Create the faa seq
			my $faa_seq;
			if ( $feat->has_tag('translation') ) {
				my ($translation) = $feat->get_tag_values('translation');
				$faa_seq = new Bio::Seq( -seq => $translation ) or throw Bio::Root::Exception("Cannot read protein sequence: $!");
			} else {
				$faa_seq =
				$ffn_seq->translate( -codontable_id => 11, -complete => 1 );
			}

			$faa_seq->id($desc);
			$faa_seq->desc($product);
			$success = $faa_out->write_seq($faa_seq);
			if ($success == 0) {
				$logger->error(".faa file is not generated successfully.");
			}

			#Print out ptt line
			my $cog = '-';
			$cog = $1 if tag( $feat, 'product' ) =~ m/^(COG\S+)/;
			my @col = (
			$start . '..' . $end,
			$strand_expand,
			( $length / 3 ) - 1,
			$gi,
			tag( $feat, 'gene' ),
			tag( $feat, 'locus_tag' ),
			'-',
			$cog,
			tag( $feat, 'product' ),
			);
				print $PTT_OUT join( "\t", @col ), "\n";

			}    #end of foreach

		# Save the details of the file we just loaded
		$self->{name} = $genome_name;
		$self->{num_proteins} = $num_proteins;
		$self->{total_length} = $total_length;
		$self->{base_filename} = $file;
		$self->{ext} = $extension;
		$self->{orig_filename} = $filename;
		$self->{formats} = $self->parse_formats($self->find_file_types());
	#	$self->{type} = 'custom';
		$self->{genome_read} = 1;

		close($PTT_OUT);
  }    #end of while
    
}    #end of read_and_convert

# REMOVE THIS FUNCTIONS?
# We really need to go back to square one and manage the files and
# formats better. Or improve the analysis pieces so they share
# formats... ugh. But for now, a tool to convert a genbank to embl
# and vise versa. Mainly needed so Sigi can deal with MicrobeDB,
# since Sigi needs Embl and we no longer generate it in MicrobeDB v2

sub convert_file {
    my $self = shift;
    my $filename = shift;
    my $outfile = shift;

    #separate extension from filename
    $filename =~ s/\/\//\//g;
    my ( $file, $extension ) = $filename =~ /(.+)\.(\w+)/;

    $outfile =~ s/\/\//\//g;
    my ( $outbase, $outextension ) = $outfile =~ /(.+)\.(\w+)/;

    $logger->debug("From filename $filename got $file, $extension");
    $logger->debug("From outfile $outfile got $outbase, $outextension");

    my $in;

    if ( $extension =~ /embl/ ) {

			$in = Bio::SeqIO->new(
					-file   => $filename,
					-format => 'EMBL'
					);
			$logger->info("The genome sequence in $filename has been read.");
    } elsif ( ($extension =~ /gbk/) || ($extension =~ /gb/) || ($extension =~ /gbff/) ) {

			$in = Bio::SeqIO->new(
					-file   => $filename,
					-format => 'GENBANK'
					);
			$logger->info("The genome sequence in $filename has been read.");
    } else {
			$logger->logdie("Can't figure out if file is genbank (.gbk) or embl (.embl)");
    }

    my $out;    
    if ( $outextension =~ /embl/ ) {

        $out = Bio::SeqIO->new(
            -file   => ">" . $outfile,
            -format => 'EMBL'
            );

    } elsif ( ($extension =~ /gbk/) || ($extension =~ /gb/) || ($extension =~ /gbff/) ) {
        $out = Bio::SeqIO->new(
            -file   => ">" . $outfile,
            -format => 'GENBANK'
            );
    } else {
        $logger->logdie("Can't figure out if file is genbank (.gbk) or embl (.embl)");
    }

    while ( my $seq = $in->next_seq() ) {

	#Create gbk or embl file
	$out->write_seq($seq);

    }

}

sub genome_stats {
    my $self = shift;
    my $base_filename = shift;

    my $in;
    # Check which file type exists
    if(-r $base_filename . '.gbk' &&
       -s $base_filename . '.gbk') {
	$logger->info("Scanning genbank file $base_filename for genome stats");

	$in = Bio::SeqIO->new(
	    -file   => $base_filename . '.gbk',
	    -format => 'Genbank'
	    );
    } elsif(-r $base_filename . '.embl' &&
	    -s $base_filename . '.embl') {
	$logger->info("Scanning embl file $base_filename for genome stats");


	$in = Bio::SeqIO->new(
	    -file   => $base_filename . '.embl',
	    -format => 'EMBL'
	    );
    } else {
	$logger->error("Error, no genbank or embl file for $base_filename exists to scan stats from");
	return 0;
    }

    # If we find more than one contig this is a problem, this should never be called
    # on an incomplete or non-assembled genome
    my $contig_count = 1;

    my $stats;

    while ( my $seq = $in->next_seq() ) {
	if($contig_count > 1) {
	    $logger->error("Error, we found a second contig in the genome from $base_filename");
	}

	#Only keep those features coding for proteins
	my @cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;

	#Remove any pseudogenes
	my @tmp_cds;
	foreach (@cds) {
	    unless ( $_->has_tag('pseudo') ) {
		push( @tmp_cds, $_ );
	    }
	}
	@cds = @tmp_cds;

	my $num_proteins = scalar(@cds);
	my $seq_length = $seq->length();

	$stats = { cds_num => $num_proteins,
		   rep_size => $seq_length
	};

        $stats->{name} = $seq->desc if($seq->desc);

	$contig_count += 1;
    }

    return $stats;
}

sub regenerate_files {
    my $self = shift;

    unless($self->{base_filename}) {
	$logger->error("Error, we can't regenerate the files unless we have a base filename");
	return 0;
    }

    if($self->{formats}->{gbk}) {
	$logger->trace("Regenerating based on genbank format");
	$self->read_and_convert($self->{base_filename} . '.gbk', $self->{name});
    } elsif($self->{formats}->{embl}) {
	$logger->trace("Regenerating based on Embl format");
	$self->read_and_convert($self->{base_filename} . '.embl', $self->{name});
    } else {
	$logger->error("Error, we don't have either genbank or embl, can't generate needed files");
	return 0;
    }

    if($cfg->{expected_exts} eq $self->find_file_types()) {
	# The regeneration was successful!
	return 1;
    } else {
	$logger->error("Error, we didn't regenerate all the files we expected to, failed, only have: " . $self->find_file_types());
	return 0;
    }

}

# Validate we have all the needed types

sub validate_types {
    my $self = shift;
    my $genome_obj = shift;

    # First are the file types the genome object thinks we
    # have correct
    my @found_types = $self->find_file_types($genome_obj->filename(), 1);

    my @formats = sort @{$genome_obj->formats()};
    if(array_diff(@formats, @found_types)) {
	$logger->warn("Genome object and file system have different sets of formats, object has [" . join(',', @formats) . '] on disk, [' . join(',', @found_types) . ']');
	$genome_obj->formats(@found_types);
	$genome_obj->update_genome();
    }

    # Next does the type of formats match what we need...
    $logger->trace("Checking formats we have: " . join(',', @formats));
    if(! $self->correct_formats(\@formats) ) {
			if('.gbk' ~~ $genome_obj->formats()) {
					$logger->info("Regenerating formats based off genbank format");
					$self->read_and_convert($genome_obj->filename() . '.gbk', $genome_obj->name());
			} elsif('.embl' ~~ $genome_obj->formats()) {
					$logger->info("Regenerating formats based off embl format");
					$self->read_and_convert($genome_obj->filename() . '.embl', $genome_obj->name());

			} else {
					$logger->error("Error, neither genbank or embl file found for " . $genome_obj->filename());
					return 0;
			}

	# We've updated formats, so we need to update the genome object with 
	# the new formats
	@formats = $self->find_file_types($genome_obj->filename(), 1);
	$logger->trace("Rechecking, now we have: " . join(',', @formats));
	$genome_obj->formats( \@formats );
	$genome_obj->update_genome();

    } else{
			# Nothing to do, we have the correct formats
			return 1;
    }

    # Do we have the correct fromats now?
    @formats = @{$genome_obj->formats()};
    if(! $self->correct_formats(\@formats )) {
			$logger->error("We still don't have all the formats we need, fail! Have: [" . join(',', @{$genome_obj->formats}) . '] Want: [' . $cfg->{expected_exts} . ']');
			return 0;
    }

    # All good, moving on...
    return 1;

}

sub correct_formats {
    my $self = shift;
    my $formats = shift;

    my @formats = sort @{$formats};
    my @expected_formats = sort(split(' ', $cfg->{expected_exts}));
    $logger->trace("Checking formats, have [" . join(',', @formats) . '] need [' . join(',', @expected_formats) . ']');
    if(array_diff(@formats, @expected_formats ) ) {
	$logger->warn("We don't have all the needed formats, have [" . join(',', @formats) . '] need [' . $cfg->{expected_exts} . ']');
	return 0;
    }

    return 1;
}

sub find_file_types {
    my $self = shift;
    my $base_filename = shift;
    my $return_array = shift;

    unless($base_filename) {
	$logger->trace("No base filename given in args, trying to use object default: " . $self->{base_filename});
	$base_filename = $self->{base_filename};
    }

    unless($base_filename) {
	$logger->error("Error, you must specify a base filename or a genome must be read before you can test the file types");
	return '';
    }

    # Fetch and parse the formats we expect to find...
    my @expected;
    foreach (split /\s+/, $cfg->{expected_exts}) { 
			$_ =~ s/^\.//;
			push @expected, $_;
    }

#    my $expected_formats = $self->parse_formats($cfg->{expected_exts});

    my @formats;
    foreach my $ext (@expected) {
	# For each format we expect to find, does the file exist?
	# And is non-zero
	if(-f "$base_filename.$ext" &&
	   -s "$base_filename.$ext") {
	    push @formats, ".$ext";
	}
    }

    # If we've been asked to return it as an array rather 
    # than a string...
    if($return_array) {
	return sort @formats;
    }

    return join ' ', @formats;
}

# Fetch and return the nucleotide sequence for
# a genome object, or the sub-sequence if requested

sub fetch_nuc_seq {
    my $self = shift;
    my $genome_obj = shift;
    my $strand = @_ ? shift : 1;
    my $start = @_ ? shift : undef;
    my $end = @_ ? shift : undef;

    my $filename = $genome_obj->filename() . '.fna';
    $logger->trace("Reading sequences from $filename");

    unless(-f $filename) {
	$logger->logdie("Error, file $filename doesn't exist");
    }

    # Grab the fna file via bioperl
    my $in = new Bio::SeqIO(-file => $filename);

    my $seqobj = $in->next_seq();

    # If we've been given a start and end, then we grab the
    # subsequence, otherwise just return the entire sequence
    if(defined($start) && defined($end)) {
	return (($strand eq '-1' || $strand eq '-') ?
		reverse_complement_as_string($seqobj->subseq($start, $end)) :
		$seqobj->subseq($start, $end));

    } else {
	return (($strand eq '-1' || $strand eq '-') ?
		reverse_complement_as_string($seqobj->seq()) :
		$seqobj->seq());

    }

}

sub fetch_protein_seq {
    my $self = shift;
    my $genome_obj = shift;
    my $strand = @_ ? shift : 1;
    my $start = @_ ? shift : undef;
    my $end = @_ ? shift : undef;

    return Bio::Seq->new(-seq => $self->fetch_nuc_seq($genome_obj, $strand, $start, $end),
		  -alphabet => 'dna')
	->translate(-codontable_id => 11, -complete => 1)
	->seq();

}

# Given a genome_obj and a list of identifiers,
# pull those sequences out of the fasta file and put
# them in a new temp file

sub make_sub_fasta {
    my $self = shift;
    my $genome_obj = shift;
    my $outfile = shift;
    my @accessions = @_;

    my $filename = $genome_obj->filename() . '.faa';
    $logger->trace("Reading sequences from $filename");
    
    unless(-f $filename) {
	$logger->logdie("Error, file $filename doesn't exist");
    }

    my $out = Bio::SeqIO->new(-file => ">$outfile" ,
				  -format => 'Fasta');
    $logger->trace("Making temporary fasta file $outfile");

    # Grab the fna file via bioperl
    my $in = new Bio::SeqIO(-file => $filename);

    # Loop through the sequences in the fasta file
    my $found = 0;
    while(my $seq = $in->next_seq()) {
        # Split the display_id in to identifier types
        my $identifiers = $self->split_header($seq->display_id);
        
        # If the sequence has a refseq accession and it's in the
        # list of accessions we're looking for
        if($identifiers->{ref} && $identifiers->{ref} ~~ @accessions) {
            $out->write_seq($seq);
            $found++;
        }
    }

    return $found;
}


sub parse_formats {
    my $self = shift;
    my $format_str = shift;

    my $formats;
    foreach (split /\s+/, $format_str) { $_ =~ s/^\.//; $formats->{$_} = 1; }

    return $formats;
}

#used to create ptt file
sub tag {
	my ( $f, $tag ) = @_;
	return '-' unless $f->has_tag($tag);
	return join( ' ', $f->get_tag_values($tag) );
}

# Find the mean of an array of values
sub mean {
    my $self = shift;
    my $result;
    foreach (@_) { $result += $_ }
    return $result / @_;
}

sub calc_gc {
    my $self = shift;
    my $seq = $_[0];
    my $g = ( $seq =~ tr/g// );
    $g += ( $seq =~ tr/G// );
    my $c = ( $seq =~ tr/c// );
    $c += ( $seq =~ tr/C// );
    return ( $g + $c ) / length($seq);
}

# Make a temp file in our work directory and return the name
sub _make_tempfile {
    my $self = shift;
    my $workdir = shift;
    my $prefix = (@_ ? shift : 'genomeutils');

    # Let's put the file in our workdir
    my $tmp_file = mktemp(File::Spec->catpath(undef, $workdir, $prefix . "tmpXXXXXXXXXX"));
    
    # And touch it to make sure it gets made
    `touch $tmp_file`;

    return $tmp_file;
}

# Split a fasta display_id line and make a hash of
# the values based on type and value

sub split_header {
    my $self = shift;
    my $id = shift;

    my @pieces = split /\|/, $id;

    my $identifiers = {};
    my $type;
    while(($type = shift @pieces) && (my $val = shift @pieces)) {
        $identifiers->{$type} = $val;
    }

    # See if we have a coordinate in the header
    if($type =~ /:c?(\d+)\.\.(\d+)/) {
        $identifiers->{start} = $1;
        $identifiers->{end} = $2;
    }

    return $identifiers;
}

1;
