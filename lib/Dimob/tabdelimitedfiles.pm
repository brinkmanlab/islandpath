=head1 DESCRIPTION

    A repackaging of Will's original code - its tried and true, why spend time rewriting it
    loosely associated subroutines used by Dimob.pm to deal with tabdelimited files

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

    Oct 24, 2016

=cut

package Dimob::tabdelimitedfiles;

#a package containing some subroutines to deal with tabdelimited files
use strict;
use warnings;
use Carp;
our ( @ISA, @EXPORT, @EXPORT_OK );    
use Exporter;
@ISA    = qw(Exporter);
@EXPORT =
  qw(extract_header extract_headerandbodyfh columns2hash combinehashes mergecolumns table2hash_rowfirst table2hash_columnfirst rows2hash get_row_by_key get_column_by_key nested_hash2twoDarray sort_twoDarray);
@EXPORT_OK =
  qw(extract_header extract_headerandbodyfh columns2hash combinehashes mergecolumns table2hash_rowfirst table2hash_columnfirst rows2hash get_row_by_key get_column_by_key nested_hash2twoDarray sort_twoDarray);
sub extract_header {

	#assume that the first line of the input file is the header and return an array of header names
	#as reference
	#optionally, user can specify which line contains the header
	my ( $inputfilename, $header_linenum ) = @_;
	if ( $header_linenum == undef ) { $header_linenum = 1; }
	my @headers;
	my $string;
	$. = 0;
	open my $inputfilehandle, '<', $inputfilename
	  or croak "Couldnot open $inputfilename";
	do { chomp( $string = <$inputfilehandle> ) }
	  until $. == $header_linenum || eof;
	@headers = split /\t/, $string;
	croak "header is empty\n" if ( ( scalar @headers ) == 0 );
	return \@headers;
}

sub extract_headerandbodyfh {

	#assume that the first line of the input file is the header if none provided
#and return an array of headers from that line
#return the rest of the file as a variable to the file handle. All lines before the header line
#will be discarded.
	my ( $inputfilename, $header_linenum ) = @_;
	if ( $header_linenum == undef ) { $header_linenum = 1; }
	my @headers;
	my $string;
	$. = 0;
	open my $inputfilehandle, '<', $inputfilename
	  or croak "Couldnot open $inputfilename";
	do { chomp( $string = <$inputfilehandle> ) }
	  until $. == $header_linenum || eof;
	@headers = split /\t/, $string;
	croak "header is empty\n" if ( scalar @headers == 0 );
	return \@headers, $inputfilehandle
	  ;    #note that at this point all lines before the header line
	       #are stripped
}

sub columns2hash {

#convert columns to hash (of arrays) using the headers as keys
#inputs: reference to an array of headers; input filehandle name (headers removed)
	my @headers = @{ shift @_ };
	my $fd_name = shift;
	my %table_content;
	while (<$fd_name>) {

		#ignore any blank lines
		next if ( is_blank($_) );
		chomp;
		my @content = split /\t/, $_;

		#make sure that @content and @headers have the same number of elements
		if ( scalar(@headers) != scalar(@content) ) {
			croak "the number of header elements do not match the number of content elements";
		}
		my $i = 0;
		foreach (@headers) {
			push @{ $table_content{$_} }, $content[$i];
			$i++;
		}
	}
	return \%table_content;
}

sub rows2hash {

	#convert rows to hash (of arrays) using a specified column as the keys
	#inputs: column number starts at 1; input filehandle name (headers removed)
	my $key_column_num = shift;
	if ( $key_column_num < 1 ) {
		croak "the column number to use as the hash keys starts at 1";
	}
	my $key_index = $key_column_num - 1;
	my $fd_name   = shift;
	my %table_content;
	while (<$fd_name>) {
		next if ( is_blank($_) );
		chomp;
		my @content = split /\t/, $_;

		#make sure that the key column exists
		croak "the column number used does not exist"
		  if ( ($key_column_num) > scalar(@content) );
		if ( exists $table_content{ $content[$key_index] } ) {
			croak "Key: $content[$key_index] is not unique!";
		}
		my $i = 0;
		foreach (@content) {
			push @{ $table_content{ $content[$key_index] } }, $content[$i];
			$i++;
		}
	}
	return \%table_content;
}

sub table2hash_columnfirst {

#convert the table to a hash based data structure using headers as the primary keys
#and a user defined column as secondary keys
#inputs: reference to an array of headers; column number starts at 1; input filehandle name
	my @headers        = @{ shift @_ };
	my $key_column_num = shift;
	if ( $key_column_num < 1 ) {
		croak "the column number to use as the hash keys starts at 1";
	}
	if ( $key_column_num > scalar(@headers) ) {
		croak "the column number used does not exist";
	}
	my $key_index = $key_column_num - 1;
	my $fd_name   = shift;

	my %table_content;
	while (<$fd_name>) {

		#ignore any blank lines
		next if ( is_blank($_) );
		chomp;
		my @content = split /\t/, $_;

		#make sure that @content and @headers have the same number of elements
		croak
"the number of header elements do not match the number of content elements"
		  if ( scalar(@headers) != scalar(@content) );
		my $i = 0;
		foreach my $header (@headers) {
			if ( exists $table_content{$header}{ $content[$key_index] } ) {
				croak "Key: $content[$key_index] is not unique!";
			}
			$table_content{$header}{ $content[$key_index] } = $content[$i];
			$i++;
		}
	}
	return \%table_content;
}

sub table2hash_rowfirst {

#convert the table to a hash based data structure using headers as the primary keys
#and a user defined column as secondary keys
#inputs: reference to an array of headers; column number starts at 1; input filehandle name
	my @headers        = @{ shift @_ };
	#print "@headers\n";
	my $fd_name   = shift;
	my @key_column_nums = @_;
	for my $key (@key_column_nums) {

	    unless ( $key >= 1 && $key <= scalar(@headers) ) {
		croak "the column number used [$key] does not exist";
	    }
	}
	my @key_indexes = map { $_ - 1 } @key_column_nums;
	my %table_content;
	while (<$fd_name>) {

		#ignore any blank lines
		next if ( is_blank($_) );
		chomp;
                # Set limit to -1 so we get rows where the last column (product) is empty
		my @content = split /\t/, $_, -1;
		#make sure that @content and @headers have the same number of elements
		croak "the number of header elements do not match the number of content elements, row: '$_'"
		  if ( scalar(@headers) != scalar(@content) );
		my $i = 0;
#		my $key = $content[$key_index];
		my $key = join('_', map { $content[$_] } @key_indexes);
#		$key =~ s/\.\./-/g;
		foreach my $field (@content) {
			if ( exists $table_content{ $key }{ $headers[$i] } )
			{
				croak "Key: $key is not unique!";
			}
			$table_content{ $key }{ $headers[$i] } = $field;
			$i++;
		}
	}
	return \%table_content;
}

sub mergecolumns {

#take an arbitary number of hashes of hashes or arrays and combine the hashes based on
#the primary hash key
	my $reference_hash  = shift;    #use this hash's key as the reference
	my @hashes_to_merge = @_;
	my %finalhash;
	foreach my $key ( keys %$reference_hash ) {
		$finalhash{$key} = $reference_hash->{$key};
		foreach my $hash (@hashes_to_merge) {
			if ( exists $hash->{$key} ) {
				my ( $k, $v );
				while ( ( $k, $v ) = each( %{ $hash->{$key} } ) ) {
					croak "key already exists in the hash"
					  if ( exists $finalhash{$key}{$k} );
					$finalhash{$key}{$k} = $v;
				}
			}

		}
	}
	return \%finalhash;
}

sub combinehashes {

	#combine two or more hashes
	#if identical keys are present in the hashes, the procedure will fail
	my @hashes_to_merge = @_;
	my %finalhash;
	foreach (@hashes_to_merge) {
		my %hash = %{$_};
		my $key;
		my $value;
		while ( ( $key, $value ) = each(%hash) ) {
			if ( exists $finalhash{$key} ) {
				croak "Duplicated key exists in the hashes being merged!";
			}
			$finalhash{$key} = $value;
		}
	}
	return \%finalhash;
}

sub nested_hash2twoDarray {

#take a nested hash and use the primary hash keys as the row lable
#return a two dimensional array
#the user provides a list of headers as an array ref and any arbitrary text to include at
#the beginning of the output
	my $hash_to_print = shift;    # a hash ref
	my $header_ref    = shift;    # an array ref
	my $titletext     = shift;    #a string
	chomp $titletext;
	my @twoDarray;
	push @twoDarray, ["$titletext"];
	push @twoDarray, [@$header_ref];

	foreach my $key ( keys %$hash_to_print ) {
		my @temparray;
		foreach my $datatype (@$header_ref) {
			push @temparray, $hash_to_print->{$key}{$datatype};
		}
		push @twoDarray, [@temparray];
	}
	return \@twoDarray;
}

sub sort_twoDarray{
#take a twoDarray and sort it according to values from a given column
#user can specify how many lines to skip before sorting begins
	my $twoDarray = shift;
	my $column = shift; #the column to use as the key in sorting
	my $linestoskip = shift || 0;
	my $row_size = scalar(@$twoDarray);
	croak "lines to skip starts at 0 and must be smaller than the total number of rows" if ($linestoskip < 0 || $linestoskip > $row_size);
	my $column_size = scalar(@{$twoDarray->[-1]});
	croak "column number starts at 1 and must be smaller than the total number of columns" if ($column < 1 || $column > $column_size);
	my $index = $column -1;
	my @finalarray;
	while ($linestoskip){
		my $line = shift @$twoDarray;
		push @finalarray, $line;
		$linestoskip--;
	}
	for my $arrayslice ( sort { $a->[$index] cmp $b->[$index] } @$twoDarray ) {
    	push @finalarray, $arrayslice;
	}
	return \@finalarray;
}

sub get_column_by_key {

	#given a hash of arrays and a key, return the array	with matching key
	my $hash_ref = shift;
	my $key      = shift;
	if (wantarray) { return @{ $hash_ref->{$key} }; }
	return $hash_ref->{$key};
}

sub get_row_by_key {

	#given a hash of arrays and a key, return the array with matching key
	my $hash_ref = shift;
	my $key      = shift;
	if (wantarray) { return @{ $hash_ref->{$key} }; }
	return $hash_ref->{$key};
}

sub array2hash {
	my $array_ref = shift;
	my %hash;
	foreach (@$array_ref) {
		$hash{$_} = 1;
	}
	return \%hash;
}

sub is_blank {

	#check if a string is empty
	my $string = shift;
	if ( $string =~ /^\s*$/ ) { return 1; }    #true
	else { return;}                            #false
}
