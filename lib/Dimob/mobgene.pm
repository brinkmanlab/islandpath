=head1 NAME

    Dimob::Mobgene

=head1 DESCRIPTION

    A repackaging of Will's original Dimob software

=head1 SYNOPSIS

    use Dimob::Mobgene;

    $mobgene_obj = Dimob::Mobgene->new();

    $mobgene_obj->parse_hmmer($hmmer_file, $parse_evaluecutoff);
    $mobgene_obj->parse_ptt($ptt_file);

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

package Dimob::Mobgene;

use strict;
use warnings;
use Moose;
use Log::Log4perl qw(get_logger :nowarn);

use Bio::SearchIO;

use Dimob::tabdelimitedfiles;

my $cfg; my $logger; my $cfg_file;

sub BUILD {
    my $self = shift;
    my $args = shift;

    $cfg = Dimob::Config->config;

    $logger = Log::Log4perl->get_logger;

    if($args->{extended_ids}) {
	$logger->trace("Using extended ids");
	$self->{extended_ids} = 1;
    }

}

#given a hmmer output file, parse it and return a list of mobility genes
#returns a {pid}{domain_name}{E_value} hash structure

sub parse_hmmer {
    my $self = shift;
    	my $hmmer_file         = shift;
	my $parse_evaluecutoff = shift;

	my %domain_hash;
	my %mobgenes;
#	open( HMMINPUT, $hmmer_file );

    my $search_pattern = 'gi\|(\d+)\|';
    if($self->{extended_ids}) {
	$search_pattern = 'gi\|(\d+)\|\:c?(\d+\.\.\d+)';
    }
	#making sure that the file is present and not empty
	if ( -s "$hmmer_file" ) {
		#print "parse_hmmer: parsing HMMER results...\n";
		my $in =
		  Bio::SearchIO->new( -file => "$hmmer_file", -format => 'hmmer' );
		while ( my $res = $in->next_result ) {
			my %domains_evalues;
			while ( my $hit = $res->next_hit ) {
				if ( $hit->significance <= $parse_evaluecutoff ) {
					$domains_evalues{ $hit->name } = $hit->significance;
				}
			}
			if ( ( scalar( keys %domains_evalues ) ) > 0 ) {
				my $id;
#				if ($res->query_name=~/gi\|(\d+)\|/){
				if ($res->query_name=~/$search_pattern/){
					$id = $1;
					if($2) {
					    my $coords = $2;
					    $coords =~ s/-/../;
					    $id .= "_$coords";
					}
				}
				else{
					$id = $res->query_name;
				}
				$mobgenes{ $id } = {%domains_evalues};
			}
		}
	}
	return \%mobgenes;
}

#given a ptt file and extract a list of genes (GI number and accession) that
sub parse_ptt {
    my $self = shift;
    #are annotated as mobility genes
    my $ptt_file    = shift;
    my $header_line = 3;       #currently the 3rd line of ptt file is the header

    my @cols = qw(4);
    if($self->{extended_ids}) {
	push @cols, 1;
    }
    my %mobgenes;
    my ( $header_arrayref, $pttfh ) =
	extract_headerandbodyfh( $ptt_file, $header_line );
    my $ptt_table_hashref = table2hash_rowfirst( $header_arrayref, $pttfh, @cols );
    #print "here's the dumping\n";
    #print Dumper $ptt_table_hashref;
    foreach my $pid (keys %{$ptt_table_hashref} ) {
	my $product = $ptt_table_hashref->{$pid}->{'Product'};
	if (   $product =~ /transposase/i
	       || $product =~ /IstB/i
	       || $product =~ /insertion element/i
	       || $product =~ /recombinase/i
	       || $product =~ /insertion sequence/i
	       || $product =~ /resolvase/i
	       || $product =~ /integrase/i
	       || $product =~ /phage/i
	       || $product =~ /transposon/i
	       || $product =~ /transposable element/i
		   || $product =~ /excisionase/i
	       )
	{
	    $mobgenes{$pid} = $product;
	}
    }
    return \%mobgenes;

}

1;
