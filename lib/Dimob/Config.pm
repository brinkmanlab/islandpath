=head1 DESCRIPTION

    A reuse of Matthew's original code for IslandViewer
    to use a config file the main variables needed for the pipeline

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

package Dimob::Config;
use MooseX::Singleton;

use strict;
use warnings;
use Config::Simple;
use Data::Dumper;

has config => (
    is     => 'ro',
    isa    => 'Ref',
    writer => '_set_config'
);

has config_file => (
    is     =>  'ro',
    isa    =>  'Str',
    writer =>  '_set_file'
);

sub initialize {
    my $self = shift;
    my $args = shift;

    my $cfg_file = $args->{cfg_file};

    die "Error, unable to read config file $cfg_file"
	unless(-f $cfg_file && -r $cfg_file);

    # Read the section [main] of the config file
    my $config = new Config::Simple($cfg_file)->param(-block => 'main');

    $self->_set_config($config);

    # Expand the variables in the config file
    $self->evaluate_parameters();

    # Save the file name so we can pass it to
    # helper scripts
    $self->_set_file($cfg_file);

    return $self;
}

# Go through the config variables and do
# substitutions as needed

sub evaluate_parameters {
    my $self = shift;

    my $config = $self->config;

    for my $param (keys %{$config}) {
	if($config->{$param} =~ /\{\{.+\}\}/) {
	    $config->{$param} =~ s/\{\{([\w_]+)\}\}/$config->{$1}/eg;
	    
	}
    }
}

# Given a path, expand out and shortened
# paths

sub expand_directory {
    my $self = shift;
    my $filename = shift;

    # Get our local config
    my $cfg = $self->config;

    # Expand filename
    if($filename =~ /\{\{.+\}\}/) {
	$filename =~ s/\{\{([\w_]+)\}\}/$cfg->{$1}/eg;
    }

    return $filename
}

# Give a file path and filename and try to
# shorten it for the user

sub shorten_directory {
    my $self = shift;
    my $path = shift;

    # Get our local config
    my $cfg = $self->config;

    $path =~ s/\/\//\//g;

    if($cfg->{rootdir} &&
	    $path =~ /$cfg->{rootdir}/) {
	$path =~ s/$cfg->{rootdir}/{{rootdir}}/;
    }

    return $path;
}

1;
