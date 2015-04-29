package modules::ConfigXML;

use strict;
use XML::Simple;
use modules::Exception;
use Data::Dumper;

sub new {
    my ($class, $xmlfile, @force_array_elements) = @_;
    my $self = bless {}, $class;
    $self->{'xml_ref'} = XMLin($xmlfile, 
			       'KeyAttr' => ['name'], 
			       'ForceArray' => \@force_array_elements,
	);

	

    return $self;
}

sub _config {
    my $self = shift;

    return $self->{'xml_ref'};
}

sub read {
    my ($self, @keys) = @_;

    my $config_ref = $self->_config;

    for (my $i = 0; $i < scalar @keys; $i++) {

		if (! defined $config_ref->{$keys[$i]}){
		    modules::Exception->throw("Requested undefined config field [" . join("-->", @keys) . "]");
		}
	
		if ($i == scalar @keys - 1){
		    return $config_ref->{$keys[$i]}
		} else {
		    $config_ref = $config_ref->{$keys[$i]};
		}	
    }
}

sub exists {
    my ($self, @keys) = @_;

    my $config_ref = $self->_config;

    for (my $i = 0; $i < scalar @keys; $i++) {

		if (! defined $config_ref->{$keys[$i]}){
		    return 0;
		}
	
		if ($i == scalar @keys - 1){
		    return 1;
		} else {
		    $config_ref = $config_ref->{$keys[$i]};
		}	
    }
}

sub empty {
	my ($self, $key) = @_;
	
	my $config_ref = $self->_config;

	if (! defined $config_ref->{$key}){
	    return 0;
	} else {
		my $value = $config_ref->{$key};
		#Returns empty hash if conf file has empty string
		if (ref $value eq 'HASH') {
			return 1;
		} else {
			return 0;
		}
	}
}

sub print {
    my ($self) = @_;

    print Dumper($self->_config);
}


return 1;
