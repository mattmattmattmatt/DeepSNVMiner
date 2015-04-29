package  modules::SystemCall;

use strict;


sub new {
    my ($class) = @_;

    my $self = bless {}, $class;

    return $self;
}

sub run {
     my ($self, $command) = @_;

    if (defined $command){
        print STDERR "Running command:\n$command\n";
        if (system($command)){
            die "Command exited with non-zero status";
        }
    } else {
        warn "No command specified to system call";
        return 0;
    }

    return 1;    
}

return 1;
