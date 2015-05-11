#!/usr/bin/perl
use strict;
use Data::Dumper;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/..";
use modules::SystemCall;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "filestub=s",
	   "min_group=i",
	   "sm_file=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{filestub} || !$OPT{min_group} || !$OPT{sm_file});


=pod

=head1 SYNOPSIS

script_name [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

script_name.pl -> One line description

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

sample -f file

=cut

my $min_group = $OPT{min_group};

my $sm_file = $OPT{sm_file};
if ( !-e $sm_file ) {
	modules::Exception->throw("File $sm_file doesn't exist");	
}

my %group_by_coord;

open(SMFILE,"$sm_file") || modules::Exception->throw("Can't open file $sm_file\n");

my $filestub = $OPT{filestub};
my $sm_single_out = $filestub . '.pass_single_supermutants.tsv';
my $sm_coord_out = $filestub . '.pass_group_supermutants.tsv';

#Load up all supermutant groups
while (<SMFILE>) {
	chomp;
    my @fields = split;
    my $unique_key = $fields[0].':'.$fields[1].':'.$fields[2].':'.$fields[3]; #chr:start:end:varbase
    $group_by_coord{$unique_key}{count}++;
    push @{$group_by_coord{$unique_key}{entries}}, $_;
}

open(SINGLE,">$sm_single_out") || modules::Exception->throw("Can't open file to write $sm_single_out\n");
open(COORD,">$sm_coord_out") || modules::Exception->throw("Can't open file to write $sm_coord_out\n");

for my $unique_key (sort {my ($a_chr,$a_coord) = split(':',$a); my ($b_chr,$b_coord) = split(':',$b); {$a_chr cmp $b_chr || $a_coord <=> $b_coord} } keys %group_by_coord) {
	if ($group_by_coord{$unique_key}{count} >= $min_group) {
		my %group_entries = ();
		for my $entry ( @{$group_by_coord{$unique_key}{entries}} ) {
			my @fields = split(" ",$entry);
			$group_entries{$fields[4] .'('.$fields[5] . '/'.$fields[6] .')'}++;
		    $entry =~ s/\s/\t/g;
		    print SINGLE $entry . "\n";
		}
		my $coord_string = join("; ",keys %group_entries);
		my ($chr,$start,$end,$base) = split(':',$unique_key);
		print COORD join("\t",$group_by_coord{$unique_key}{count},$chr,$start,$end,$base,$coord_string) ."\n";
	}
}
close SINGLE;
close COORD;










