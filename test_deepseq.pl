#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;
use vars qw(%OPT);
use FindBin;
use lib "$FindBin::Bin";
use modules::SystemCall;
use modules::Exception;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "graph"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});



=pod

=head1 SYNOPSIS

test_deepseq.pl -graph generate_graphs (requires R) [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

test_deepseq.pl -> Test the installation of deepseq using sample data

=head1 DESCRIPTION

Feb 5th, 2015

Wrapper script for testing Deepseq installatino

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

#Test after running ./configure_deepseq.pl
>test_deepseq.pl

=cut

my (undef,$base) = fileparse(abs_path($0));
my $conf_file = $base.'/deepseq_test.conf';

if (!-e $conf_file) {
	modules::Exception->throw("Please run 'configure_deepseq.pl -test' first");
}


my $test_conf = $base.'/deepseq_test.conf';


my $workdir = "$base/sample/results";
my $out = $workdir.'/test.pass_group_supermutants.tsv';

my $default_command = "./run_deepseq.pl -conf_file $test_conf -filename_stub test -working_dir $workdir -read1_fastq $base/sample/read1.fq -read2_fastq $base/sample/read2.fq -coord_bed $base/sample/ref.bed";
my $sys_call = modules::SystemCall->new();

if ($OPT{graph}) {
	$default_command .= " -graph";
}

print "Testing Deepseq default....\n";
&test_deepseq($out,$default_command,82);

print "\nInstallation Success!\n";

#subroutine to check output
sub test_deepseq {
	my ($out,$command,$expected_lines) = @_;
	
	$sys_call->run($command);
	
	if ( !-e $out ) {
		modules::Exception->throw("ERROR: Installation failed; Output file $out wasn't generated");	
	} 

	if ( !-s $out) {
		modules::Exception->throw("ERROR: Installation failed; Output file $out is empty");
	}
	
	open(FILE,"$out") || modules::Exception->throw("Can't open file $out\n");
	my @lines = <FILE>;

	if (@lines != $expected_lines) {
		modules::Exception->throw("ERROR: Installation failed; Output file $out doesn't have correct number of lines");
	}
	close FILE;
}


