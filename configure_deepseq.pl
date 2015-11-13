#! /usr/bin/perl -w

use strict;
use FindBin;
use lib "$FindBin::Bin";
use modules::Deepseq;
use modules::Exception;
use modules::SystemCall;
use Data::Dumper;
use File::Basename;
use Cwd 'abs_path';
use Pod::Usage;
use Getopt::Long;
use vars qw(%OPT);

GetOptions(\%OPT, 
				"help|h",
	   			"man|m",
				"conf_file=s",
				"test"
			);
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});


=pod

=head1 SYNOPSIS

configure_deepseq.pl -test generate_test_conf -conf_file generate_specific conf_file

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

configure_deepseq.pl -> Wrapper script for configure deepseq run

=head1 DESCRIPTION

date

a script that creates deepseq.conf containing path info

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

sample -f file

=cut


my (undef,$base) = fileparse(abs_path($0));
my $conf_file;

if ($OPT{test}) {
	$conf_file = defined $OPT{conf_file}?$OPT{conf_file}:$base.'/deepseq_test.conf';
} else {
	$conf_file = defined $OPT{conf_file}?$OPT{conf_file}:$base.'/deepseq.conf';
}

if (-e $conf_file) {
	modules::Exception->throw("Conf file $conf_file already exists, please remove this file before running this script");
}

#Set the path variable and generate a conf file 
print "What is the path to samtools [default=/usr/bin/samtools]?";
my $samtools = <STDIN>;
chomp $samtools;
if ( !-e $samtools && -e '/usr/bin/samtools' ) {
	$samtools =  '/usr/bin/samtools';
} elsif (!-e $samtools) {
	modules::Exception->throw("Samtools $samtools doesn't exist");	
}

print "What is the path to bwa [default=/usr/bin/bwa]?";
my $bwa = <STDIN>;
chomp $bwa;
if ( !-e $bwa && -e '/usr/bin/bwa' ) {
	$bwa =  '/usr/bin/bwa';
} elsif (!-e $bwa) {
	modules::Exception->throw("Bwa $bwa doesn't exist");	
}

my $sys_call = modules::SystemCall->new();
my $ref_fasta_abs;

if ($OPT{test}) {
	$ref_fasta_abs = abs_path('./sample/ref.fa');
} else {
	print "What is the path to single reference fasta file?";
	my $ref_fasta = <STDIN>;
	chomp $ref_fasta;
	if ( !-e $ref_fasta ) {
		modules::Exception->throw("Ref fasta file $ref_fasta doesn't exist");
	} 
	
	$ref_fasta_abs = abs_path($ref_fasta);
}

#Check the bwa index is present....
my $bwa_index_file = $ref_fasta_abs . '.sa';

if (!-e $bwa_index_file) {
	#Sometimes it ref.fa.sa but sometimes it's ref.sa so check both
	print "Can't find the index file $bwa_index_file\n" unless $OPT{test};
	($bwa_index_file = $ref_fasta_abs) =~ s/\.fa$/\.sa/;
	if (!-e $bwa_index_file) {
		print "Can't find index file $bwa_index_file\n" unless $OPT{test};
		print "Can't find any bwa index files, so we will create it. This may take a while..." unless $OPT{test};
		my $bwa_index_command = join(" ",
								$bwa,
								'index -a bwtsw',
								$ref_fasta_abs
							);
		$sys_call->run($bwa_index_command);
		$bwa_index_file = $ref_fasta_abs;		
	}
}

#bwa only requires the file prefix name
$bwa_index_file =~ s/\.sa$//;

open(CONF,">$conf_file") || modules::Exception->throw("Can't open xml file $conf_file for writing\n");
print CONF "bwa=$bwa\n";
print CONF "samtools=$samtools\n";
print CONF "ref_fasta=$ref_fasta_abs\n";
print CONF "bwa_index=$bwa_index_file\n";
close CONF;

