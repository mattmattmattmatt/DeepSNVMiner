#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/..";
use Pod::Usage;
use modules::Exception;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "var_dir=s",
	   "group_file=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{var_dir} || !$OPT{group_file});



=pod

=head1 SYNOPSIS

var_summary.pl -var_dir directory_with_files_for_each_coord -group_file file_containing_group_counts 

Required flags: -var_dir -group_file

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

my $var_dir = $OPT{var_dir};

if (!-d $var_dir) {
	modules::Exception->throw("ERROR: Directory $var_dir doesn't exist");
}

my %group_counts = ();

my $group_file = $OPT{group_file};


if ( !-e $group_file ) {
	modules::Exception->throw("File $group_file doesn't exist");	
}

open(GROUP,"$group_file") || modules::Exception->throw("Can't open file $group_file\n");

while (<GROUP>) {
	chomp;
	my @fields = split();
	my ($barcode) = $fields[1] =~ /([ACGT]+)/;
	$group_counts{$barcode} = $fields[0]; 
#	print STDOUT "Line:$_\tField0:".$fields[0]."\tField1:".$fields[1]."\tBarcode:$barcode\n";
#	die "Argh\n";
}

opendir (DIR,$var_dir) || modules::Exception->throw("ERROR: Can't open directory");
my @files = grep {/\d+/} readdir DIR;
closedir DIR;

open(SUMMARY,">variant_summary.txt") || modules::Exception->throw("ERROR: Can't open variant_summary file");

for my $file (sort {my ($a_chr,$a_coord) = $a =~ /(\S+)_(\d+)/;my ($b_chr,$b_coord) = $b =~ /(\S+)_(\d+)/; $a_chr cmp $b_chr || $a_coord <=> $b_coord } @files) {
	#my ($chr,$coord) = $file =~ /([0-9XYM]+):(\d+)/;
	next if $file =~ /fillmd/;
	open(FILE,"$var_dir/$file") || modules::Exception->throw("Can't open file $var_dir/$file\n");
	my %var_data = ();
	
	#Keep track of duplicate lines which sometimes appear with fillmd..... -> couldn't duplicate problem and breaks handling of read pairs containing same variant
	#my %lines = (); 
	
	while (<FILE>) {
		chomp;
		#next if exists $lines{$_};
		#$lines{$_}++;
		my @fields = split("\t");

		#Only report variants in the block
		my ($file_start,$file_end) = $file =~ /_(\d+)_(\d+)/;
		next unless $file_start < $fields[1] && $file_end > $fields[1];
		
		my @tag_fields = split(":",$fields[5]);
		my $barcode = $tag_fields[-1];
		$barcode =~ s/UID=//; 



		next if $fields[3] eq ''; # added by DA - this field is frequently empty
		
		
		
		my $non_match_count = () = $fields[7] =~ /[ATCG]/g; #skip reads with 10 or more mismatches
		
		if ($non_match_count >= 10) {
			next;
		}

		next if $fields[7] =~ /[ATCG]{5}/; #Skip reads with 5 or more consecutive mismatches

		#Finally skip variants reported within the first or last 3 bases
		if ($fields[3] <= 3 || (length($fields[7])-$fields[3]) <= 3) {
			next;
		}
		
		
		$var_data{$fields[0]}{$fields[1]}{$fields[2]}{$barcode}{$fields[4]}++;	
		
	}
	
	close FILE;
	
	for my $chr ( sort keys %var_data ) {
	    for my $start ( sort {$a<=>$b} keys %{$var_data{$chr}} ) {
	    	for my $end (keys %{$var_data{$chr}{$start}}) {
				for my $barcode (sort keys %{$var_data{$chr}{$start}{$end}}) {
					for my $base (sort {$a cmp $b} keys %{$var_data{$chr}{$start}{$end}{$barcode}}) {
						my $var_count = $var_data{$chr}{$start}{$end}{$barcode}{$base};
						my $group_count = $group_counts{$barcode};
						
						#Secondary alignments sometimes get counted twice
						$var_count = $group_count if $var_count > $group_count;
						
						if (!exists $group_counts{$barcode}) {
							#modules::Exception->warning("ERROR: Problem with barcode $barcode, skipping");
							next;
						} 
						
						my $percent = sprintf("%.2f",($var_count / $group_count) * 100) if $group_count;
						print SUMMARY "$chr $start $end $base $barcode $var_count $group_count $percent%\n";
					}
				}
	    		
	    	}
		}
	}
	
	
}

#These files are often large and not needed after this
#system("rm $var_dir/*fillmd");

close SUMMARY;




