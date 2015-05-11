#!/usr/bin/perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/..";
use modules::Exception;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);

	GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "read1_file=s",
	   "read2_file=s",
	   "uid_len1=i",
	   "uid_len2=i",
	   "adaptor_seq=s",
	   "min_seqlen=i",
	   "no_uid"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{read1_file} || !$OPT{read2_file} || !$OPT{uid_len1});



=pod

=head1 SYNOPSIS

group_by_barcode.pl -read1_file read1_file -read2_file read2_file -uid_len1 uid_length_from_5`_end -uid_len2 uid_length_from_3`_end -adaptor filter_out_these_adaptor_seqs(divided_by_comma) -no_uid no_uid_present_in_sequence 

Required flags: -read1_file -read2_file

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

my $read1 = $OPT{read1_file};
my $read2 = $OPT{read2_file};

my $uid_len1 = $OPT{uid_len1};
my $uid_len2 = $OPT{uid_len2};

my $no_uid = defined $OPT{no_uid}?1:0;
my $min_seqlen = defined $OPT{min_seqlen}?$OPT{min_seqlen}:0;

#If there adaptor we need to filter out
my $adaptor = 0;
my @adaptors;
if (defined $OPT{adaptor_seq}) {
	$adaptor = 1;
	@adaptors = split(",",$OPT{adaptor_seq});
}

#Flag to tell whether we need 
my $split_barcode = $uid_len2 > 0?1:0;

open(READ1,"$read1") || die "Can't open file $read1\n";
open(READ2,"$read2") || die "Can't open file $read2\n";

my $count = 0;
my $any_good_data = 0;
my $tag1;
my $tag2;
my $match_seq_length;

my $match_count = 0;
my $no_match = 0;

my %barcode_count = ();

open(FILTERREAD1,">$read1.filter") || die "Can't open $read1.filter";
open(FILTERREAD2,">$read2.filter") || die "Can't open $read2.filter";

#First create a cleaned up fastq file w/o barcodes and the barcode in the tag
while (!eof(READ1) && !eof(READ2)) {
	my $line1 = <READ1>;
	my $line2 = <READ2>;
	chomp $line1;
	chomp $line2;
	
	if ($count%4 == 0) { #tag line
	     # save sequence tag and ignore anything else 
		($tag1) = $line1 =~ /(\S+)/;
		($tag2) = $line2 =~ /(\S+)/;
		# reset global flags
		$any_good_data = 0;
		$match_seq_length = 0;
	} 
	
	if ($count%4 == 1) { #sequence line

		#Skip all N lines
		if ($line1 =~ /^N.*N$/ ||  $line2 =~ /^N.*N$/) {
			$count++;
			next;
		}
		
		#Get the revcom
		my $revcom_line2 = revcom($line2);

		#First remove any adaptors if present
		if ($adaptor) {
			for my $adaptor_seq (@adaptors) {
				if ($line1 =~ /^$adaptor_seq/) {
					my $filtered_line1 = &remove_adaptor($line1,$adaptor_seq);
					$line1 = $filtered_line1;
				}

				if ($revcom_line2 =~ /^$adaptor_seq/) {
					my $filtered_line2 = &remove_adaptor($revcom_line2,$adaptor_seq);
					$revcom_line2 = $filtered_line2;
				}
			}
			
			
			#if (length($line1) != length($revcom_line2)) {
				#modules::Exception->throw("ERROR: Filtering made read lengths different\n$line1\n$revcom_line2");
			#}
			
		}

		my $final_seq1 = my $final_seq2 = my $barcode1 = my $barcode2;
		if ($split_barcode) {
			#uid from both ends
			my ($bar1_part1,$seq1,$bar1_part2) = $line1 =~ /^(\S{$uid_len1})(\S+)(\S{$uid_len2})/;
			my ($bar2_part1,$seq2,$bar2_part2) = $revcom_line2 =~ /^(\S{$uid_len1})(\S+)(\S{$uid_len2})/;
			$barcode1 = $bar1_part1.$bar1_part2;
			$barcode2 = $bar2_part1.$bar2_part2;
			
			if ($no_uid) { #Without uids just use the sequence
				$final_seq1 = $line1;
				$final_seq2 = $line2;
			} elsif ($barcode1 eq $barcode2) { #With uids need to ensure barcodes match
				$final_seq1 = $seq1;
				$final_seq2 = $seq2;
			}
			
		} else {
			#uid from single end
			my ($bar1,$seq1) = $line1 =~ /^(\S{$uid_len1})(\S+)/;
			my ($bar2,$seq2) = $revcom_line2 =~ /^(\S{$uid_len1})(\S+)/;
			$barcode1 = $bar1;
			$barcode2 = $bar2;
			
			if ($no_uid) {
				$final_seq1 = $line1;
				$final_seq2 = $line2;
			} elsif ($barcode1 eq $barcode2) {
				$final_seq1 = $seq1;
				$final_seq2 = $seq2;
			} 
			
		}

		
	

	    # On with counting a match, presuming there was one

	    if (defined $final_seq1 && defined $final_seq2) { # if these are defined, the barcodes between read1 and read2 matched or -no_uid flag is used
			if (length($final_seq1) >= $min_seqlen && length($final_seq2) >= $min_seqlen) {
				#Check the sequence content left is long enough
				$match_count++;
				$match_seq_length = length($final_seq1);
				$any_good_data = 1;
				print FILTERREAD1 "$tag1:UID=$barcode1\n$final_seq1\n";
				print FILTERREAD2 "$tag2:UID=$barcode2\n$final_seq2\n";
				
			}
		} else {
			$no_match++;
			#print "NO MATCH:\n$line1 ($length1)\n$line2 ($length2)\n\n";
		}
			
	} 
	
	if ($any_good_data) {
		if ($count%4 == 2) {
			#'+' line
			print FILTERREAD1 "$line1\n";
			print FILTERREAD2 "$line2\n";
		} elsif ($count%4 == 3) {
			#Quality string; we need to shorten this to match the sequence bases
			if ($no_uid) {
				print FILTERREAD1 "$line1\n";
				print FILTERREAD2 "$line2\n";
			} else {
				my ($readqual1) = $line1 =~ /^\S{$uid_len1}(\S{$match_seq_length})/;
				my ($readqual2) = $line2 =~ /^\S{$uid_len2}(\S{$match_seq_length})/;
				print FILTERREAD1 "$readqual1\n";
				print FILTERREAD2 "$readqual2\n";
			}
		}
	}
	
	$count++;
	if ($count%1000000 ==0) {
		print "Match $match_count\n" unless $no_uid;
		print "No match: $no_match\n\n" unless $no_uid;
	}
	
	#last if $count %100000 == 0;
}


print "Match $match_count\n";
print "No match: $no_match\n";



### Subroutines

sub revcom {
	my $string = shift;

	# Reverse-complement now
	$string =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	$string = reverse $string;

	return $string;
}

#Remove whatever portion of the adaptor we find
sub remove_adaptor {
	my ($line,$adaptor_seq) = @_;
	
	my $adaptor_len = length($adaptor_seq);
	my $return_line;
	#Iterate over the length of the adaptor to find longest match
	for ( my $count = 0 ; $count < $adaptor_len ; $count++ ) {
	    my ($adaptor_chunk) = substr($adaptor_seq,$count,$adaptor_len);
	    if ($line =~ /^$adaptor_chunk/) {
	    	$line =~ s/$adaptor_chunk//;
	    	return $line;
	    }
	}
	return $line;
}

