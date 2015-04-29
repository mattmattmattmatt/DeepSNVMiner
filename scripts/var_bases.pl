#! /usr/bin/perl -w

use strict;

#Use samtools to generate variant calls over a range of coordinates in bed files
my ($bed_file, $bam_file, $samtools_bin, $genome_fasta_file, $output_dir) = @ARGV;

mkdir($output_dir) if !-d $output_dir;

open(BED,"$bed_file") || modules::Exception->throw("Can't open file $bed_file\n");

while (<BED>) {
	my ($chr,$start_base,$end_base) = split();
	next unless $chr =~ /^chr/ || $chr =~ /^[0-9XYM]/;
	#($start_base, $end_base) = sort {$a <=> $b} ($start_base, $end_base);
	
	my $awk_command = q('BEGIN {OFS = FS = "\t" } ; {n=split($10,a,"") ; if(a[(pos-$4)+1] != "=" ) print $3, pos,(pos-$4)+1, a[(pos-$4)+1], $1, $4, $10 }');
	
	
	for (my $coord = $start_base; $coord <= $end_base; $coord++) {
	
	    print "Processing coord : $chr:$coord\n";
	
	    my $command = "$samtools_bin view -b $bam_file $chr:$coord\-$coord \| $samtools_bin fillmd -e - $genome_fasta_file \| grep -v \"\^\@\" \| awk -v pos=$coord $awk_command \| grep \'==\' \> $output_dir\/var_${chr}:${coord}\.out ";
	    
	    print STDERR $command . "\n";
	
	    `$command`;
	
	
	}
}

