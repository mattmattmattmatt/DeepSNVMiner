#! /usr/bin/perl -w

use strict;
use FindBin;
use lib "$FindBin::Bin/..";
use modules::Deepseq;
use modules::Exception;
use modules::SystemCall;
use Data::Dumper;
use File::Basename;
use Cwd 'abs_path';
use Pod::Usage;
use Getopt::Long;
use vars qw(%OPT);


# Options from command line:

GetOptions(
			\%OPT, 
	   		"coord_bed=s",
	   		"ref_fasta=s",
	   		"outdir=s",
	   		"samtools=s",
	   		"bam=s",
	   		) || modules::Exception->throw("ERROR: Problem with command line arguments");;

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || (!$OPT{coord_bed} || !$OPT{samtools} || !$OPT{outdir} || !$OPT{coord_bed}) || !$OPT{bam}) ;



=pod

=head1 SYNOPSIS

var_bases.pl -samtools full_samtools_path -ref_fasta full_path_to_reference_fasta -coord_bed bed_file_containing_target_coordinates -bam full_bam_path -outdir output_directory  

Required flags: -samtools -bam -ref_fasta -outdir -coord_bed

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

var_bases.pl -> Script to run samtools calmd and parse results

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

sample -f file

=cut


#Use samtools to generate variant calls over a range of coordinates in bed files and parse the results
my $bed_file = $OPT{coord_bed};
my $bam_file = $OPT{bam};
my $samtools_bin = $OPT{samtools};
my $genome_fasta_file = $OPT{ref_fasta};
my $output_dir = $OPT{outdir};

mkdir($output_dir) if !-d $output_dir;
my $sys_call = modules::SystemCall->new();

open(BED,"$bed_file") || modules::Exception->throw("Can't open file $bed_file\n");

while (<BED>) {
	my ($chr,$start_base,$end_base) = split();
	next unless $chr =~ /^chr/ || $chr =~ /^[0-9XYM]/;
	
	#my $awk_command = q('BEGIN {OFS = FS = "\t" } ; {n=split($10,a,"") ; if(a[(pos-$4)+1] != "=" ) print $3, pos,(pos-$4)+1, a[(pos-$4)+1], $1, $4, $10 }');
	my $outfile = $output_dir. "/${chr}_${start_base}_${end_base}.fillmd";
    my $command = "$samtools_bin view -b $bam_file $chr:$start_base\-$end_base \| $samtools_bin calmd -e - $genome_fasta_file > $outfile";
	
	print STDERR $command . "\n";
	$sys_call->run($command);
	
	if (!-s $outfile) {
		modules::Exception->throw("ERROR: Can't open output file $outfile\n");
	}
	
	open(FILE,"$outfile") || modules::Exception->throw("Can't open file $outfile\n");
	
	
	my $var_file = $output_dir. "/${chr}_${start_base}_${end_base}.variants";
	
	open(VAR,">$var_file") || modules::Exception->throw("Can't open file to write $var_file\n");
	
	my %mutations = ();
	
	while (<FILE>) {
		chomp;
		next if /^@/; #skip headers
		next if /\t==+=\t/; #Skip perfect matches
		
		my @fields = split("\t");
		next if $fields[5] =~ /[SHP]/; #Only handle snvs, insertions, and deletions
		
		my @cigar_lengths = $fields[5] =~ /\d+/g;
		my @cigar_types = $fields[5] =~ /\D+/g;
		
		my $mutant_coord = $fields[3]; #Keep track of genomic coordinate
		my $base_index = 0; #Keep track of string index
		my @bases = split("",$fields[9]);
		
		for ( my $block = 0 ; $block < @cigar_types ; $block++ ) {
		    my $block_count = $cigar_lengths[$block];
		    if ($cigar_types[$block] =~ /D/) {
		    	my $start_coord = $mutant_coord;
		    	my $end_coord = $mutant_coord + $cigar_lengths[$block];
		    	print VAR join("\t",$fields[2],$start_coord,$end_coord,$base_index+1,'-'.$cigar_lengths[$block],$fields[0],$fields[3],$fields[9],$fields[5])."\n";
		    	
		    	#Here we don't change the base index counter only the genome coordinate
		    	$mutant_coord += $cigar_lengths[$block];
		    	
		    } elsif ($cigar_types[$block] =~ /I/) {
		    	my $inserted_bases;
		    	my $insertion_length = $cigar_lengths[$block];

		    	
		    	while ($insertion_length > 0) {
		    		$inserted_bases .= $bases[$base_index];
		    		$insertion_length--;
		    	}
		    	
				print VAR join("\t",$fields[2],$mutant_coord,$mutant_coord,$base_index+1,'+'.$inserted_bases,$fields[0],$fields[3],$fields[9],$fields[5])."\n";											    	

		    	#Here we don't change the genomic coordinate only the base index counter
		    	$base_index += $cigar_lengths[$block];

		    	
		    } elsif ($cigar_types[$block] =~ /M/) {
		    	while ($block_count > 0) {
		    		if ($bases[$base_index] ne '=') {
						my $mutant_base = $bases[$base_index];
						print VAR join("\t",$fields[2],$mutant_coord,$mutant_coord,$base_index+1,$mutant_base,$fields[0],$fields[3],$fields[9],$fields[5])."\n";
		    		} 
		    		
		    		$base_index++;
		    		$mutant_coord++;
		    		$block_count--;
		    	}
		    } elsif ($cigar_types[$block] =~ /N/) {
		    	#Here we don't change the base index counter only the genome coordinate
		    	$mutant_coord += $cigar_lengths[$block];
		    }
		    
		}		
		
		
		
		
		
		
			
	}
	
	
	
	#for (my $coord = $start_base; $coord <= $end_base; $coord++) {
	
	    #print "Processing coord : $chr:$coord\n";
	
		
		#my $command = "$samtools_bin view -b $bam_file $chr:$coord\-$coord \| $samtools_bin fillmd -e - $genome_fasta_file \| grep -v \"\^\@\" \| awk -v pos=$coord $awk_command \| grep \'==\' \> $output_dir\/var_${chr}:${coord}\.out ";
	    
	    
	    
	    
	
	    #`$command`;
	
	
	#}
}

