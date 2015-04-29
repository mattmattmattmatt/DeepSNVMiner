package modules::Deepseq;
use File::Basename;
use Data::Dumper;
use strict;


#Contains the steps to run grouped into 6 main blocks: check_fastq,bwa,samtools,call_variants,reports,and graph
my %COMMANDS = (1 => 
				  	 {
				  	 	'block' => 'check_fastq',
				  	 	'commands' => [
				  	 					q(sed -n '2~4p' FASTQFILE1 |wc -l), # Total reads in FASTQ file 1
				  	 					q(sed -n '2~4p' FASTQFILE2 |wc -l),	# Total reads in FASTQ file 2	
				  	 					q(sed -n '2~4p' FASTQFILE1 | grep ^N.*N$ | wc -l),# Count all-N reads in FASTQ file 1
				  	 					q(sed -n '2~4p' FASTQFILE2 | grep ^N.*N$ | wc -l),# Count all-N reads in FASTQ file 2
				  	 					q(sed -n '2~4p' FASTQFILE1 | grep -v ^N.*N$ > FASTQFILE1.no_Ns_R1),# Produce FASTQ file with all-N reads removed (read 1)
				  	 					q(sed -n '2~4p' FASTQFILE2 | grep -v ^N.*N$ > FASTQFILE2.no_Ns_R1),# Produce FASTQ file with all-N reads removed (read 2)
				  	 					q(cat FASTQFILE1.no_Ns_R1 FASTQFILE2.no_Ns_R1 | awk '{ print length($0); }' > FILENAMESTUB.no_Ns_sizes),# Produce a file containing sizes of all non-N reads
				  	 					],
				  	 	#Reports to stdouts (1 = STDOUT; 0 = STDERR)
				  	 	'stdout' => [
				  	 				1,
				  	 				1,
				  	 				1,
				  	 				1,
				  	 				0,
				  	 				0,
				  	 				0
				  	 				],
						
				    },
							  
				2 => 						    
					{
				  	 	'block' => 'pool_reads',
				  	 	'commands' => [
										"GROUPBYBARCODE -read1_file FASTQFILE1 -read2_file FASTQFILE2 ADAPTOR NO_UID -uid_len1 UID_LEN1 -uid_len2 UID_LEN2 MIN_SEQLEN",# Filter the reads, remove the barcode sequence content and add barcodes to the tags
										#q(perl -e '$wcout = `wc -l  FASTQFILE1.filter`; $wcout =~ /^(\d+)/; $len = $1 - 4; `head -$len FASTQFILE1.filter > temp.fastq.filter`; `mv temp.fastq.filter FASTQFILE1.filter`;'),# monumental hack - to remove last line of filtered fastq for Yogesh's data, which seems to always have a bizarre character or summat at the last read of the file
				  	 					#q(perl -e '$wcout = `wc -l  FASTQFILE2.filter`; $wcout =~ /^(\d+)/; $len = $1 - 4; `head -$len FASTQFILE2.filter > temp.fastq.filter`; `mv temp.fastq.filter FASTQFILE2.filter`;'),
				  	 					],
				  	 	#Reports to stdouts (1 = STDOUT; 0 = STDERR)
				  	 	'stdout' => [
				  	 				0,
				  	 				#0,
				  	 				#0
				  	 				],
						
				    },		   
							    
				3 => 
				  	 {
				  	 	'block' => 'bwa',
				  	 	'commands' => [
				  	 					q(BWA aln -t NUMTHREADS -O 5 -M 10 -i 0 BINDEX FASTQFILE1.filter > FASTQFILE1.sai),# Align reads 1
				  	 					q(BWA aln -t NUMTHREADS -O 5 -M 10 -i 0 BINDEX FASTQFILE2.filter > FASTQFILE2.sai),# Align reads 2
				  	 					q(BWA sampe BINDEX FASTQFILE1.sai FASTQFILE2.sai FASTQFILE1.filter FASTQFILE2.filter > FILENAMESTUB.sam),# Merge alignment ...
				  	 					q(grep -v XT:A:R FILENAMESTUB.sam | SAMTOOLS view -S -b - | SAMTOOLS sort -m1800000000 - FILENAMESTUB),# Massage and sort alignment
				  	 					q(SAMTOOLS index FILENAMESTUB.bam )# Generate alignment index *.bai file
				  	 					],
				  	 	#Reports to stdouts (1 = STDOUT; 0 = STDERR)
				  	 	'stdout' => [
				  	 				0,
				  	 				0,
				  	 				0,
				  	 				0,
				  	 				0
				  	 				],
						
				    },			 
				
				4 => 
				  	 {
				  	 	'block' => 'call_variants',
				  	 	'commands' => [
				  	 					q(VARBASES BEDFILE FILENAMESTUB.bam SAMTOOLS GENOMEFASTA vars_FILENAMESTUB), # Get the variant bases for each read
				  	 					],
				  	 	#Reports to stdouts (1 = STDOUT; 0 = STDERR)
				  	 	'stdout' => [
				  	 				0
				  	 				],
						
				    },	
				    		
				5 => 
				  	 {
				  	 	'block' => 'report',
				  	 	'commands' => [
				  	 					q(grep UID= FASTQFILE1.filter FASTQFILE2.filter | cut -d':' -f9 | sed -e 's/\/1//' -e 's/UID=//' | sort | uniq -c |  perl -ne 's/^\s+//; s/[\t\s]+/ /g; print $_ . "\n"' > FILENAMESTUB_R1.group_counts_final),# Generate stats - number of reads in groups
				  	 					q(VARSUMMARY -var_dir vars_FILENAMESTUB -group_file FILENAMESTUB_R1.group_counts_final),# Group by barcode for each coordinate
										q(cut -d' ' -f1 FILENAMESTUB_R1.group_counts_final |sort | uniq -c | perl -e 'while(<>){chomp; /(\d+)\s+(\d+)/; $tally += $1;} print $tally . "\n";'),#Total reads in groups
										q(grep -w ^1 FILENAMESTUB_R1.group_counts_final |wc -l),#Total reads in singleton groups
										q(perl -ne '@cols = split /\s/; print if $cols[5] >= SM_COUNT && $cols[4]/$cols[5] >= SM_PORTION;' variant_summary.txt > FILENAMESTUB.single_supermutants.txt),#Commands to filter all supermutants and count them
										q(FINALSUMMARY -sm_file FILENAMESTUB.single_supermutants.txt -filestub FILENAMESTUB -min_group MIN_GROUP),#Commands to group supermutants by coord and apply min uid filter 
										q(wc -l FILENAMESTUB.pass_group_supermutants.tsv) #Number of supermutant groups (defaults: 10 or more members per UID; 40% or more of bases mutant; 2 barcode groups at least)
				  	 					],
				  	 	#Reports to stdouts (1 = STDOUT; 0 = STDERR)
				  	 	'stdout' => [
				  	 				0,
				  	 				0,
				  	 				1,
				  	 				1,
				  	 				0,
				  	 				0,
				  	 				1
				  	 				],
						
				    },	
				6 => 
				  	 {
				  	 	'block' => 'graph',
				  	 	'commands' => [
										q(cut -d' ' -f1 FILENAMESTUB_R1.group_counts_final | grep -v -w ^1 > no_Ns_group_sizes;  Rscript R_DIR/group_size.R WORKING_DIR FILENAMESTUB_R1.group_counts_final; ls group_sizes.graph.pdf),#Graph of group sizes
										q(wc -l vars_FILENAMESTUB/var_* | perl -ne 'chomp; s/^\s+//; s/[\t\s]+/ /g; /(\d+).+var_(\S+):(\d+).out/; print join("\t", $2, $3, $1) . "\n";' > variant_coord;),
										q(GRAPH -r_dir R_DIR -variants WORKING_DIR/variant_coord -coord_bed BEDFILE -super_mutants WORKING_DIR/FILENAMESTUB.pass_single_supermutants.tsv) 
										#perl -e 'while(<>) {@cols = split /\s+/; $coord_tally{$cols[0]}++;} foreach $coord ( (10..20) ) {print join ("\t", $coord, defined($coord_tally{$coord})?$coord_tally{$coord}:0) ."\n"; }' FILENAMESTUB.passing_supermutants.txt > supermutant_coords.txt; R CMD BATCH R_DIR/supermut_coord.R; ls supermut_coords.graph.pdf) # Distribution of supermutants graph
										# R CMD BATCH R_DIR/mut_coord.R; ls mut_coords.graph.pdf),# Distribution of all mutations graph
				  	 					#q(R CMD BATCH R_DIR/supermut_coord.R; ls supermut_coords.graph.pdf) # Distribution of supermutants graph
				  	 					],
				  	 	#Reports to stdouts (1 = STDOUT; 0 = STDERR)
				  	 	'stdout' => [
				  	 				1,
				  	 				1,
				  	 				0
				  	 				],
						
				    },		
				);



sub new {
	my ($class, %args) = @_;
	
	my $self = bless {}, $class;

    my @required_args = (
    					-filename_stub,
	   					-read1_fastq,
	   					-read2_fastq,
	   					-working_dir,
	   					-uid_len1,
	   					-uid_len2,
	   					-coord_bed,
	   					-bwa,
	   					-samtools,
	   					-ref_fasta,
	   					-base
    					);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set $args{$required_arg}");
		} 
    }
   	
	my $num_threads = defined $args{-num_threads}?$args{-num_threads}:1;
	   	
   	#Set the regexes for replacement
   	my %mapping = ();
	my $bin_dir = $args{-base} .'/bin';
	my $r_dir = $args{-base} .'/R';
   	my $scripts_dir = $args{-base} . '/scripts';
   	
   	#Required args first ->  these variables replace regexes in the commad blocks above
   	$mapping{BWA} = $args{-bwa};
   	$mapping{WORKING_DIR} = $args{-working_dir};
	$mapping{FILENAMESTUB} = $args{-filename_stub};
	$mapping{FASTQFILE1} = $args{-read1_fastq};
	$mapping{FASTQFILE2} = $args{-read2_fastq};
   	$mapping{UID_LEN1} = $args{-uid_len1};
   	$mapping{UID_LEN2} = $args{-uid_len2};
   	$mapping{BEDFILE} = $args{-coord_bed};
   	$mapping{SAMTOOLS}= $args{-samtools};
   	$mapping{GENOMEFASTA} = $args{-ref_fasta};
	(my $bwa_index = $args{-ref_fasta}) =~ s/\.fa//;
	$mapping{BINDEX} = $bwa_index;
	
	
   	$mapping{R_DIR}  =  $r_dir;
   	$mapping{NUMTHREADS} = defined $args{-num_threads}?$args{-num_threads}:1;
   	$mapping{SM_COUNT} = defined $args{-sm_count}?$args{-sm_count}:10;
   	$mapping{SM_PORTION} = defined $args{-sm_portion}?$args{-sm_portion}:0.9;
   	$mapping{MIN_GROUP} = defined $args{-min_group}?$args{-min_group}:2;
   	$mapping{NO_UID} = defined $args{-no_uid}?'-no_uid':'';
   	$mapping{MIN_SEQLEN} = defined $args{-min_seqlen}?"-min_seqlen ".$args{-min_seqlen}:'';
   	$mapping{GROUPBYBARCODE} = "$scripts_dir/group_by_barcode.pl";
   	$mapping{VARBASES}  = "$scripts_dir/var_bases.pl";
   	$mapping{VARSUMMARY}  = "$scripts_dir/var_summary.pl";
   	$mapping{FINALSUMMARY} = "$scripts_dir/final_summary.pl";
   	$mapping{GRAPH} = "$scripts_dir/graph.pl";
   	$mapping{ADAPTOR} = ''; #Gets reset later if we find adaptor sequence
   
   	$self->{mapping} = \%mapping;
   	
   	return $self;
}

#Get the block number to start from
sub get_step_number {
	my ($self,$step) = @_;
	for my $command_count (sort {$a<=>$b} keys %COMMANDS) {
		if ($COMMANDS{$command_count}{block} eq $step) {
			return $command_count;
		}
	}
	return 1;
}

#Check the step to block to start from exists
sub check_step {
	my ($self,$step) = @_;
	for my $command_count (sort {$a<=>$b} keys %COMMANDS) {
		if ($COMMANDS{$command_count}{block} eq $step) {
			return 1;
		}
	}
	return 0;
}


#If there's an adaptor sequence we need to pass it in to group_by_barcode
sub set_adaptor {
	my ($self,$adaptor_seq) = @_;
	$self->{mapping}->{ADAPTOR} = '-adaptor_seq '.$adaptor_seq;
}

#Get the adaptor from the sequence
sub get_adaptor {
	#TODO: Search for adaptor
	#return 'ATTT';
}


sub get_commands {
	my ($self) = @_;

	my %mapping = %{$self->{mapping}};

	# Fiddle with a local copy of commands
	my %local_commands = %COMMANDS;

	# Substitute configuration into commands
	for my $command_count (sort {$a<=>$b} keys %local_commands) {

		for (my $i = 0; $i < @{$local_commands{$command_count}{commands}}; $i++) {
			for my $replace_key (keys %mapping) {
		
				if ($local_commands{$command_count}{commands}->[$i] =~ /$replace_key/) {
					$local_commands{$command_count}{commands}->[$i] =~ s/$replace_key/$mapping{$replace_key}/g;
				}
			}
		}
	} 

	return \%local_commands;
}


return 1;





