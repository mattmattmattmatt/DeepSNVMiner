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


# Options from command line:

GetOptions(
			\%OPT, 
			"help|h",
	   		"man|m",
	   		"filename_stub=s",
	   		"read1_fastq=s",
	   		"read2_fastq=s",
	   		"start_command=s",
	   		"working_dir=s",
	   		"no_adaptor",
	   		"no_uid",
	   		"uid_len1=i",
	   		"uid_len2=i",
	   		"coord_bed=s",
	   		"ref_fasta=s",
	   		"bwa=s",
	   		"samtools=s",
	   		"sm_count=i",
	   		"sm_portion=f",
	   		"threads=i",
	   		"graph",
	   		"min_seqlen=i",
	   		"min_group=i",
	   		"conf_file=s",
	   		"no_pair_match",
	   		"cut_length=i",
	   		"no_revcom"
	   		) || modules::Exception->throw("ERROR: Problem with command line arguments");;

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{filename_stub} || !$OPT{read1_fastq} || !$OPT{read2_fastq} || !$OPT{coord_bed}) ;



=pod

=head1 SYNOPSIS

run_deepseq.pl -filename_stub unique_samplename -read1_fastq fastq1 -read2_fastq fastq2 -start_command command_to_start_from(default=first_command) -working_dir working_dir(default=pwd/filename_stub_RAND) -no_adaptor no_adaptor_sequence(default=present) -uid_len1 length_of_5prime_uid_length(default=10) -uid_len2 length_of_3prime_uid_length(default=0) -bwa full_bwa_path(default=/usr/bin/bwa) -samtools full_samtools_path(default=/usr/bin/samtools) -ref_fasta full_path_to_reference_fasta(must contain bwa index files as well) -coord_bed bed_file_containing_target_coordinates -sm_count min_number_of_reads_for_supermutant(default=10) -sm_portion min_fraction_of_variant_bases_within_UID_reads(default=0.90) -threads num_threads_for_bwa(default=1) -config create_config_file -graph generate_variant_graphs_for_each_genomic_region(requires R) -min_seqlen minimum_sequence_length_to_keep_seq(default=0) -min_group min_num_of_passing_groups_to_qualify_for_supermutant(default=2) -conf_file deepseq_conf_fiel(default=deepseq.conf) -no_pair_match don't_require_barcodes_in_read_pairs_to_match -cut_length remove_this_many_bases(default=uid1+uid2) -no_revcom no_revcom_for_read2  

Required flags: -filename_stub -read1_fastq -read2_fastq -coord_bed

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

run_deepseq.pl -> Wrapper script for deepseq run

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

sample -f file

=cut


my $bwa;
my $bwa_index;
my $samtools;
my $ref_fasta;
my $PRINT_TO_LOGFILE = 1;
my $PRINT_TO_STDOUT = 1;
my (undef,$base) = fileparse(abs_path($0));
my $conf_file = defined $OPT{conf_file}?$OPT{conf_file}:$base.'/deepseq.conf';

#Set the variables either using a conf file or parse them from the command line
if (-e $conf_file) {
 	#Read the config file to get paths needed
	open(FILE,"$conf_file") || modules::Exception->throw("Can't open file $conf_file\n");
	
	while (<FILE>) {
		chomp;
		if (/samtools=(.+)$/) {
			$samtools = $1;
			if ( !-e $samtools ) {
				modules::Exception->throw("Samtools $samtools doesn't exist");	
			}
		} elsif (/bwa=(.+)$/) {
			$bwa = $1;
			if ( !-e $bwa ) {
				modules::Exception->throw("BWA $bwa doesn't exist");	
			}
		} elsif (/ref_fasta=(.+)$/) {
			$ref_fasta = $1;
			if ( !-e $ref_fasta ) {
				modules::Exception->throw("Ref fasta $ref_fasta doesn't exist");	
			}
		} elsif (/bwa_index=(.+)$/) {
			$bwa_index = $1;
			if ( !-e $bwa_index ) {
				modules::Exception->throw("BWA index file $bwa_index doesn't exist");	
			}
		} 
	}
} else {
	#otherwise all the variables need to be defined on the command line or available in /usr/bin
	if (!$OPT{ref_fasta}) {
		modules::Exception->throw("ERROR: Without config file must use -ref_fasta");
	} else {
		$ref_fasta = $OPT{ref_fasta};		
	}
	
	my $ref_fasta_abs = abs_path($ref_fasta);
	
	#Check the bwa index is present....
	my $bwa_index_file = $ref_fasta_abs . '.sa';
	if (!-e $bwa_index_file) {
		#Sometimes it ref.fa.sa but sometimes it's ref.sa so check both
		print "Can't find the index file $bwa_index_file\n";
		($bwa_index_file = $ref_fasta_abs) =~ s/\.fa$/\.sa/;
		print "Checking for index file $bwa_index_file\n";
		if (!-e $bwa_index_file) {
			modules::Exception->throw("ERROR: Need to generate the index files for bwa in the referece directory; either 1) run bwa beforehand with 'bwa index -a bwtsw ref.fa' OR 2) run configure_deepseq.pl first");
		}
	}
	#bwa only requires the file prefix name
	($bwa_index = $bwa_index_file) =~ s/\.sa$//;
	
	if ($OPT{bwa}) {
		$bwa = $OPT{bwa};		
	} elsif (-e '/usr/bin/bwa') {
		$bwa = '/usr/bin/bwa';			
	} else {
		modules::Exception->throw("ERROR: Without config file must use -bwa or have bwa at /usr/bin/bwa");
	}
	
	if ($OPT{samtools}) {
		$samtools = $OPT{samtools};		
	} elsif (-e '/usr/bin/samtools') {
		$samtools = '/usr/bin/samtools';			
	} else {
		modules::Exception->throw("ERROR: Without config file must use -samtools or have samtools at /usr/bin/samtools");
	}
}

#Check the fasta is available
if (!-e $ref_fasta) {
	modules::Exception->throw("ERROR: Problem with $ref_fasta");
}

my $sys_call = modules::SystemCall->new();

my $ref_fasta_abs = abs_path($ref_fasta);


my $pwd = `pwd`;
chomp $pwd;

#now get the sequence specific variables
my $filename_stub = $OPT{filename_stub};
my $read1_fastq = $OPT{read1_fastq};
my $read2_fastq = $OPT{read2_fastq};

if ($read1_fastq =~ /bz2$/ || $read1_fastq =~ /gz$/ || $read1_fastq =~ /zip$/) {
	modules::Exception->throw("ERROR: Can't work with compressed files; please uncompress first");
}

my $working_dir;
if (defined $OPT{working_dir}) {
	$working_dir = $OPT{working_dir};
} else {
	$working_dir = $pwd . '/' . $filename_stub . '_' . $$;
	if (! -e $working_dir){
		`mkdir $working_dir`;
    }
}

if ($OPT{graph}) {
	if (!-e "$working_dir/graphs") {
		`mkdir $working_dir/graphs`;
	}
}

my $adaptor = defined $OPT{no_adaptor}?0:1; #default is that there is adaptor sequence info
my $uid = defined $OPT{no_uid}?0:1; #default is that there is uid info
my $uid_len1 = defined $OPT{uid_len1}?$OPT{uid_len1}:10;
my $uid_len2 = defined $OPT{uid_len2}?$OPT{uid_len2}:0;
my $coord_bed = $OPT{coord_bed};
if ( !-e $coord_bed ) {
	modules::Exception->throw("File $coord_bed doesn't exist");	
}

#Change all directories files to absolute path
my $working_dir_abs = abs_path($working_dir);
my $bwa_abs = abs_path($bwa);
my $samtools_abs = abs_path($samtools);
my $coord_bed_abs = abs_path($coord_bed);
my $read1_fastq_abs = abs_path($read1_fastq);
my $read2_fastq_abs = abs_path($read2_fastq);

#Check they all exist
if ( !-e $ref_fasta_abs ) {
	modules::Exception->throw("File $ref_fasta_abs doesn't exist");	
}

if ( !-e $read1_fastq_abs || !-e $read2_fastq_abs ) {
	modules::Exception->throw("File $read1_fastq_abs || !-e $read2_fastq_abs doesn't exist");	
}

if ( !-e $working_dir_abs ) {
	modules::Exception->throw("File $working_dir_abs doesn't exist");	
}

if ( !-e $bwa_abs ) {
	modules::Exception->throw("File $bwa_abs doesn't exist");	
}

if ( !-e $samtools_abs ) {
	modules::Exception->throw("File $samtools_abs doesn't exist");	
}

if ( !-e $coord_bed_abs ) {
	modules::Exception->throw("File $coord_bed_abs doesn't exist");	
}


my $read1_fastq_short = my $read2_fastq_short;
unless (chdir $working_dir_abs) {
    modules::Exception->throw("Didn't manage to change to working dir [$working_dir_abs].  Still in this directory [$pwd]");
    #unless $pwd eq $working_dir;
} else {
	#Create symlinks for read files (don't copy; too time consuming)
	`ln -s $read1_fastq_abs .`;
	`ln -s $read2_fastq_abs .`;
	#pass these as we've created symlinks now
	($read1_fastq_short) = basename($read1_fastq_abs);
	($read2_fastq_short) = basename($read2_fastq_abs);
}




#Set the required arguments
my %args = (
			-filename_stub => $filename_stub,
	   		-read1_fastq => $read1_fastq_short,
	   		-read2_fastq => $read2_fastq_short,
	   		-working_dir => $working_dir_abs,
	   		-uid_len1 => $uid_len1,
	  		-uid_len2 => $uid_len2,
	   		-coord_bed => $coord_bed_abs,
	   		-ref_fasta => $ref_fasta_abs,
	   		-bwa => $bwa_abs,
	   		-bwa_index => $bwa_index,
	   		-samtools => $samtools_abs,
	   		-base=>$base
			);

#Optional arguments
if (!$uid) {
	$args{-no_uid} = 1;
}

if ($OPT{no_pair_match}) {
	$args{-no_pair_match} = 1;
}

if ($OPT{no_revcom}) {
	$args{-no_revcom} = 1;
}

if ($OPT{cut_length}) {
	$args{-cut_length} = $OPT{cut_length};
}

if (defined $OPT{sm_count}) {
	$args{-sm_count} = $OPT{sm_count};
}

if (defined $OPT{sm_portion}) {
	$args{-sm_portion} = $OPT{sm_portion};
}

if (defined $OPT{threads}) {
	$args{-num_threads} = $OPT{threads};
}

if (defined $OPT{min_seqlen}) {
	$args{-min_seqlen} = $OPT{min_seqlen};
}

if (defined $OPT{min_group}) {
	$args{-min_group} = $OPT{min_group};
}

my $start_command = defined $OPT{'start_command'}?$OPT{start_command}:'check_fastq';
my $start_number = 1;


# Start a log/results file 

my $logfile = $filename_stub . '.log'; 
my $OUT;
if (! -e $logfile){
	#$logfile = $working_dir_abs . '/' . $filename_stub . '.log';
	open($OUT, ">$logfile")
	    or modules::Exception->throw("Unable to start logfile [$logfile] : $!\n");
} else {
	open($OUT, ">>$logfile")
	    or modules::Exception->throw("Unable to start logfile [$logfile] : $!\n");
}


### select $OUT; # write as much blurb to logfile instead of STDOUT (turn this back to STDOUT with 'select STDOUT;')

print $OUT "\n### Log edited : " . &low_rent_timestamp() . " ###\n\n" if $PRINT_TO_LOGFILE;

print STDOUT "Logfile (check here also for error messages) : $logfile\n" if $PRINT_TO_STDOUT;
print $OUT "Working directory : $working_dir\n" if $PRINT_TO_LOGFILE;

# Create deepseq object
my $deepseq = modules::Deepseq->new(
									%args
									);

if ($start_command ne 'check_fastq') {
	if (!$deepseq->check_step($start_command)) {
		modules::Exception->throw("ERROR: Can't start at step $start_command; options are check_fastq,pool_reads,bwa,call_variants,report,or graph");
	} 
	$start_number = $deepseq->get_step_number($start_command);
}


#If there is adaptor sequence(s) get it from the fastq and mark it in each sequence to avoid creating supergroup of adaptors
if ($adaptor) {
	if (my $adaptor_seqs = $deepseq->get_adaptor()) {
		my $adaptor_str = join(",",@{$adaptor_seqs});
		print "Will remove the following adaptors:\n$adaptor_str\n";
		$deepseq->set_adaptor($adaptor_seqs);
	}	
}

# Retrieve commands from config file
my $commands = $deepseq->get_commands();

# Try to run each command in order
for my $command_count (sort keys %{$commands}) { 
	my @commands = @{$commands->{$command_count}{commands}};
	my @outs = @{$commands->{$command_count}{stdout}};
	
	if (@commands != @outs) {
		modules::Exception->throw("ERROR: Different number of commands and stdout flags");
	}
	
	my $block = $commands->{$command_count}{block};
	if ($start_number > 1 && $command_count < $start_number) {
		next;
	}
	
	if ($block eq 'graph') {
		next unless $OPT{graph};
	}
	
	print "\n\nRunning block $block\n\n";
	
	for (my $i = 0; $i < @commands; $i++) {
		print $OUT "\tRun " . $commands[$i] . "\n" if $PRINT_TO_LOGFILE;
		if ($outs[$i] == 1) {
			my $output = `$commands[$i]`;
			print $OUT "$output" if $PRINT_TO_LOGFILE;			
		} else {
			$sys_call->run($commands[$i]);
		}
	}
	
}


# Close logfile
close($OUT);

# Return to previous directory
chdir $pwd;


sub low_rent_timestamp() {

    my ($sec,
	$min,
	$hour,
	$mday,
	$month,
	$year,
	$wday,
	$yday,
	$isdst) = localtime(time);

    $year += 1900;
    $month += 1;

    my $timestamp_str = $year  . '-' . $month . '-' . $mday . ' ' .  $hour . ':' . $min . ':' . $sec;

    return $timestamp_str;
}
