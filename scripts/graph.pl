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
	   "coord_bed=s",
	   "super_mutants=s",
	   "var_dir=s",
	   "r_dir=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{coord_bed} || !$OPT{super_mutants} || !$OPT{var_dir} || !$OPT{r_dir});



=pod

=head1 SYNOPSIS

graph.pl -coord_bed coord_bed_file -super_mutants super_mutants_file -var_dir variants_sub_directory -r_dir R_directory

Required flags: -coord_bed -super_mutants -var_dir -r_dir

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

graph.pl -> Script to generate graphs of supermutants for every range in the bed file

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

sample -f file

=cut

#while(<>) {@cols = split /\s+/; $coord_tally{$cols[0]}++;} foreach $coord ( (10..20) ) {print join ("\t", $coord, defined($coord_tally{$coord})?$coord_tally{$coord}:0) ."\n"; }' FILENAMESTUB.passing_supermutants.txt > supermutant_coords.txt; R CMD BATCH R_DIR/supermut_coord.R; ls supermut_coords.graph.pdf) # Distribution of supermutants graph



my $coord_bed = $OPT{coord_bed};
if ( !-e $coord_bed ) {
	modules::Exception->throw("File $coord_bed doesn't exist");	
}

my $super_mutants = $OPT{super_mutants};
if ( !-e $super_mutants ) {
	modules::Exception->throw("File $super_mutants doesn't exist");	
}

my $r_dir = $OPT{r_dir};
if ( !-d $r_dir ) {
	mkdir($r_dir);
}

my (undef,$workdir) = fileparse($super_mutants);
if ( !-d $workdir ) {
	modules::Exception->throw("File $workdir doesn't exist");	
}

my %var_tally = ();
my %sm_tally = ();

open(SUPER,"$super_mutants") || modules::Exception->throw("Can't open file $super_mutants\n");
while (<SUPER>) {
    chomp;
    next unless /\S/;
    my ($chr,$coord) = split;
    $sm_tally{"$chr:$coord"}++;
}

my $var_dir = $OPT{var_dir};
opendir(DIR,"$var_dir") || modules::Exception->throw("Can't open directory $var_dir");
my @variant_files = grep {/variants$/} readdir DIR;

for my $var_file (@variant_files) {
	open(FILE,"$var_dir/$var_file") || modules::Exception->throw("Can't open file $var_file\n");
	while (<FILE>) {
		chomp;
    	next unless /\S/;
    	my ($chr,$coord) = split("\t");
    	$var_tally{"$chr:$coord"}++;
	}
	close FILE;
}

my @r_commands = ();

if ($coord_bed =~ /ref.bed/) {
	$coord_bed =~ s/ref/ref_graph/;	
}

print "Open bed $coord_bed\n";

open(BED,"$coord_bed") || modules::Exception->throw("Can't open file $coord_bed\n");
while (<BED>) {
	my ($chr,$start,$end) = split;	
	next unless $chr =~ /^chr/ || $chr =~ /^[0-9XYM]/; #skip headers
	if ($end - $start > 1000) {
		#Skip graphs larger than 1000bp
		print "Skip graph $chr:$start-$end; too large a genomic region\n";
	}
	#Generate the R commands for each block
	my $sm_out = "$chr:$start"."_supermutant_coords.txt";
	my $var_out = "$chr:$start"."_variants_coords.txt";
	open(SMOUT,">$sm_out") || modules::Exception->throw("Can't open file to write $sm_out\n");
	open(VAROUT,">$var_out") || modules::Exception->throw("Can't open file to write $var_out");
	my $sm_flag = 0;
	my $var_flag = 0;
	for my $coord ($start..$end) {
		if (exists $sm_tally{"$chr:$coord"}) {
			$sm_flag = 1;
		}
		if (exists $var_tally{"$chr:$coord"}) {
			$var_flag = 1;
		}
		
		print SMOUT join ("\t", $coord, defined($sm_tally{"$chr:$coord"})?$sm_tally{"$chr:$coord"}:0) . "\n";
		print VAROUT join ("\t", $coord, defined($var_tally{"$chr:$coord"})?$var_tally{"$chr:$coord"}:0) . "\n";
	}
	close SMOUT;
	close VAROUT;
	#Rscript R_DIR/group_size.R WORKING_DIR WORKING_DIR/FILENAMESTUB_R1.group_counts_final
	my $sm_jpeg = "$chr:$start"."_supermutant.jpeg";
	my $var_jpeg = "$chr:$start"."_allvariants.jpeg";
	#Now create the R commands
	if ($sm_flag) {
		my $r_super = join(" ", 'Rscript', $r_dir .'/supermut_coord.R',$workdir,$sm_out,$sm_jpeg); #don't bother with empty files
		push @r_commands, $r_super;	
	}
	if ($var_flag) {
		my $r_variants = join(" ", 'Rscript', $r_dir .'/mut_coord.R',$workdir,$var_out,$var_jpeg);
		push @r_commands, $r_variants;
	}
	
}

my $sys_call = modules::SystemCall->new();

for my $command ( @r_commands ) {
	$sys_call->run($command);
}
#$sys_call->run("rm *_coords.txt");
$sys_call->run("mv *jpeg graphs");
