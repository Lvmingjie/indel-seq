#! usr/bin/perl -w
use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2019.12.16	Counting FASTQ data amount	
# Modified:			2022.11.17	Counting FASTA 
# Modified:			2023.02.17	Counting AnnoStat_tables
############################################################

my $in = shift;
my $out = shift;

die "
Usage:	
perl   xxx.pl   file_in   out_prefix

INPUT:
indels_out_2023Apr181120_table.xls
t1_stack_cols_Rtable

" if !defined $out;

my $date = `date +%Y%b%d%H%M`;
$date =~ s/[\r\n]+//g;

my($flie_name, $tag, $count, @head, @t, $i);
my(%all_gene_number, $gene_tag, %count_gene_class, %count_gene_ori);
my($file_number, $out_number);
############################################################


## Input
open (IN, "< $in") or die $!;
open (OUT, "> $out"."_".$date.".xls") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t = split ("\t", $_);
	
	for $i(7..$#t){
		if($t[$i] =~ /\//){
			$t[$i] =~ s/\//_/;
			print OUT $t[0]."_".$t[1]."\t".$t[$i]."\t".$t[$i+1]."\n";
		}
	}
}
close IN;
close OUT;
print "\n";
