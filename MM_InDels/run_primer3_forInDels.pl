#! usr/bin/perl -w
use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2020.05.03	运行primer3; 
#					primer3参数设置：http://primer3.ut.ee/primer3web_help.htm
# Modified:			2023.04.15	For InDel context

############################################################

my $fasta_file = shift;
my $out = shift;

die  "
Usage:
perl   xxx.pl   fasta_file   output_fullname

Note:
indels_out_member3_2023Apr151519_context.fa
indels_out_member3_2023Apr151519_context_primer3out.xls

" if !defined $out;

my $date = `date +%Y%b%d%H%M`;
$date =~ s/[\r\n]+//g;

my(@t, $count, $tag, $i, $seq, $SRR_length, $file_name);
my($primer_left_0_seq, $primer_right_0_seq);
my($primer_left_0, $primer_right_0, $primer_pair_0_product_size);
############################################################


## record fasta_file
if($fasta_file =~ /\.gz$/){
	open (IN, "gunzip -c $fasta_file |") or die $!;
}else{
	open (IN, "< $fasta_file") or die $!;
}
open (OUT2, "> $out") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	if($_ =~ /^>(\S+)/){
		$tag = $1;
		print "Dealing with: $tag\r";
		
		undef $SRR_length;
		if($tag =~/\d+_context_\d+_\d+_(\d+)bp/){
			$SRR_length = $1 -200;
		}
	}else{

		$seq = $_;
		
		#$count ++;
		open (OUT, "> temp_primer3_input_".$date) or die $!;
		print OUT "SEQUENCE_ID=".$tag."\n";
		print OUT "SEQUENCE_TEMPLATE=".$seq."\n";
		print OUT "SEQUENCE_TARGET=101,".$SRR_length."\n";
		print OUT "SEQUENCE_EXCLUDED_REGION=101,".$SRR_length."\n";
		print OUT "PRIMER_TASK=generic"."\n";
		print OUT "PRIMER_PICK_LEFT_PRIMER=1"."\n";
		print OUT "PRIMER_PICK_INTERNAL_OLIGO=0"."\n";
		print OUT "PRIMER_PICK_RIGHT_PRIMER=1"."\n";
		print OUT "PRIMER_OPT_SIZE=20"."\n";
		print OUT "PRIMER_MIN_SIZE=18"."\n";
		print OUT "PRIMER_MAX_SIZE=22"."\n";
		print OUT "PRIMER_EXPLAIN_FLAG=1"."\n";
		#print OUT "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/mnt/new15T/1_biosofts/primer3-2.4.0/src/primer3_config/"."\n=\n";
		#print OUT "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/mnt/freeNAS_50T/1_biosofts/primer3-2.4.0/src/primer3_config/"."\n=\n";
		print OUT "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/aglab200T/aglab/cr/1_biosofts/primer3-2.4.0/src/primer3_config/"."\n=\n";
		
		close OUT;
		
		$file_name = "temp_primer3_input_".$date;
		open (IN2, "cat $file_name | primer3_core |") or die $!;
		while(<IN2>){
			$_ =~ s/[\r\n]+//g;
			$primer_left_0_seq			= $1 if $_ =~ /PRIMER_LEFT_0_SEQUENCE=(\w+)/;
			$primer_right_0_seq			= $1 if $_ =~ /PRIMER_RIGHT_0_SEQUENCE=(\w+)/;
			$primer_left_0 				= $1 if $_ =~ /PRIMER_LEFT_0=(\d+,\d+)/;
			$primer_right_0				= $1 if $_ =~ /PRIMER_RIGHT_0=(\d+,\d+)/;
			$primer_pair_0_product_size	= $1 if $_ =~ /PRIMER_PAIR_0_PRODUCT_SIZE=(\d+)/;
		}
		close IN2;
		
		$primer_left_0_seq = "0" 			if !defined $primer_left_0_seq;
		$primer_right_0_seq = "0" 			if !defined $primer_right_0_seq;
		$primer_left_0 = "0" 				if !defined $primer_left_0;
		$primer_right_0 = "0" 				if !defined $primer_right_0;
		$primer_pair_0_product_size = "0" 	if !defined $primer_pair_0_product_size;

		print OUT2 $tag."\t".$seq."\t".$primer_left_0_seq."\t".$primer_left_0."\t".$primer_right_0_seq."\t".$primer_right_0."\t".$primer_pair_0_product_size."\n";
		
	}
}
close IN;
close OUT2;
print "\n";

