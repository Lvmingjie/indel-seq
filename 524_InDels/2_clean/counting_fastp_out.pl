#! usr/bin/perl -w
use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2019.05.05
# Modified: 		2020.01.26
# Modified:
############################################################

my $list = shift;
#my $out = shift;

die "
perl   xxx.pl   fastp_log_all

Input:
fastp_log
Output:
fastp_log_countout_xxxxxxxx.xls

" if !defined $list;


my($tag, @t, $i, $j, $file_name, $k);
my(%sample_group, %target_group);
my(%total_reads, %uniq_reads, %uniq_number);
my($label, %label_group, %value, $num);
my($judge, $sample);
############################################################

`date` =~ /\S+\s+(\S+)\s+(\d+)\s+(\d+):(\d+):\d+\s+\S+\s+(\d+)/;
my $date = $5.$1.$2.$3.$4;

#`date` =~ /^(20\d+).+\s+(\d+).+\s+(\d+).+\s+\S+\s+(\d+):(\d+):\d+\s+\S+/;
#my $date = $1.$2.$3.$4.$5;
############################################################

open (OUT, "> $list"."_countout_".$date.".xls") or die $!;
#print OUT "#Sample";

open (IN, "< $list") or die $!;
while(<IN>){
    $_ =~ s/[\r\n]+//g;
	
	if($_ !~ /,/ && $_ =~ /(.+):(\s+(\S+))?/ && $_ !~ /ERROR/){
		$label = $1;
		#$num = $2;
		$label =~ s/ /_/g;
		
		print OUT $label."\t" if !defined $judge;
		
		#if($num){
		#	print OUT "\t".$num;
		#}else{
		#	print OUT "\t-";
		#}
		
	}elsif(/fastp\s+-i\s+/){
		#$sample = $1;
		#print OUT "\t".$sample;
		
		print OUT "#Sample\n" if !defined $judge;
		$judge ++;
		
	#}elsif(/fastp\s+v\d+.+time\s+used:.+seconds/){
		
	}
}
close IN;


open (IN, "< $list") or die $!;
while(<IN>){
    $_ =~ s/[\r\n]+//g;
	
	if($_ !~ /,/ && $_ =~ /(.+):(\s+(\S+))?/ && $_ !~ /ERROR/){
		$label = $1;
		$num = $2;
		$label =~ s/ /_/g;
		
		if($num){
			print OUT $num."\t";
		}else{
			print OUT "-\t";
		}
		
	}elsif(/fastp\s+-i\s+(\S+)_R1\.fq\.gz\s+/ || /fastp\s+-i\s+.+\/(\S+)_1\.fastq\.gz\s+/){
		$sample = $1;
		print OUT $sample."\n";
		
	#}elsif(/fastp\s+v\d+.+time\s+used:.+seconds/){
		
	}

}
close IN;
close OUT;

