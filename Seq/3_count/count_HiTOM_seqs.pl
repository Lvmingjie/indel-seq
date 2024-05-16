#! usr/bin/perl -w
#use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2020.04.23	
# Modified:			2022.09.24	Hi-TOM
############################################################


my $merged_fq = shift;
my $primer_table = shift;
my $out = shift;

die  "
Usage:
perl   xxx.pl   merged.fq.gz   primer_table   out_prefix

Input:
all_merged_12samples.fq.gz
primer_table_jgy20230928
all_merged_12samples_out

Output:
all_merged_12samples_out_xxxxxx_.xls
all_merged_12samples_out_xxxxxx_top.fa

" if !defined $out;

my %reverse_comp = ( 
"A" => "T",
"T" => "A",
"C" => "G",
"G" => "C",
"N" => "N",
"a" => "t",
"t" => "a",
"c" => "g",
"g" => "c",
"n" => "n"
);

my $date = `date +%Y%b%d%H%M`;
$date =~ s/[\r\n]+//g;

my (%list, %head1, %head2, %check_duplicates, $rc_head1, $rc_head2);
my ($count, @t, $tail1, $tail2, $i, $j, $seq_head, $seq_tail);
my (%total_count, $rc_seq, $order, $percent, $top1_percent);

my (%forward, %reverse, $forward_rc, $reverse_rc);
############################################################


## Record primer_table
open (IN, "< $primer_table") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t = split("\t", $_);
	
	if($t[0] !~ /^#/){
		
		$list{$t[0]} ++;
		$forward{$t[0]} = $t[2];
		$reverse{$t[0]} = $t[4];
		
		$check_duplicates{$t[2]} ++;
		$check_duplicates{$t[4]} ++;
		$forward_rc = base_reverse ($t[2]);
		$check_duplicates{$forward_rc}++;
		$reverse_rc = base_reverse ($t[4]);
		$check_duplicates{$reverse_rc}++;
	}
}
close IN;


foreach $i (keys %check_duplicates){
	if($check_duplicates{$i} > 1){
		print "WARNING! Repeat Primers!\t".$i."\n";
	}
}


## Count top3 abundant target seqs
if($merged_fq =~ /\.fq\.gz$/ || $merged_fq =~ /\.fastq\.gz$/){		## ".fq.gz"	or	".fastq.gz"
	open (IN, "gunzip -c $merged_fq |") or die $!;
}elsif($merged_fq =~ /\.fq$/ || $merged_fq =~ /\.fastq$/){			## ".fq"  	or	".fastq"
	open (IN, "< $merged_fq") or die $!;
}else{
	die "Check the fastq file name!\t".$merged_fq."\n"; 
}
while (<IN>){
    $_ =~ s/[\r\n]+//g;
	
	if(/^@\w+:\d+:\w+/){
		$i = 0;
	}else{
		$i ++;
	}
	$t[$i] = $_;
	
	if($i == 1){
		$count ++;
        print "Analysing merged reads: $count\r" if $count % 1000 == 0;
		
		foreach $j (sort keys %list){
			$forward_rc = base_reverse ($forward{$j});
			$reverse_rc = base_reverse ($reverse{$j});
			
			if($t[$i] =~ /$forward{$j}(\w+)$reverse_rc/){
				$total_count{$j}++;
				${"count_".$j}{$1} ++;
			}elsif($t[$i] =~ /$reverse{$j}(\w+)$forward_rc/){
				$total_count{$j}++;
				$rc_seq = base_reverse ($1);
				${"count_".$j}{$rc_seq} ++;
			}
		}
	}
}
close IN;
print "\n";
undef $count;


## Output
open (OUT, "> $out"."_".$date.".xls") or die $!;
open (OUT2, "> $out"."_".$date."_top.fa") or die $!;
print OUT "#\tForward\tReverse\tTotal_Number\tPercent\tTop1_seq\tTop1_Num\tPercent\tTop2_seq\tTop2_Num\tPercent\tTop3_seq\tTop3_Num\tPercent\n";
foreach $i (sort keys %list){
	
	if ($total_count{$i} && $total_count{$i} >= 1000){ ##设置基因型最低频次阈值, >= 1000 的基因型全部输出;
	#if ($total_count{$i} && $total_count{$i} >= 100){ ##设置基因型最低频次阈值, >= 1000 的基因型全部输出;
		print OUT $i."\t".$forward{$i}."\t".$reverse{$i}."\t".$total_count{$i}."\t100%";
		
		undef $order;
		foreach $j (sort {${"count_".$i}{$b} <=> ${"count_".$i}{$a}} keys %{"count_".$i}){
			$order ++;
			
			undef $percent;
			$percent = int((${"count_".$i}{$j} / $total_count{$i})*10000) /100;
			
			if($order == 1){
				print OUT "\t".$j."\t".${"count_".$i}{$j}."\t".$percent."%";
				print OUT2 ">".$i."_top".$order."_".$percent."%\n".$j."\n";
				#$top1_percent = $percent;
			}elsif($order > 1){
				#if($percent >= 1){	##设置基因型最低频次阈值, >= 1% 的基因型全部输出;
				if($percent >= 0.5){	##设置基因型最低频次阈值, >= 0.5% 的基因型全部输出;
					print OUT "\t".$j."\t".${"count_".$i}{$j}."\t".$percent."%";
					print OUT2 ">".$i."_top".$order."_".$percent."%\n".$j."\n";	
				}
				#if( ($percent/$top1_percent) >= 0.75){
				#	print OUT2 ">".$i."_top".$order."_".$percent."%\n".$j."\n";	
				#}
			}
		}
		print OUT "\n";
	}
}
close OUT;
close OUT2;
print "\n";


sub base_reverse{
	my @bre = reverse(split ('', $_[0]));
	my $baseout;
	foreach $i(@bre){
	     $baseout .= $reverse_comp{$i};
	}
	$baseout;
}
