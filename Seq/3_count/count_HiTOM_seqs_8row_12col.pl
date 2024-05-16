#! usr/bin/perl -w
#use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2020.04.23	
# Modified:			2022.09.24	Hi-TOM
############################################################


my $merged_fq = shift;
my $primer_table = shift;
my $primer_ID = shift;
my $barcode_row = shift;
my $barcode_col = shift;
my $out = shift;

die  "
Usage:
perl  xxx.pl  merged.fq.gz  primer_table  primer_ID  barcode_row  barcode_col  out_prefix

Input:
indel1_1R_merged.fq.gz
primer_table
Indel1
barcode_row
barcode_col
indel1_1R_out



Output:
indel1_1R_out_xxxxxx_.xls
indel1_1R_out_xxxxxx_top.fa

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

my(%list, %head1, %head2, %check_duplicates, $rc_head1, $rc_head2);
my($count, @t, $tail1, $tail2, $i, $j, $seq_head, $seq_tail);
my(%total_count, $rc_seq, $order, $percent, $top1_percent);

my(%forward, %reverse, $forward_rc, $reverse_rc);
my(%barcode_row_group, %barcode_col_group, $rc_tag, %col_and_row);
my($x, $y, $m, $n);
############################################################


## Record primer_table
open (IN, "< $primer_table") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t = split("\t", $_);
	
	#if($t[0] !~ /^#/){
	$t[0] =~ s/_//g;
	if($t[0] eq $primer_ID){	##只考虑一套引物	
		$list{$t[0]} ++;
		$forward{$t[0]} = $t[2];
		$reverse{$t[0]} = $t[4];
		
		#$check_duplicates{$t[2]} ++;
		#$check_duplicates{$t[4]} ++;
		$forward_rc = base_reverse ($t[2]);
		#$check_duplicates{$forward_rc}++;
		$reverse_rc = base_reverse ($t[4]);
		#$check_duplicates{$reverse_rc}++;
		print $_."\n";
	}
}
close IN;


#foreach $i (keys %check_duplicates){
#	if($check_duplicates{$i} > 1){
#		print "WARNING! Repeat Primers!\t".$i."\n";
#	}
#}


## Record barcode_row
open (IN, "< $barcode_row") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t = split("\t", $_);
	
	if($t[1]){
		$barcode_row_group{$t[0]} = uc($t[1]);
		print "barcode_row:\t".$t[0]."\t".$t[1]."\n";
	}
}
close IN;

## Record barcode_col
open (IN, "< $barcode_col") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t = split("\t", $_);
	
	if($t[1]){
		$barcode_col_group{$t[0]} = uc($t[1]);
		print "barcode_col:\t".$t[0]."\t".$t[1]."\n";
	}
}
close IN;


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
			
			##8 row 与 reverse primer 绑定, 12 col 与 forward primer绑定;
			if($t[$i] =~ /([ATCG]{5})$forward{$j}(\w+)$reverse_rc([ATCG]{5})/){			
				$total_count{$j}++;
				
				$rc_tag = base_reverse ($3);
				$col_and_row{$1."_".$rc_tag} ++;
				#print $1."\t".$rc_tag."\t".$col_and_row{$1."_".$rc_tag}."\n";
				
				${"count_".$1."_".$rc_tag}{$2} ++;
			}elsif($t[$i] =~ /([ATCG]{5})$reverse{$j}(\w+)$forward_rc([ATCG]{5})/){
				$total_count{$j}++;
				
				$rc_tag = base_reverse ($3);
				$col_and_row{$rc_tag."_".$1} ++;
				#print $rc_tag."\t".$1."\t".$col_and_row{$rc_tag."_".$1}."\n";
				
				$rc_seq = base_reverse ($2);
				${"count_".$rc_tag."_".$1}{$rc_seq} ++;
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
#print OUT "#\tForward\tReverse\tTotal_Number\tPercent\tTop1_seq\tTop1_Num\tPercent\tTop2_seq\tTop2_Num\tPercent\tTop3_seq\tTop3_Num\tPercent\n";
print OUT "#Row_Col\tTotal_Number\tPercent\tTop1_seq\tTop1_Num\tPercent\tTop2_seq\tTop2_Num\tPercent\tTop3_seq\tTop3_Num\tPercent\n";
foreach $i (sort keys %list){
	
	if ($total_count{$i} && $total_count{$i} >= 100){
		#print OUT $i."\t".$forward{$i}."\t".$reverse{$i}."\t".$total_count{$i}."\t100%";
		
		foreach $m(sort keys %barcode_row_group){
			foreach $n(sort {$a <=> $b} keys %barcode_col_group){
				print OUT $m."_".$n;
				
				$x = $barcode_row_group{$m};
				$y = $barcode_col_group{$n};
				
				if(defined $col_and_row{$y."_".$x}){
					print OUT "\t".$col_and_row{$y."_".$x}."\t100%";
					
					undef $order;
					foreach $j (sort {${"count_".$y."_".$x}{$b} <=> ${"count_".$y."_".$x}{$a}} keys %{"count_".$y."_".$x}){
						$order ++;
							
						undef $percent;
						$percent = int((${"count_".$y."_".$x}{$j} / $col_and_row{$y."_".$x})*10000) /100;
							
						if($order == 1){
							print OUT "\t".$j."\t".${"count_".$y."_".$x}{$j}."\t".$percent."%";
							print OUT2 ">".$i."_".$m."_".$n."_top".$order."_".$percent."%\n".$j."\n";
							#$top1_percent = $percent;
						}elsif($order > 1){
							if($percent >= 5){
								print OUT "\t".$j."\t".${"count_".$y."_".$x}{$j}."\t".$percent."%";
								print OUT2 ">".$i."_".$m."_".$n."_top".$order."_".$percent."%\n".$j."\n";
							}
							#if( ($percent/$top1_percent) >= 0.75){
							#	print OUT2 ">".$i."_top".$order."_".$percent."%\n".$j."\n";	
							#}
						}
					}
				}
					print OUT "\n";
			}
		}
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
