#! usr/bin/perl -w
use strict;

# Author:			Rui Chen <chenrui_taas@126.com; chenrui.taas@gmail.com>
# Program Date:		2021.02.07	For Indel marker 
# Modified:			2021.06.23	提取两个群体之间差异的InDels, geno文件貌似有问题
# Modified:			2023.04.14	filtering grapre high polymorphism InDel markers
# Modified:
############################################################


my $vcf_InDels		= shift;
my $genome_fa 		= shift;
my $out				= shift;

die "
perl   xxx.pl   vcf_InDels   genome.fa   out_prefix

Input:
final_merged_InDels_grape524_2021Aug08.vcf.gz
grape_12X_genome_ChrID.fa
indels_out

Output:
indels_out_xxxxxx.vcf
indels_out_xxxxxx_table.xls
indels_out_xxxxxx_context.fa

" if !defined $out;

my $date = `date +%Y%b%d%H%M`;
$date =~ s/[\r\n]+//g;


my($count, $count1, $count2, %sample_all, %popu1, %popu2);
my($tag, %name, %seq);
my(@t, $i, %record_gene); 
my($chr, $site, $ref, $target, $var, @m);
my($left, $right, $length, $length_ref, $length_var, $diff);
my($total_1, $total_2, $max_length);
my(@n, %p1, %p2, @x, @y, $tag1, $tag2);
my(@mark1, @mark2, @w, $p1_out, $p2_out);
my(%record, $end, $member, $min_length, $diff_length, %indel_length, $length_banch, $judge, $seq_length, $seq_context);
my($longest, $genotype_cluster, %genotype, %genotype2, $genotype_missing, $genotype_max, $times);
############################################################


## Record SNPs
#if ($vcf_SNPs =~ /\.vcf\.gz$/ ){			
#	open (IN, "gunzip -c $vcf_SNPs |") or die $!;
#}elsif($vcf_SNPs =~ /\.vcf$/ 	){
#	open (IN, "< $vcf_SNPs") or die $!;
#}
#while (<IN>) {
#	$_ =~ s/[\r\n]+//g;
#	$count ++;
#	print "Counting lines:\t$count\r" if $count%100000 == 0;
	
#	unless(/^#/){
#		undef @t;
#		@t = split("\t", $_);
		
#		$record{$t[0]."_".$t[1]} ++;
#	}
#}
#close IN;
#print "\n";
#undef $count;


## Record InDels
#if ($vcf_InDels =~ /\.vcf\.gz$/ ){			
#	open (IN, "gunzip -c $vcf_InDels |") or die $!;
#}elsif($vcf_InDels =~ /\.vcf$/ 	){
#	open (IN, "< $vcf_InDels") or die $!;
#}
#while (<IN>) {
#    $_ =~ s/[\r\n]+//g;
#	$count ++;
#	print "Counting lines:\t$count\r" if $count%100000 == 0;
	
#	unless(/^#/){
#		undef @t;
#		@t = split("\t", $_);
		
#		$record{$t[0]."_".$t[1]} ++;
		
#		$length = length($t[3]);
#		if($length > 1){
#			$end = $t[1] + $length -1;
#			$record{$t[0]."_".$end} ++;
#		}
#	}
#}
#close IN;
#print "\n";
#undef $count;


##Record genome
if($genome_fa =~ /\.gz$/){
	open (IN, "gunzip -c $genome_fa |") or die $!;
}else{
	open (IN, "< $genome_fa") or die $!;
}
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	if($_ =~ /^>(\S+)/){
		$tag = $1;
		$name{$tag} = $_;
		print "Loading Genome:\t$tag\r";
	}else{
		$seq{$tag} .= $_;
	}
}
close IN;
print "\n";



## Output
if ($vcf_InDels =~ /\.vcf\.gz$/ ){			
	open (IN, "gunzip -c $vcf_InDels |") or die $!;
}elsif($vcf_InDels =~ /\.vcf$/ 	){
	open (IN, "< $vcf_InDels") or die $!;
}
open (OUT1, "> $out"."_".$date.".vcf") or die $!;
open (OUT2, "> $out"."_".$date."_table.xls") or die $!;
open (OUT3, "> $out"."_".$date."_context.fa") or die $!;
while (<IN>) {
    $_ =~ s/[\r\n]+//g;
	$count ++;
	print "Recording and outputting vcf lines: $count\r" if $count % 100000 == 0;
	
	if(/^#/){
		print OUT1 $_."\n";
		
	}else{
		undef @t;
		@t = split("\t", $_);
		
		$chr 	= $t[0];
		$site 	= $t[1];
		$ref 	= $t[3];	#VCF
		$var	= $t[4];	#VCF
		
		if($var =~ /,/ && $var !~ /\*/){
			
			undef @n;
			@n = split(",", $var);
			
			push @n, $ref;
			
			#统计多态性的个数;
			$member = scalar @n;
			
			#计算多态性的最长与最短marker;
			$max_length =1;
			$min_length =100;
			foreach $i(@n){
				$max_length = length($i) if $max_length < length($i);
				$min_length = length($i) if $min_length > length($i);
			}
			
			#计算多态性的最长与最短 的长度差;
			$diff_length = $max_length - $min_length +1;
			
			#记录每一种多态性的长度;
			undef %indel_length;
			foreach $i(@n){
				$indel_length{$i} = length($i);
			}
			undef $length_banch;
			foreach $i(sort {$indel_length{$a} <=> $indel_length{$b}} keys %indel_length){
				$length_banch .= $indel_length{$i}."nt;";
			}
			
			#记录长度差>20的次数;
			undef @x;
			@x = sort {$indel_length{$a} <=> $indel_length{$b}} keys %indel_length;
			
			$times = 0;
			for $i(0..$#x-1){
				$times ++ if (length($x[$i+1]) - length($x[$i])) >= 20;
			}
			
			#记录最长的那个InDel序列
			undef $longest;
			foreach $i( sort {$indel_length{$b} <=> $indel_length{$a}} keys %indel_length){
				$longest = $i;
				last;
			}
			
			#上下游截取100bp, 并排除存在SNP或InDel的序列;
			$left = $t[1] -100;
			$right = $t[1] + length($t[3]) -1 +100;
			#undef $judge;
			#for $i($left..$right){
			#	$judge ++ if $record{$t[0]."_".$i};
			#}
			
			#$seq_length = $right - $left +1;
			#$seq_context = substr($seq{$chr}, $left-1, $seq_length);
			
			$seq_length = 200 + $max_length;
			$seq_context = substr($seq{$chr}, $left-1, 100);
			#$seq_context .= "N" x $max_length;
			$seq_context .= $longest;
			$seq_context .= substr($seq{$chr}, ($left + 100 + length($t[3]) -1), 100);
			
			#print $chr." ".$site." ".$member." ".$max_length." ".$min_length." ".$diff_length." ".$length_banch." ".$seq_length."\n";
			
			#if($judge <= 2 && $member >= 3 && $diff_length >= 20 && $diff_length <= 100){
			if($member >= 4 && $diff_length >= 20 && $diff_length <= 100 && $times >= 2){	
				if($seq_context !~ /N/){
					
					#统计各genotype的占比数量
					undef %genotype;
					for $i(10..$#t){
						undef @m;
						@m = split(":", $t[$i]);
						$m[0] =~ s/\|/\//;
						$genotype{$m[0]} ++;
					}
					
					undef $genotype_cluster;
					$genotype_missing = 0;
					$genotype_max = 0;
					undef %genotype2;
					foreach $i (sort keys %genotype){
						#$genotype_cluster .= $i.":".$genotype{$i}.";";
						$genotype_cluster .= "\t".$i."\t".(int(($genotype{$i}/524)*10000)/100);
						
						if($i =~ /\.\/\./){
							$genotype_missing = int(($genotype{$i}/524)*10000)/100;
						}else{
							$genotype_max = int(($genotype{$i}/524)*10000)/100 if $genotype_max <= int(($genotype{$i}/524)*10000)/100;
						}
					}
					#print $genotype_missing."\t".$genotype_max."\n";
					
					if($genotype_missing <= 10 && $genotype_max <= 50){
						print OUT1 $_."\n";
						print OUT2 $t[0]."\t".$t[1]."\t".$t[2]."\t".$t[3]."\t".$t[4]."\t".$member."\t".$length_banch."\t".$genotype_cluster."\n";
						print OUT3 ">".$chr."_".$site."_context_".$left."_".$right."_".$seq_length."bp\n".$seq_context."\n";
					}
				}
			}
		}
	}
}
close IN;
close OUT1;
close OUT2;
close OUT3;
print "\n";

