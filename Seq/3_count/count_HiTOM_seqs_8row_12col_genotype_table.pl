#! usr/bin/perl -w
#use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2020.04.23	
# Modified:			2022.09.24	Hi-TOM
# Modified:			2022.10.31	Hi-TOM 直接输出 genotype table

############################################################


my $merged_fq_file_list = shift;
my $primer_table = shift;
my $genotype_table = shift;
my $barcode_row = shift;
my $barcode_col = shift;
my $out = shift;

die  "
Usage:
perl  xxx.pl  merged_fq_file_list  primer_table  genotype_table  barcode_row  barcode_col  out

Input:
merged_fq_file_list1
primer_table_19
all_merged_out_19_2022Oct292247.xls
barcode_row
barcode_col
genotype_out

Output:
genotype_out_xxxxxx_.xls

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

my(%forward, %reverse, %forward_rc, %reverse_rc);
my(%barcode_row_group, %barcode_col_group, $rc_tag, %col_and_row);
my($x, $y, $z, $m, $n, $merged_fq);

my($indel_tag, %indel_group, $plate_tag, %plate_group);
my($col_tag, $insert_seq, $row_tag, $tag, $number);
my($threshold1, $top1, $top2);
############################################################


## Record primer_table
open (IN, "< $primer_table") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t = split("\t", $_);
	print "indel_primer:\t".$t[0]."\n";
	
	if($t[0] =~ /[Ii]n[Dd]el_?(\d+)/){
		$indel_tag = $1;
	
		$indel_group{$indel_tag} ++;
		$forward{$indel_tag} = $t[2];
		$reverse{$indel_tag} = $t[4];
	}
}
close IN;


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


## Record genotype_table
open (IN, "< $genotype_table") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t = split("\t", $_);
	
	if($t[0] !~ /^#/){
		if($t[0] =~ /[Ii]n[Dd]el_?(\d+)/){
			$indel_tag = $1;
			
			$x = scalar @t;
			$y = ($x-5)/3;
			print $indel_tag."\t".$x."\t".$y."\n";
			
			for $z(1..$y){
				print $t[5+($z-1)*3]."\n";
				${"genotype_".$indel_tag}{$t[5+($z-1)*3]} = $z;
			}
		}
	}
}
close IN;


## Record merged_fq_file_list
open (IN, "< $merged_fq_file_list") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	if($_ =~ /[Ii]n[Dd]el(\d+)_(\d+)R_merged\.fq\.gz/){
		print "Analysing merged file: $_\n";
		
		$indel_tag = $1;
		#$indel_group{$indel_tag}++;
		
		$plate_tag = $2;
		$plate_group{$plate_tag}++;
		
		print $indel_tag."\t".$plate_tag."\n";
		
		$merged_fq = $_;
		if($merged_fq =~ /\.fq\.gz$/ || $merged_fq =~ /\.fastq\.gz$/){		## ".fq.gz"	or	".fastq.gz"
			open (IN2, "gunzip -c $merged_fq |") or die $!;
		}elsif($merged_fq =~ /\.fq$/ || $merged_fq =~ /\.fastq$/){			## ".fq"  	or	".fastq"
			open (IN2, "< $merged_fq") or die $!;
		}else{
			die "Check the fastq file name!\t".$merged_fq."\n"; 
		}
		while (<IN2>){
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
				
				$forward_rc{$indel_tag} = base_reverse ($forward{$indel_tag});
				$reverse_rc{$indel_tag} = base_reverse ($reverse{$indel_tag});
				
				##8 row 与 reverse primer 绑定, 12 col 与 forward primer 绑定, insert seqs 也与 forward primer 绑定;
				if($t[$i] =~ /([ATCG]{5})$forward{$indel_tag}(\w+)$reverse_rc{$indel_tag}([ATCG]{5})/){			
					
					$col_tag = $1;
					$insert_seq = $2;
					$row_tag = base_reverse ($3);
					
					${"total_".$indel_tag}{$plate_tag."_".$row_tag."_".$col_tag} ++;
					${"count_".$indel_tag."_".$plate_tag."_".$row_tag."_".$col_tag}{$insert_seq} ++;
					#print $indel_tag."\t".$plate_tag."\t".$col_tag."\t".$row_tag."\t"."111111111111111111111\n";
						
				#}elsif($t[$i] =~ /([ATCG]{5})$reverse{$indel_tag}(\w+)$forward_rc{$indel_tag}([ATCG]{5})/){ ##此种情况不存在! 
					
					#$col_tag = base_reverse ($3);
					#$insert_seq = base_reverse ($2);
					#$row_tag = $1;
					
					#${"total_".$indel_tag}{$plate_tag."_".$row_tag."_".$col_tag} ++;
					#${"count_".$indel_tag."_".$plate_tag."_".$row_tag."_".$col_tag}{$insert_seq} ++;
					#print $indel_tag."\t".$plate_tag."\t".$col_tag."\t".$row_tag."\t"."222222222222222222222\n";
					
				}
			}
		}
		close IN2;
		print "\n";
		undef $count;
	}
}
close IN;


## Output
open (OUT, "> $out"."_".$date.".xls") or die $!;
print OUT "#Sample";
foreach $i (sort {$a <=> $b} keys %indel_group){
	print OUT "\t".$i;
}
print OUT "\n";


foreach $j (sort {$a <=> $b} keys %plate_group){
	foreach $m(sort keys %barcode_row_group){
		foreach $n(sort {$a <=> $b} keys %barcode_col_group){
		
			print OUT $j."_".$m."_".$n;
			$x = $barcode_row_group{$m};
			$y = $barcode_col_group{$n};
			
			foreach $i (sort {$a <=> $b} keys %indel_group){
				#print "total_".$i."_".$j."_".$x."_".$y."\n";
				
				if (${"total_".$i}{$j."_".$x."_".$y} && ${"total_".$i}{$j."_".$x."_".$y} >= 1000){ #每个孔的reads总数 >= 1000;
					undef $tag;
					undef $number; #循环数, 只记录Top1, Top2;
					undef $top1;
					undef $top2;
					$threshold1 = ${"total_".$i}{$j."_".$x."_".$y} *0.3;
					
					foreach $z ( sort {${"count_".$i."_".$j."_".$x."_".$y}{$b} <=> ${"count_".$i."_".$j."_".$x."_".$y}{$a}} keys %{"count_".$i."_".$j."_".$x."_".$y} ){
						$number ++;
						
						if($number == 1){
							$top1 = ${"count_".$i."_".$j."_".$x."_".$y}{$z};
							if ($top1 >= $threshold1 && ${"genotype_".$i}{$z}){
								$tag = ${"genotype_".$i}{$z};
							}
							#print $number."\t".$threshold1."\t".$top1."\t".$tag."\n";
							
						}elsif($number == 2){
							$top2 = ${"count_".$i."_".$j."_".$x."_".$y}{$z};
							if ($tag && $top1 && $top2/$top1 >= 0.3 && ${"genotype_".$i}{$z}){
								$tag = $tag."_".${"genotype_".$i}{$z};
							}
							#print $number."\t".$threshold1."\t".$top2."\t".$tag."\n";
						}
					}
					
					if($tag){
						print OUT "\t".$tag;
					}else{
						print OUT "\tNG"; ##less_genotype
					}
					
				}else{
					print OUT "\tNN"; ##less_number
				}
				
			}
			print OUT "\n";
		}
	}	
}
close OUT;
print "\n";


sub base_reverse{
	my @bre = reverse(split ('', $_[0]));
	my $baseout;
	foreach $i(@bre){
	     $baseout .= $reverse_comp{$i};
	}
	$baseout;
}
