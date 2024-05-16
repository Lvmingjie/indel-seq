#! usr/bin/perl -w
#use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2020.04.23	
# Modified:			2022.09.24	Hi-TOM
# Modified:			2022.10.30	Hi-TOM fragments

############################################################


my $input_list = shift;
my $out = shift;

die  "
Usage:
perl   xxx.pl   input_list   out_prefix

Input:
input_list1
fragments_out

Output:
fragments_out_xxxxxx_.xls

" if !defined $out;


my $date = `date +%Y%b%d%H%M`;
$date =~ s/[\r\n]+//g;


my($count, @t, $i, $j, );
my($input_file, $tag, %count_number, %count_percent);
my($fragments, $top1, $top2, $diff, $numbers, $percents, $length);
############################################################


## Record input_list
open (IN, "< $input_list") or die $!;
open (OUT, "> $out"."_".$date.".xls") or die $!;
print OUT "Plate\tRow_Col\tTotal_Number\tPercent\tFragments\tNumbers\tPercents\tLength\tTop1-Top2\n";
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	$input_file = $_;
	
	if($_ =~ /indel(\d+_\d+R)_out_[\w\d]+\.xls/ || $_ =~ /InDel(\d+_\d+R)_out_[\w\d]+\.xls/){
		$tag = $1;
        print "Analysing file: $input_file\n";
		
		open (IN2, "< $input_file") or die $!;
		while (<IN2>) {
			$_ =~ s/[\r\n]+//g;
			
			unless(/^#Row_Col/){
			
				undef @t;
				@t = split("\t", $_);
				
				if($t[1] && $t[1] >= 100){ ##针对Indel_17中, 部分孔产出的序列极少;
					
					$t[2] =~ s/%//;
					undef %count_number;
					undef %count_percent;
					for $i(1..4){
						if($t[$i*3]){
							$count_number{length($t[$i*3])} += $t[$i*3+1];
							
							$t[$i*3+2] =~ s/%//;
							$count_percent{length($t[$i*3])} += $t[$i*3+2];
						}
					}
					
					undef $count;
					undef $top1;
					undef $top2;
					undef $diff;
					undef $numbers;
					undef $percents;
					undef $length;
					
					$fragments = scalar(keys %count_number);
					
					foreach $i(sort {$count_number{$b} <=> $count_number{$a}} keys %count_number){
						$numbers .= $count_number{$i}.";";
						$percents .= $count_percent{$i}.";";
						$length .= $i."bp;";
						
						$count ++;
						if($count == 1){
							$top1 = $i;
						}elsif($count == 2){
							$top2 = $i;
						}
						
						if($top1 && $top2){
							$diff = abs($top1 - $top2);
						}
					}
					
					print OUT $tag."\t".$t[0]."\t".$t[1]."\t".$t[2]."\t".$fragments."\t".$numbers."\t".$percents."\t".$length;
					
					if($diff){
						print OUT "\t".$diff."\n";
					}else{
						print OUT "\t-\n";
					}
				}else{
					#print OUT $tag."\t".$t[0]."\t".$t[1]."\t".$t[2]."\ttoo_less\n";
					print OUT $tag."\t".$t[0]."\ttoo_less\n";
				}
			}

		}
		close IN2;
		
	}
}
close IN;
close OUT;
print "\n";

