#! usr/bin/perl
use warnings;
use strict;

my $soapout = shift;
my $namelist = shift;
my $out = shift;

my (%b, %n, $tag, @t, %record);

die "

perl xxx.pl    in_file    primer3out_ePCR_verified    outfile_name

" if !defined $out;


open (IN, "< $namelist") or die $!;
while (<IN>) {
    $_ =~ s/[\r\n]+//g;
	undef @t;
	@t = split("\t", $_);
	
	if($t[0] =~ /(Chr\d+)_(\d+)_context/){
		$record{$1."_".$2} ++;
	}
}
close IN;

open (IN, "< $soapout") or die $!;
open (OUT, "> $out") or die $!;
while (<IN>){
    $_ =~ s/[\r\n]+//g;
	
	if(/^#/){
		print OUT $_."\n"; 
	}else{
		undef @t;
		@t = split("\t", $_);
		if($record{$t[0]."_".$t[1]}){
			print OUT $_."\n"; 
		}
	}
}
close IN;
close OUT;

