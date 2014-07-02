#!/usr/bin/perl
use strict;

if (@ARGV!=3){
	print "Usage: cat <input_file> | $0 <output_file1> <output_file2> <output_file3>\n";
	print "Function: split sam by its supporting allele\n";
	print "    <output_file1>    read IDs which support the ref-allele\n";
	print "    <output_file2>    read IDs which support the alt-allele\n";
	print "    <output_file3>    positions of ref-/alt-allele in each read\n";
	exit 1;
}

my ($output_file1,$output_file2,$output_file3)=@ARGV;

open OUTPUT_FILE1, ">$output_file1" or die "Can't open file:$!";
open OUTPUT_FILE2, ">$output_file2" or die "Can't open file:$!";
open OUTPUT_FILE3, ">$output_file3" or die "Can't open file:$!";

my %ref_pos=();
my %alt_pos=();

while (<STDIN>){
	chomp;
	my @tab=split /\s+/;
	my $id=$tab[0];
	my $chr=$tab[1];
    my $v_pos=$tab[3];
	my $ref_nt=$tab[4];
	my $alt_nt=$tab[5];
	my $r_start=$tab[6];
	my $cigar=$tab[7];
	my @seq=split //,$tab[8];

	my @cigar_num=split(/[MNDSI]/,$cigar);
	my @cigar_patt=split(/[0-9]+/,$cigar);
	
	my $offset=0;
	my $sum=0;
	my $r_pos=$v_pos-$r_start-1;
	
	for my $i (0 .. $#cigar_num)
	{
		if ($cigar_patt[$i+1] eq "M")
		{
			$sum+=$cigar_num[$i];
		}
		elsif ($cigar_patt[$i+1] eq "D")
		{
			$offset+=$cigar_num[$i];
		}
		elsif ($cigar_patt[$i+1] eq "S")
		{
			$sum+=$cigar_num[$i];
		}
		elsif ($cigar_patt[$i+1] eq "I")
		{
			$sum+=$cigar_num[$i];
			$offset-=$cigar_num[$i];
		}
		elsif ($cigar_patt[$i+1] eq "N")
		{
			$offset+=$cigar_num[$i];
		}
		
		if ($r_pos<=$sum)
		{
			last;
		}
	}
	
	$r_pos-=$offset;

	if ($seq[$r_pos+1] eq $ref_nt)
	{
		print OUTPUT_FILE1 $id,"\n";
		if (!exists $ref_pos{$chr."_".$v_pos})
		{
			$ref_pos{$chr."_".$v_pos}=$r_pos+2;
		}
		else
		{
			$ref_pos{$chr."_".$v_pos}.=",".($r_pos+2);
		}
	}
	elsif ($seq[$r_pos+1] eq $alt_nt)
	{
		print OUTPUT_FILE2 $id,"\n";
		if (!exists $alt_pos{$chr."_".$v_pos})
		{
			$alt_pos{$chr."_".$v_pos}=$r_pos+2;
		}
		else
		{
			$alt_pos{$chr."_".$v_pos}.=",".($r_pos+2);
		}
	}
	else
	{
		print STDERR $_,"\n";
	}
}
close FH;

foreach(keys %ref_pos)
{
	print OUTPUT_FILE3 $_,"\t","ref","\t",$ref_pos{$_},"\n";
}
foreach(keys %alt_pos)
{
	print OUTPUT_FILE3 $_,"\t","alt","\t",$alt_pos{$_},"\n";
}