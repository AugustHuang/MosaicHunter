#!/usr/bin/perl

use strict;
use Getopt::Long;

my $trim_end=0;
my $trim_intron=0;
my $trim_indel=0;
my $help=0;

$help=1 if $#ARGV!=2;

GetOptions("trim_end=i" => \$trim_end,
           "trim_intron=i" => \$trim_intron,
           "trim_indel=i" => \$trim_indel,
           "help" => \$help );

if($help)
{
	print "Usage: samtools view -h | trimBamByBlock.pl --trim_num=<int> | samtools view -Sb -\n";
	print "    --trim_end       number of bases need to be trimmed from both ends of read\n";
	print "    --trim_intron    number of bases need to be trimmed closer to splicing boundaries\n";
	print "    --trim_indel     number of bases need to be trimmed closer to alignment indels\n";
	print "    --help           display help info\n";
	exit;
}

while(my $line=<STDIN>)
{
	chomp $line;
	if ($line =~ /^@/)
	{
		print $line,"\n";
	}
	else
	{
		my @tab=split /\s+/,$line;
		my $cigar=$tab[5];
		my @seq=split //,$tab[9];
		my @qual=split //,$tab[10];
		if ($cigar =~ /[^[0-9MDNSI]/)
		{
			print STDERR $_,"\n";
			last;
		}
		my @cigar_num=split /[MNDSI]/,$cigar;
		my @cigar_patt=split /[0-9]+/,$cigar;
		
		my @intron_pos;
		my @indel_pos;
		my $pos=0;
		
		for my $i (0 .. $#cigar_num)
		{
			if ($cigar_patt[$i+1] eq "M")
			{
				$pos+=$cigar_num[$i];
			}
			elsif ($cigar_patt[$i+1] eq "D")
			{
				push (@indel_pos,$pos);
			}
			elsif ($cigar_patt[$i+1] eq "S")
			{
				$pos+=$cigar_num[$i];
			}
			elsif ($cigar_patt[$i+1] eq "I")
			{
				push (@indel_pos,$pos);
				$pos+=$cigar_num[$i];
				push (@indel_pos,$pos);
			}
			elsif ($cigar_patt[$i+1] eq "N")
			{
				push (@intron_pos,$pos);
			}
		}
		
		my $start;
		my $end;
		my $len=$#seq;  ## $len = sequence length - 1

		foreach my $pos (@intron_pos)
		{
			if ($pos-$trim_intron<0)
			{
				$start=0;
			}
			else
			{
				$start=$pos-$trim_intron;
			}
			if ($pos+$trim_intron-1>$len)
			{
				$end=$len;
			}
			else
			{
				$end=$pos+$trim_intron-1;
			}
			@seq[$start .. $end]=("N")x($end-$start+1);
			@qual[$start .. $end]=("!")x($end-$start+1);
		}
		
		foreach my $pos (@indel_pos)
		{
			if ($pos-$trim_indel<0)
			{
				$start=0;
			}
			else
			{
				$start=$pos-$trim_indel;
			}
			if ($pos+$trim_indel-1>$len)
			{
				$end=$len;
			}
			else
			{
				$end=$pos+$trim_indel-1;
			}
			@seq[$start .. $end]=("N")x($end-$start+1);
			@qual[$start .. $end]=("!")x($end-$start+1);
		}
		
		@seq[0 .. $trim_end-1]=("N")x($trim_end);
		@qual[0 .. $trim_end-1]=("!")x($trim_end);
		
		@seq[$len-$trim_end+1 .. $len]=("N")x($trim_end);
		@qual[$len-$trim_end+1 .. $len]=("!")x($trim_end);
		
		print $tab[0],"\t",$tab[1],"\t",$tab[2],"\t",$tab[3],"\t",$tab[4],"\t",$cigar,"\t",$tab[6],"\t",$tab[7],"\t",$tab[8],"\t",@seq,"\t",@qual,"\n";
	}
}