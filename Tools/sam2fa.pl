#!/usr/bin/perl
use strict;

if (@ARGV!=1){
	print "Usage: samtools view <bam file> | $0 - \n";
	print "Function: convert bam into fasta\n";
	exit 1;
}

my $file = shift @ARGV;
open (FH, "<$file");

while (<FH>){
	if (/^@/){ print $_;next;}
	chomp;
	my @tab = split /\s+/;
	my $id    = $tab[0];
        my $flag  = $tab[1];
	my $seq   = $tab[9];
	if (($flag & 16 )==16) {$seq = &opposite ($seq); }
	my $end;
	if (($flag & 128 )== 128) { $end = 2; }
	if (($flag & 64  )== 64 ) { $end = 1; }

	print ">",$id,"/",$end,"\n",$seq,"\n";
}
close FH;

sub opposite{
        chomp;
        my $final_seq;
        my $last_nucleotide=chop($_[0]);

        while ($last_nucleotide ne ""){
                (my $match_nucleotide=$last_nucleotide)=~tr/atcgATCG/tagcTAGC/;
                $final_seq.=$match_nucleotide;
                $last_nucleotide=chop($_[0]);
        }
        return $final_seq;
}

