#!/rd/gentoo/usr/bin/perl -w
use strict;

if (@ARGV != 4){
	print "This program is to pick up the alignment of the highest score from BLAT result; All these result has passed percent identity and length coverage criterion. The 22th field is the score.\n";
	print "    The first parameter is the input file of calculated BLAT result.\n";
	print "    The second parameter is the input file of calculated BLAT result.\n";
	print "    The third parameter is the output file with single highest-score hit.\n";
	print "    The fourth parameter is the output file with more than one highest-score hits.\n";
	exit 1;
}

my $file1=shift @ARGV;
my $file2=shift @ARGV;
my $file3=shift @ARGV;
my $file4=shift @ARGV;
my %highest;

open (FH1,"<$file1");
while (<FH1>){
	chomp;
	my @row = split /\t/;
	my $score= $row[21];
	my $qName = $row[9];
	my $block =$row[17];

	if (!exists $highest{$qName}){
		$highest{$qName}[0]=$score;
		$highest{$qName}[1]=$block;
		$highest{$qName}[2]=0; #record for highest hit number, if only one highest result, here is 0; else here is 1;
	}
	else{
		if($score > $highest{$qName}[0]){
			$highest{$qName}[2]=0;
			$highest{$qName}[0]=$score;	
			$highest{$qName}[1]=$block;
		}
		elsif($score == $highest{$qName}[0]){
			if($block >$highest{$qName}[1]){
				$highest{$qName}[1]=$block;
				$highest{$qName}[2]=0;
			}
			elsif($block == $highest{$qName}[1]){
				$highest{$qName}[2]=1;
			}
		}
	}
}
close FH1;

open (SI,">$file3");
open (MU,">$file4");
open (FH2,"<$file2");
while (<FH2>){
	chomp;
	my @row = split /\t/;
    my $qName = $row[9];
	my $score=$row[21];
	my $block=$row[17];

	if($score == $highest{$qName}[0] && $block == $highest{$qName}[1]){
		if($highest{$qName}[2]==0){
			print SI $_,"\n";
		}
		if($highest{$qName}[2]==1){
			print MU $_,"\n";
		}
	}
}
close FH2;
close SI;
close MU;
