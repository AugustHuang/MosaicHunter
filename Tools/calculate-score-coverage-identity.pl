#!/rd/gentoo/usr/bin/perl -w
use strict;

if( @ARGV != 1 ) {
    print "Usage: $0 <psl file> \n";
    exit 1;
}

my $fh1=shift @ARGV;
my $score;
my $lc;
my $pid;
open FH1,"<$fh1";

while (<FH1>)
{
  chomp($_);
  my @rec=split(/\t/,$_);
  $score=$rec[0]+$rec[2]-$rec[1]-$rec[4]-$rec[6];
  $lc=($rec[12]-$rec[11])/$rec[10];
  if ((($rec[12] - $rec[11]) - ($rec[16] - $rec[15])) <0)
  {
    $pid=(100-(1000*($rec[1]+$rec[4]+ 3*log(1)))/($rec[0] + $rec[2] + $rec[1]) *0.1)/100;    
  }
  else
  {
    $pid=(100-(1000*($rec[1]+$rec[4]+ 3*log(1+  (($rec[12] - $rec[11]) - ($rec[16] - $rec[15]))             )))/($rec[0] + $rec[2] + $rec[1]) *0.1)/100; 
  }
  for (my $i = 0; $i<=$#rec;$i++)
  {
   print $rec[$i]."\t";
  }
  print $score."\t".$lc."\t".$pid."\n";
}
close FH1;


