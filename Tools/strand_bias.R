#!/bin/env Rscript

args<-commandArgs(TRUE)
verbose=TRUE
input_file=args[1]
output_file=args[2]

my.sb.phred <- function(x)
{
	ref1=as.numeric(x[6])
	ref2=as.numeric(x[7])
	alt1=as.numeric(x[8])
	alt2=as.numeric(x[9])
	a=matrix(c(ref1,ref2,alt1,alt2),nrow=2)
	round(-log10(fisher.test(a,alternative="two.sided")$p.value)*10,digits=6)
}

input=read.delim(input_file,header=F)
output=cbind(input,apply(input,1,my.sb.phred))

write.table(output,output_file,quote=F,sep="\t",row.names=F,col.names=F)
