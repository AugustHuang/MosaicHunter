#!/bin/env Rscript

args<-commandArgs(TRUE)
verbose=TRUE
input_file=args[1]
output_file=args[2]
phred_cutoff=as.numeric(args[3])

my.wilcox.phred <- function(x)
{
	ref_pos=as.numeric(unlist(strsplit(x[3],",")))
	alt_pos=as.numeric(unlist(strsplit(x[4],",")))
	round(-log10(wilcox.test(ref_pos,alt_pos)$p.value)*10,digits=6)
}

input=read.delim(input_file,header=F,stringsAsFactors=F)
phred=apply(input,1,my.wilcox.phred)
filtered=input[phred>phred_cutoff,]
output=data.frame(filtered[,1],filtered[,2]-1,filtered[,2])

write.table(output,output_file,quote=F,sep="\t",row.names=F,col.names=F)