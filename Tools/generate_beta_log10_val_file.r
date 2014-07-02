#!/bin/env Rscript

## Version: 0.3.0 (2013-10-26)
## Author: Adam Yongxin Ye @ CBI

usage <- "No enough command line parameters\n"
usage <- paste0(usage, "\n")
usage <- paste0(usage, "Usage: Rscript generate_beta_log10_val_file.R <max_depth> <output_filename>\n")
usage <- paste0(usage, "\n")
usage <- paste0(usage, " * This program is used to generate log10 beta function value table\n")
usage <- paste0(usage, "   which is read and used in Yyx_genotype_log10lik_with_precalc_beta\n")
usage <- paste0(usage, "\n")
usage <- paste0(usage, "Version: 0.3.1 (2014-01-10)\n")
usage <- paste0(usage, "Author: Adam Yongxin Ye @ CBI\n")
usage <- paste0(usage, "\n")

args <- commandArgs(TRUE)
if(length(args)<2){
	stop(usage);
}
max_depth <- args[1]
output_filename <- args[2]

depth <- 1
#tmp_vec = log10(beta(1+(0:depth),1+(depth:0)))
tmp_vec = lbeta(1+(0:depth),1+(depth:0)) / log(10)
write.table(t(tmp_vec),output_filename,append=F,quote=F,sep=",",row.names=F,col.names=F)
for(depth in 2:max_depth){
	tmp_vec = lbeta(1+(0:depth),1+(depth:0)) / log(10)
	write.table(t(tmp_vec),output_filename,append=T,quote=F,sep=",",row.names=F,col.names=F)
}
