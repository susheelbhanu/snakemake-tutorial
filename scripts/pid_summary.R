#!/usr/bin/env R

library("dplyr")

# read the gene clusters summary table which is the output of
# anvi-summarize command run prior:
PID <- read.table(file = snakemake@input[["PID"]], header = FALSE,
                       sep = "", quote = "", skip=1) 

# setting first column as rownames
rownames(PID) <- PID[,1]
PID[,1] <- NULL        

# getting summary of the PIDs
PID_summary <- summary(PID)

write.table(PID_summary, file = snakemake@output[["summary"]], 
                sep="\t", row.names=TRUE, quote=FALSE)
