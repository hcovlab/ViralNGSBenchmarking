#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)

library(magrittr,quietly = TRUE, warn.conflicts=FALSE)
library(janitor,quietly = TRUE, warn.conflicts=FALSE)
library(readr,quietly = TRUE, warn.conflicts=FALSE)
library(stringr,quietly = TRUE, warn.conflicts=FALSE)
library(dplyr,quietly = TRUE, warn.conflicts=FALSE)

scenario<-data.frame(Scenario=str_sub(args[1],1,-2),Replicate=str_sub(args[2],1,-2))

basic<-row_to_names(t(read.table("bm_stats.txt", header = FALSE, sep = ":", dec = ".")),row_number=1)
report<-readr::read_tsv("./quast_results/latest/transposed_report.tsv",show_col_types = FALSE)
misassemblies<-readr::read_tsv("./quast_results/latest/contigs_reports/transposed_report_misassemblies.tsv",show_col_types = FALSE)
unaligned<-row_to_names(t(read_tsv("./quast_results/latest/contigs_reports/unaligned_report.tsv",show_col_types = FALSE)),row_number=1)
reads<-row_to_names(t(read_tsv("./quast_results/latest/reads_stats/reads_report.tsv",show_col_types = FALSE)),row_number=1)
mapping<-read.table("MappingAnalysis.txt", header = T, sep = "", dec = ".")
variants<-read.table("VariantAnalysis.txt", header = T, sep = "", dec = ".")

shiver_tm<-data.frame(row_to_names(t(read.table("shiver_tm.txt", header = FALSE, sep = ":", dec = ".")),row_number=1))
smaltalign_tm<-data.frame(row_to_names(t(read.table("smaltalign_tm.txt", header = FALSE, sep = ":", dec = ".")),row_number=1))
viralngs_tm<-data.frame(row_to_names(t(read.table("viralngs_tm.txt", header = FALSE, sep = ":", dec = ".")),row_number=1))
vpipe_tm<-data.frame(row_to_names(t(read.table("vpipe_tm.txt", header = FALSE, sep = ":", dec = ".")),row_number=1))


QUAST<-rbind(scenario,scenario,scenario,scenario)%>%
  cbind(rbind(basic,basic,basic,basic))%>%
  cbind(report)%>%
  mutate(ID=paste(Scenario,"_",Replicate,"_",Assembly,sep=""))%>%
  cbind(reads,misassemblies,unaligned,mapping[2:5,],variants)%>%
  cbind(rbind(shiver_tm,smaltalign_tm,viralngs_tm,vpipe_tm))

print(ncol(QUAST))

colnames(QUAST)<-colnames(QUAST)%>%
  str_replace_all("\\.","")%>%
  str_replace_all("# ","")%>%
  str_replace_all(" ","_")

rm(basic,report,misassemblies,unaligned,reads,mapping,variants,shiver_tm,smaltalign_tm,viralngs_tm,vpipe_tm)

write.table(QUAST,file=paste("../../../",str_sub(args[1],1,-2),"_",str_sub(args[2],1,-2),"_results.csv",sep=""),sep=",",dec=".")
