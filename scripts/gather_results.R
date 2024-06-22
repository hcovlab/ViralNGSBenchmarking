#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)

library(magrittr,quietly = TRUE, warn.conflicts=FALSE)
library(janitor,quietly = TRUE, warn.conflicts=FALSE)
library(readr,quietly = TRUE, warn.conflicts=FALSE)
library(stringr,quietly = TRUE, warn.conflicts=FALSE)
library(dplyr,quietly = TRUE, warn.conflicts=FALSE)

if(str_sub(args[1],1,3)=="SIM"){
  scenario<-data.frame(Scenario=str_sub(args[1],1,-2),Replicate=str_sub(args[2],1,-2))
}else{
  scenario<-data.frame(Dataset=str_sub(args[1],1,-2),Sample=str_sub(args[2],1,-2)) # sgs and ss_ngs
}

if(str_sub(args[1],1,3)!="NGS"){
  basic<-row_to_names(t(read.table("bm_stats.txt", header = FALSE, sep = ":", dec = ".")),row_number=1)
  print(paste("basic:",ncol(basic)))
  report<-readr::read_tsv("./quast_results/latest/transposed_report.tsv",show_col_types = FALSE)
  print(paste("report:",ncol(report)))
  misassemblies<-readr::read_tsv("./quast_results/latest/contigs_reports/transposed_report_misassemblies.tsv",show_col_types = FALSE)
  print(paste("misassemblies:",ncol(misassemblies)))
  unaligned<-row_to_names(t(read_tsv("./quast_results/latest/contigs_reports/unaligned_report.tsv",show_col_types = FALSE)),row_number=1)
  print(paste("unaligned:",ncol(unaligned)))
  reads<-row_to_names(t(read_tsv("./quast_results/latest/reads_stats/reads_report.tsv",show_col_types = FALSE)),row_number=1)
  print(paste("reads:",ncol(reads)))
  mapping<-read.table("MappingAnalysis.txt", header = T, sep = "", dec = ".")
  print(paste("mapping:",ncol(mapping)))
  variants<-read.table("VariantAnalysis.txt", header = T, sep = "", dec = ".")
  print(paste("variants:",ncol(variants)))

  # additonal lines for old and new SIM compatibility
  droplist=c("N90","NG90","auN","auNG","L90","LG90","# structural variations","NA90","NGA90","auNA","auNGA","LA90","LGA90","N75","NG75","L75","LG75","NA75","NGA75","LA75","LGA75")
  excludelist=intersect(droplist,colnames(report))
  report<-select(report,-excludelist)
  #-------
  droplist2=c("# structural variations")
  excludelist2=intersect(droplist2,colnames(misassemblies))
  misassemblies<-select(misassemblies,-excludelist2)

}else{
  size<-read.table("readsize.txt", header = T, sep = "", dec = ".")
  print(paste("dataset size:",ncol(size)))
}


shiver_tm<-data.frame(row_to_names(t(read.table("shiver_tm.txt", header = FALSE, sep = ":", dec = ".")),row_number=1))
smaltalign_tm<-data.frame(row_to_names(t(read.table("smaltalign_tm.txt", header = FALSE, sep = ":", dec = ".")),row_number=1))
viralngs_tm<-data.frame(row_to_names(t(read.table("viralngs_tm.txt", header = FALSE, sep = ":", dec = ".")),row_number=1))
vpipe_tm<-data.frame(row_to_names(t(read.table("vpipe_tm.txt", header = FALSE, sep = ":", dec = ".")),row_number=1))


if(str_sub(args[1],1,3)=="SIM"){
  QUAST<-rbind(scenario,scenario,scenario,scenario)%>%
  cbind(rbind(basic,basic,basic,basic))%>%
  cbind(report)%>%
  mutate(ID=paste(Scenario,"_",Replicate,"_",Assembly,sep=""))%>% # sim
  cbind(reads,misassemblies,unaligned,mapping[2:5,],variants)%>%
  cbind(rbind(shiver_tm,smaltalign_tm,viralngs_tm,vpipe_tm))
  rm(report,misassemblies,unaligned,reads,mapping,variants,shiver_tm,smaltalign_tm,viralngs_tm,vpipe_tm)
}else if(str_sub(args[1],1,3)!="NGS"){
  QUAST<-rbind(scenario,scenario,scenario,scenario)%>%
  cbind(rbind(basic,basic,basic,basic))%>%
  cbind(report)%>%
  mutate(ID=paste(Dataset,"_",Sample,"_",Assembly,sep=""))%>% # sgs and ss_ngs
  cbind(reads,misassemblies,unaligned,mapping[2:5,],variants)%>%
  cbind(rbind(shiver_tm,smaltalign_tm,viralngs_tm,vpipe_tm))
  rm(report,misassemblies,unaligned,reads,mapping,variants,shiver_tm,smaltalign_tm,viralngs_tm,vpipe_tm)
}else{
  QUAST<-rbind(scenario,scenario,scenario,scenario)%>%
  cbind(cbind(size))%>%
  cbind(rbind(shiver_tm,smaltalign_tm,viralngs_tm,vpipe_tm))
  rm(shiver_tm,smaltalign_tm,viralngs_tm,vpipe_tm,size)
}

print(ncol(QUAST))

colnames(QUAST)<-colnames(QUAST)%>%
  str_replace_all("\\.","")%>%
  str_replace_all("# ","")%>%
  str_replace_all(" ","_")

write.table(QUAST,file=paste("../../../",str_sub(args[1],1,-2),"_",str_sub(args[2],1,-2),"_results.csv",sep=""),sep=",",dec=".")
