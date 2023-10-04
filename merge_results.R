#!/usr/bin/env Rscript

filenames <- list.files(".", pattern="*.csv", full.names=TRUE)

counter=0
for (f in filenames) {
  if(counter==0){
    results<<-read.csv(file=f,header=T,sep=",",dec=".")
    counter=counter+1
  }else{
    obj<-read.csv(file=f,header=T,sep=",",dec=".")
    results<<-rbind(results, obj)
  }
}

write.table(results,file="GABM_results.csv",sep=",",dec=".")