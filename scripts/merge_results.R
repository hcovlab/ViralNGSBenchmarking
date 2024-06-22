#!/usr/bin/env Rscript

filenames <- list.files(".", pattern="*.csv", full.names=TRUE)

counter=0
for (f in filenames) {
  if(counter==0){
    results<<-read.csv(file=f,header=T,sep=",",dec=".")
    #deletecols<-c("N90","N75","NG90","NG75","auN","auNG","L90","L75","LG90","LG75","structural_variations","structural_variations.1","NA90","NA75","NGA90","NGA75","auNA","auNGA","LA90","LA75","LGA90","LGA75")
    #deletecols<-c()
    #print(colnames(results))
    #print(colnames(results)[names(results) %in% deletecols])
    #results<<-results[ , -which(names(results) %in% deletecols) ]
    counter=counter+1
    print(ncol(results))
  }else{
    obj<-read.csv(file=f,header=T,sep=",",dec=".")
    #obj2<-obj[ , -which(names(obj) %in% deletecols) ]
    #print(ncol(obj2))
    #print(colnames(obj)[names(obj) %in% deletecols])
    results<<-rbind(results, obj)
    print(paste(f,nrow(results)))
  }
}
print(results)
write.table(results,file="GABM_results.csv",sep=",",dec=".") # all