library(BSgenome.Celegans.UCSC.ce11)
library(rtracklayer)

#read in the public files' table, always set header so that the header is defined
publicfiles<-read.table("./kranz_files.txt", header=TRUE, stringsAsFactors = FALSE)
publicfiles
for (i in 1:dim(publicfiles)[1]) {
  download.file(publicfiles$Sample[i],destfile=paste0(publicfiles$FileName[i],"_ce10.wig"))

  #import the bigwig files
  pub<-import(paste0(publicfiles$FileName[i],"_ce10.wig"))
  pub

  #import the downloaded chain file
  chainFile<-import.chain("/home/mdas/public_data/ce10ToCe11.over.chain")

  #lift the files from ce10 to ce11
  pub_ce11<-unlist(liftOver(pub, chain = chainFile))
  #to make sure all chromosomes are present in the seqlevels and seqlengths (get from BSgenome)
  seqlevels(pub_ce11)<-seqlevels(Celegans)
  seqlengths(pub_ce11)<-seqlengths(Celegans)

  #to export the filess
  export(pub_ce11, paste0(publicfiles$FileName[i], "_ce11.bw"), "bw")
}
