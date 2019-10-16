library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(rtracklayer)

#modEncode blacklist https://sites.google.com/site/anshulkundaje/projects/blacklists
blackListed<-"https://github.com/Boyle-Lab/Blacklist/raw/master/lists/ce11-blacklist.v2.bed.gz"
download.file(blackListed,destfile="./ce11-blacklist.v2.bed.gz")
system("gunzip ./ce11-blacklist.v2.bed.gz")
