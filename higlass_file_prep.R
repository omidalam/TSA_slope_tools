# bw_path
# chrom.sizes.path: Negspy chrom sizes.
higlass_file_prep<-function(bw_path,chrom.sizes.path){
  library(rtracklayer)
  bw<-import(bw_path)
  
  
}
library(rtracklayer)
library(GenomeInfoDb)
bw_path<-"~/Box_Sync/Andy_lab/TSA_BigWig/tsa_bw_hg38/k562_c1r1_20k_mw20k_hg38_slop.bw"
chrom.sizes.path<-"~/git/negspy/negspy/data/hg38/chromInfo.txt"
gnm<-read.table(chrom.sizes.path)
bw<-import(bw_path)
valids<-(bw[seqnames(bw)%in%as.character(gnm$V1)])

