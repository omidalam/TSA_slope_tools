# bw_path
# chrom.sizes.path: Negspy chrom sizes.
higlass_file_prep<-function(bw_path,chrom.sizes.path){
  library(rtracklayer)
  bw<-import(bw_path)
  
  
}
bw_path<-"~/Box_Sync/Andy_lab/TSA_BigWig/tsa_bw_hg38/k562_c1r1_20k_mw20k_hg38_slop.bw"
bw<-import(bw_path)
chroms<-c(paste("chr",1:22,sep=''),"chrX")

valids<-Rle(bw[seqnames(bw)==chroms])
seqnames(valids)
seqnames(bw)
bw
valids
1:10
