#' BW_bnr 
#'
#' This function bins your BigWig files (In form of GRanges) into the desired bins. The default is 20 kb.
#' @param BW input BigWig file - GRanges object
#' @param bin bin window. Defaults to 20000 bp
#' @return binned GRanges
#' @seealso 
#' @export
#' @example 
#' x <- BW_bnr(bigwig, bin=2e4)
BW_bnr<-function(BW,bin=2e4){
  #gets a BW file and bin it over fixed distances
  require(rtracklayer)
  require(GenomicRanges)
  bins<-tileGenome(seqinfo(BW),tilewidth = bin,cut.last.tile.in.chrom = TRUE)
  BW_Rle<-coverage(BW,weight="score")
  BW_bin<-binnedAverage(bins,BW_Rle,"score")
  return(BW_bin)
}

