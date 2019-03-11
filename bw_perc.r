# make percentile out of bigwig files
bw_perc<-function(bigwig_path){
  bw<-import(bw_file_path)
  perc.rank <- function(x) trunc(rank(x))/length(x)
  bw$score<-perc.rank(bw$score)
  export(object = ,con = paste(tools::file_path_sans_ext(bw_file_path),"_pct.bw",sep=''))
}