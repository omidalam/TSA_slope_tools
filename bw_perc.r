# make percentile out of bigwig files
#Arguments
# bigwig_path
# cut_off <-percentile cutoff between 0 to 1
bw_perc<-function(bigwig_path,cut_off=FALSE,bigwig_export=FALSE){
  library(rtracklayer)
  library(GenomeInfoDb)
  bw<-import(bigwig_path)
  perc.rank <- function(x) trunc(rank(x))/length(x)
  temp<-bw
  temp$score<-perc.rank(bw$score)
  cutoff_file_ext<-""
  if (cut_off){
    temp<-temp[temp$score>cut_off]
    cutoff_file_ext<-trunc(cut_off*100)
  }
  if (bigwig_export){
      export(object = temp,con = paste(tools::file_path_sans_ext(bigwig_path),"_",
                                       cutoff_file_ext,"pct.bw",sep='')
      )
    }else return(temp)
    
    
  }