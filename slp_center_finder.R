## Center finder
## This function finds the centers of high slope regions.
slp_center_find<-function(bw_slp_path,quantile_cutoff=0.9,method="local_center",center_width=2e4,export=TRUE){
  #

  library(rtracklayer)
  
  slp<-import(bw_slp_path)
  cut_off<-quantile(slp$score,quantile_cutoff)
  slp.high<-slp[slp$score>cut_off]
  slp.high.bed<-reduce(slp.high)
  
  if (method=="local_center"){
    slp$max<-0
    slp$max[(which(diff(sign(diff(slp$score)))==-2)+1)]<-1
    slp.sub<-slp[(slp$score>cut_off)&slp$max==1]
    
    ol<-as.data.frame(findOverlaps(slp.high.bed,
                                   slp.sub,
                                   minoverlap=1L,
                                   type="any",
                                   select ="all" ))
    
    ol$subscore<-slp.sub$score[ol$subjectHits] # enquire the slope values
    ol<-ol[order(ol$queryHits,-ol$subscore),]
    ol<-ol[!duplicated(ol$queryHits),]
    slp.center<-resize(slp.sub[ol$subjectHits],width = center_width, fix="center")
  }
  if (method=="simple_center"){
    slp.center<-resize(slp.high.bed,width = center_width, fix="center")
  }
  
  cutoff_file_ext<-paste("top",(1-quantile_cutoff)*100,sep="")
  
  if (export){
    export(object=slp.center, con=paste(tools::file_path_sans_ext(basename(bw_slp_path)),"_",
                                        cutoff_file_ext,method,".bed",sep=''))
    export(object = slp.high,con=paste(tools::file_path_sans_ext(basename(bw_slp_path)),"_",
                                       cutoff_file_ext,".bw",sep=''))
  }
  return(list(slp,slp.high,slp.high.bed))
  
  
}
