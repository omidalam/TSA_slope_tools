# bw_path
# human genome assembly: "hg19" or "hg38"
higlass_file_prep<-function(bw_path,genome){
  library(rtracklayer)
  library(GenomeInfoDb)
  
  # Importing the proper negspy genome file
  gnm_file<-tempfile(pattern = "",fileext = ".txt")
  if (genome=="hg38"){
    chrom.sizes.path<-"https://raw.githubusercontent.com/omidalam/negspy/master/negspy/data/hg38/chromInfo.txt"
    download.file(url = chrom.sizes.path,destfile = gnm_file)
  }else if(genome=="hg19"){
    chrom.sizes.path<-"https://raw.githubusercontent.com/omidalam/negspy/master/negspy/data/hg19/chromInfo.txt"
    download.file(url = chrom.sizes.path, destfile = gnm_file)
  }
  gnm<-read.table(gnm_file)
  
  # Import the original bw
  bw<-import(bw_path)
  
  # Select the valid chromosomes in the original BW
  valids<-(bw[seqnames(bw)%in%as.character(gnm$V1)])
  valids<-valids[!is.na(valids$score)]
  # Clip the bedgraph based on the negspy genopme
  bedgraph<-tempfile(pattern = 'bg',fileext = ".bedGraph")
  bedgraph_clipped<-tempfile(pattern = 'bg',fileext = ".bedGraph")

  export(valids,con = bedgraph)
  clip_command<-paste("bedClip",bedgraph,gnm_file,bedgraph_clipped)
  system(clip_command)
  
  # Convert the bedGraph to bigwig
  bw_exp_path<-paste(tools::file_path_sans_ext(bw_path),"_negpy_",genome,".bw",sep='')
  bw_command<-paste("bedGraphToBigWig",bedgraph_clipped,gnm_file,bw_exp_path)
  system(bw_command)
}




