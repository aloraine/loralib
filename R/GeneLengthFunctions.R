# download the file from the given url if not already present in the current working directory
# read the file and return a data frame with named columns
getDataForOneSpecies = function(url) {
  fname = tail(unlist(strsplit(url,"/")),1)
  species = ifelse(length(grep("Araport11",fname))>0,"A_thaliana",paste(unlist(strsplit(fname,"_"))[1:2],collapse = "_"))
  if (! file.exists(fname) ) {
    download.file(url,fname)
  }
  data = read.delim(fname,as.is=T,sep="\t",header=FALSE,quote="")
  genes = data[,c(1,2,3,4,6,10,13,14)]
  column_names=c("chr","start","end","id","strand","num_exons",
                 "locus","description")
  names(genes)=column_names
  genes$species = factor(species)
  genes$log10length=log10(genes$end-genes$start)
  genes=genes[,c("species","log10length",column_names)]
  return(genes)
}

# return a vector containing URLs of gene model (bed-detail) files 
getGeneModelUrls = function() {
  urls = c("http://igbquickload.org/quickload/H_sapiens_Dec_2013/H_sapiens_Dec_2013_ncbiRefSeqCurated.bed.gz",
            "http://igbquickload.org/quickload/M_musculus_Dec_2011/M_musculus_Dec_2011_refGene.bed.gz",
            "http://igbquickload.org/quickload/D_melanogaster_Jul_2014/D_melanogaster_Jul_2014.bed.gz",
            "http://quickload.bioviz.org/quickload/X_tropicalis_Jul_2016/X_tropicalis_Jul_2016_ncbiRefSeq.bed.gz",
            "http://quickload.bioviz.org/quickload/S_cerevisiae_Apr_2011/S_cerevisiae_Apr_2011_ncbiRefSeq.bed.gz",
            "http://igbquickload.org/quickload/A_thaliana_Jun_2009/Araport11.bed.gz",
            "http://igbquickload.org/quickload/O_sativa_japonica_Oct_2011/O_sativa_japonica_Oct_2011.bed.gz",
            "http://igbquickload.org/quickload/P_trichocarpa_Aug_2012/P_trichocarpa_Aug_2012.bed.gz",
            "http://quickload.bioviz.org/quickload/D_rerio_May_2017/D_rerio_May_2017_refGene.bed.gz"
            )
  return(urls)
}

# get a giant data frame with gene models for all the species
getDataForAllSpecies = function() {
  urls = getGeneModelUrls()
  lst = lapply(urls,getDataForOneSpecies)
  df = Reduce(rbind,lst)
  return(df)
}

# accepts a data frame with gene model lengths, output from getDataForAllSpecies
# get a sorted, named vector of median lengths by species
getMediansBySpecies = function(dat) {
  species = levels(dat$species)
  medians = sapply(species,function(x){median(dat[dat$species==x,"log10length"])})
  o = order(medians,decreasing=FALSE)
  medians = medians[o]
  return(medians)
}

getGenomeStructure = function(u) {
  genome_size_url = paste(c(head(unlist(strsplit(u,"/")),-1),"genome.txt"),collapse="/")
  fname = tail(unlist(strsplit(u,"/")),1)
  species = ifelse(length(grep("Araport11",fname))>0,"A_thaliana",paste(unlist(strsplit(fname,"_"))[1:2],collapse = "_"))
  dat = read.delim(genome_size_url,header=F,as.is = T)
  names(dat)=c("chr","size")
  # UCSC includes many alternative or poorly assembled sequences with names
  # like chr1_alt
  dat = dat[grep("_",dat$chr,invert=T),]
  # Flybase has many chromosomes with names starting 211
  dat = dat[grep("^211",dat$chr,invert=T),]
  dat$species=species
  dat=dat[,c("species","chr","size")]
  dat$species=factor(dat$species)
  return(dat)
}

getAllGenomeStructures = function() {
  urls = getGeneModelUrls()
  lst = lapply(urls,getGenomeStructure)
  df = Reduce(rbind,lst)
  return(df)
}

# accepts a data frame of chromosome sizes, the output from getAllGenomeStructures
# returns a named vector with genome size
getGenomeSizes = function(d) {
  species = levels(d$species)
  sizes = unlist(lapply(species,function(x)sum(d[d$species==x,]$size)))
  names(sizes)=species
  return(sizes)
}