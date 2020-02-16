#' Get gene model data for a species
#'
#' For BINF 3121 Markdown "Compare gene lengths across species"
#'
#' @param url URL of BED-detail file to read
#' @return Data frame with gene model data, including log10-transformed gene lengths and species
#' @export
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

#' Get gene model annotation BED-detail file URLs for species panel.
#'
#' For BINF 3121 Markdown "Compare gene lengths across species"
#' Includes species: H_sapiens (human), M_musculus (mouse),
#' D_melanogaster (fruit fly), X_tropicalis (frog), S_cerevisiae (yeast),
#' A_thaliana (mouse-ear cress), O_sativa (rice), P_trichocarpa (cottonwood),
#' D_reri (zebrafish)
#'
#' @return Character vector with gene model BED-detail file URLs
#' @export
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

#' Get gene model data for species panel
#'
#' Used in BINF 3121 Markdown "Compare gene lengths across species"
#'
#' @return Data frame with gene model data for all species in the panel. Includes columns with species (a factor) and log10-transformed gene lengths.
#' @export
getGeneLengthsForPanel = function() {
  urls = getGeneModelUrls()
  lst = lapply(urls,getDataForOneSpecies)
  df = Reduce(rbind,lst)
  return(df)
}

#' Get gene model median lengths for species panel
#'
#' Used in BINF 3121 Markdown "Compare gene lengths across species"
#'
#' @param dat - Output of getGeneLengthsForPanel
#' @return Named numeric vector with median gene lengths, sorted by size
#' @export
getMedianGeneLengthsForPanel = function(dat) {
  species = levels(dat$species)
  medians = sapply(species,function(x){median(dat[dat$species==x,"log10length"])})
  o = order(medians,decreasing=FALSE)
  medians = medians[o]
  return(medians)
}

#' Loads chromosomes and their sizes for one genome.
#'
#' Used in BINF 3121 Markdown "Compare gene lengths across species"
#' Reads the genome.txt files from an IGB Quickload site.
#' Removes poorly-assembled or alternative chromosomes from assemblies.
#'
#' @param u URL of gene model annotation file
#' @return Data frame with chromosome names, sizes, and species (a factor)
#' @export
getGenomeStructureForOneSpecies = function(u) {
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

#' Loads chromosomes and their sizes for all genomes represented in getGeneModelUrls
#'
#' Used in BINF 3121 Markdown "Compare gene lengths across species"
#' Uses getGenomeStructureForOneSpecies.
#'
#' @return Data frame with chromosome names, sizes, and species (a factor)
#' @export
getGenomeStructuresForPanel = function() {
  urls = getGeneModelUrls()
  lst = lapply(urls,getGenomeStructureForOneSpecies)
  df = Reduce(rbind,lst)
  return(df)
}

#' Creates a named vector with species panel genome sizes
#'
#' Used in BINF 3121 Markdown "Compare gene lengths across species"
#'
#' @param dat - Output from getGenomeStructuresForAllSpecies (optional)
#' @return Named numeric vector with genome sizes
#' @export
getGenomeSizesForPanel = function(dat=NULL) {
  if (is.null(dat)) {
    dat = getGenomeStructuresForAllSpecies()
  }
  species = levels(dat$species)
  sizes = unlist(lapply(species,function(x)sum(dat[dat$species==x,]$size)))
  names(sizes)=species
  return(sizes)
}
