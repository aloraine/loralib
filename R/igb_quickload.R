
#' Get data frame with supported genome versions from an IGB Quickload site
#'
#' This function loads the contents.txt file from an IGB Quickload site into
#' a data frame with two columns:
#'
#' genome_version - an IGB compatible genome version, e.g., H_sapiens_Jul_2013
#' description - human-friendly, free-text description of the genome version
#'
#' @param u Url of the IGB Quickload site
#' @return Data frame read from contents.txt
#' @export
getGenomeVersionForQuickloads = function(url="http://igbquickload.org/quickload") {
  u = paste(url,"contents.txt",sep="/")
  dat = read.delim(u,header=F,as.is = T,sep="\t")
  names(dat) = c("genome_version","description")
  return(dat)
}


#' Get name(s) of file(s) IGB loads by default when a user opens a genome.
#'
#' Typically there is only one such file, and it usually contains curated,
#' canonical gene models for a genome. It is usually BED-Detail format.
#'
#' Note that the file name could be relative to the genome version
#' root or it could be an absolute URL.
#'
#' If no file is loaded by default, returns NA.
#'
#' @param u Url of the the genome directory
#' @export
getDefaultAnnotationsUrl = function(genome_version,quickload_root="http://igbquickload.org/quickload") {
  u = paste(quickload_root,genome_version,"annots.xml",sep="/")
  doc = xml2::read_xml(u)
  nodes = xml_find_all(doc,xpath="//file[@load_hint='Whole Sequence']")
  names = xml_attr("name")
  return(names)
}
