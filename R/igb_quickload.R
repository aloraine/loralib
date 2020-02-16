
#' Get data frame with genome version directories from an IGB Quickload site
#'
#' This function loads the contents.txt file from an IGB Quickload site.
#'
#' @param u Url of the IGB Quickload site
#' @return Data frame with contents.txt
#' @export
getGenomeVersionUrls = function(url="http://igbquickload.org/quickload") {
  u = paste(url,"contents.txt",sep="/")
  dat = read.delim(u,header=F,as.is = T,sep="\t")
  names(dat) = c("genome_dir","description")
  return(dat)
}
