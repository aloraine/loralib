% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/igb_quickload.R
\name{getGenomeVersionForQuickloads}
\alias{getGenomeVersionForQuickloads}
\title{Get data frame with supported genome versions from an IGB Quickload site}
\usage{
getGenomeVersionForQuickloads(url = "http://igbquickload.org/quickload")
}
\arguments{
\item{url}{Url of the IGB Quickload site. Optional. Default is http://igbquickload.org/quickload.}
}
\value{
Data frame read from contents.txt
}
\description{
This function loads the contents.txt file from an IGB Quickload site into
a data frame with two columns:
}
\details{
genome_version - an IGB compatible genome version, e.g., H_sapiens_Jul_2013
description - human-friendly, free-text description of the genome version
}
