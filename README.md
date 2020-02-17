# loralib

Functions used in teaching Statistics for Bioinformatics at UNC Charlotte.
Inspired by [rafalib](https://github.com/rafalab/rafalib).

To install from github:

```
install.packages("devtools")
library(devtools)
install_github("aloraine/loralib")
```

## Example

Draw a scatter plot showing the relationship between gene length and genome size.

```
library(loralib)
genes=getGeneLengthsForPanel()
medians=getMedianGeneLengthsForPanel(genes)
sizes=getGenomeSizesForPanel()/10**6
main="Gene length and genome size"
xlab="genome sizes (Mb)"
ylab="log10(median gene length)"
xlim=c(0,4000)
plot(medians~sizes,pch=16,xlab=xlab,ylab=ylab,las=1,col="lightblue",main=main,xlim=xlim)
text(medians~sizes,labels=names(medians),cex=0.9,font=2,pos=4)
```
