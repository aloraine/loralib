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
Then draw a boxplot showing distribution of gene lengths by species.

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
old.par=par(no.readonly=TRUE)
par(mar=c(5.1,7.5,4.1,2.1))
boxplot(log10length~species,data=genes,las=1,horizontal=TRUE,xlab=ylab)
par(old.par)
```


