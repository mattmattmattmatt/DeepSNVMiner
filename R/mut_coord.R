args <- commandArgs(trailingOnly = TRUE)
input <- paste(args[1],args[2],sep="/")
data <-read.table(file=input)
out <- paste(args[1],args[3],sep="/")
jpeg(out)
genomecoords <- c()
for (i in data[,1]) {if (i/10 == as.integer(i/10)) {genomecoords[length(genomecoords)+1] = i} else {genomecoords[length(genomecoords)+1] = ''}}
barplot(data[,2],main="Location of Mutations", las = 2, names.arg=genomecoords,cex.names=0.8,cex.axis=0.8)
dev.off()