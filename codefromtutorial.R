#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("airway")
#biocLite("AnnotationHub")

library(DESeq2)
library(airway)
dir <- system.file("extdata", package="airway", mustWork=TRUE)
gtffile <- file.path(dir, "Homo_sapiens.GRCh37.75_subset.gtf")

library("GenomicFeatures")
txdb <- makeTxDbFromGFF(gtffile, format="gtf")

ebg <- exonsBy(txdb, by="gene")
ebg[[1]]



# making a txdb from an AnnotationHub
library(AnnotationHub)
ah <- AnnotationHub()
qu <- query(ah, c("Ensembl","gtf","Saccharomyces cerevisiae","release-80"))
stopifnot(length(qu) == 1) # check there was only one match
gtf <- qu[[1]] # this actually downloads the gene model
txdbFromAHub <- makeTxDbFromGRanges(gtf)


list.files(dir)
csvfile = file.path(dir, "sample_table.csv")
(sampleTable <- read.csv(csvfile,row.names=1))
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)


seqinfo(bamfiles[1])
library("GenomicAlignments")
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

assay(se)
colSums(assay(se))
colData(se)
colData(se) = DataFrame(sampleTable)
rowRanges(se)
str(metadata(rowRanges(se)))


##########################
data("airway")
se <- airway
round(colSums(assay(se))/1e6,1)
library(DESeq2)
dds = DESeqDataSet(se, design = ~ cell + dex)

countdata = assay(se)
dim(countdata); head(countdata); tail(countdata)
coldata <- colData(se)
(ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ cell + dex))

rld = rlog(dds)
head(assay(rld))

# calculate normalization factor
dds <- estimateSizeFactors(dds)

# show raw counts
head(counts(dds))
# now as for the noralized counts
head(counts(dds, normalized=TRUE))

par(mfrow=c(1,2))
# plot of log2 normalized counts before rlog correction
plot(log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ), pch=16, cex=0.3)
# plot of rlog corrected counts
plot(assay(rld)[ , 1:2], pch=16, cex=0.3)


sampleDists <- dist( t( assay(rld) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
library("pheatmap")
library("RColorBrewer")
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


plotPCA(rld, intgroup = c("dex", "cell"))
(data <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))

library("ggplot2")
ggplot(data, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

mds <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mds, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)

