library(DESeq2)
library(tidyverse)

txi <- readRDS("RObjects/txi.rds")
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv", col_types = "cccc")

all(colnames(txi$counts)==sampleinfo$SampleName)

simple.model <- as.formula(~ Status)
model.matrix(simple.model, data = sampleinfo)

sampleinfo <- mutate(sampleinfo, Status = fct_relevel(Status, "Uninfected"))
model.matrix(simple.model, data = sampleinfo)

# build deseq2 object

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = simple.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

# deseq2 workflow

ddsObj <- estimateSizeFactors(ddsObj.filt)

normalizationFactors(ddsObj.filt)
normalizationFactors(ddsObj)

logcounts <- log2(counts(ddsObj, normalized = FALSE) + 1)
limma::plotMA(logcounts, array = 5, ylim = c(-5,5))
abline(h=0, col = "red")

logNormalizedCounts <- log2(counts(ddsObj, normalized = TRUE) + 1)
limma::plotMA(logNormalizedCounts, array = 5, ylim = c(-5,5))
abline(h=0, col = "red")

ddsObj <- estimateDispersions(ddsObj)
plotDispEsts(ddsObj)

ddsObj <- nbinomWaldTest(ddsObj)

# deseq function

ddsObj <- DESeq(ddsObj.filt)

# results

results.simple <- results(ddsObj, alpha = 0.05)
results.simple

# Exercise 1

sum(results.simple$padj < 0.05)
sum(is.na(results.simple$padj))

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange > 0, na.rm = TRUE)

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange < 0, na.rm = TRUE)


