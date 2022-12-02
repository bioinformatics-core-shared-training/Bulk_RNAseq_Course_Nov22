
# load libraries
library( tximport )
library(DESeq2)
library( tidyverse)


# read metadata
sampleinfo <- read_tsv("data/samplesheet.tsv")
View(sampleinfo)

# read counts data
files <- file.path( "salmon", sampleinfo$SampleName, "quant.sf" )
# add names
files <- set_names(files, sampleinfo$SampleName)

# read tx2gene
tx2gene <- read_tsv( "references/tx2gene.tsv" )


# import salmon data
txi <- tximport( files = files, tx2gene = tx2gene, type = "salmon")

class(txi)
txi$abundance %>% head()
txi$counts %>% head()
txi$length %>% head()
txi$countsFromAbundance

saveRDS( txi, file = "salmon_output/txi.rds" )

# Exercise 1
tpm <- tximport( files = files, type = "salmon", 
                 tx2gene = tx2gene, 
                 countsFromAbundance = "lengthScaledTPM")

tpm$counts

###########################################################################################
# create raw counts
rawCounts <- round( txi$counts, digits = 0)

# dims
dim(rawCounts)

class(rawCounts)
# filter low exp. genes
keep <- rowSums(rawCounts ) > 5
sum(keep)

# filter raw counts
filtCounts <- rawCounts[ keep, ]
dim(filtCounts)

# Counts distribution 
summary( filtCounts)

# bob-plot
boxplot(filtCounts, las = 2)

# raw counts mean - var relationship

plot( rowMeans(filtCounts),
      rowSds(filtCounts),
      xlim = c(1, 1000),
      ylim =  c(1, 500)
      ) 

# data transformation
# log2
logCounts <- log2( filtCounts + 1 ) 

# create color code
statusCols <- str_replace_all( sampleinfo$Status, 
                              c(Infected = "red", Uninfected = "orange"))

boxplot( logCounts, col = statusCols, las = 2)
abline( h = median( logCounts), col = "blue")

# log2 counts mean var relationship
plot( rowMeans(logCounts),
      rowSds(logCounts),
      main = "log2 trans"
      )

# VST
vst_counts <- vst(filtCounts)
boxplot( vst_counts, col = statusCols, las = 2)
abline( h = median( vst_counts), col = "blue")

plot( rowMeans(vst_counts),
      rowSds(vst_counts),
      main = "VST"
)

# Exercise 2

rlogcounts <- rlog(filtCounts)

boxplot( rlogcounts, col = statusCols, las = 2)
abline( h = median( rlogcounts), col = "blue")


# PCA

dim(filtCounts)

# get PCA data
pcDat <- prcomp( t(rlogcounts ))

# PCA plots
library(ggfortify)

autoplot(pcDat)

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5
         )

# Exercise 3

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5,
         x = 2,
         y = 3
)

# PCA labels
library(ggrepel)

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5,
         x = 1,
         y = 2
) +
  geom_text_repel( aes( x = PC1, y = PC2, label = SampleName), box.padding = 0.6)


# correct sample labels

sampleinfo <- mutate( sampleinfo, 
        Status = case_when(
          SampleName == "SRR7657882" ~ "Uninfected",
          SampleName == "SRR7657873" ~ "Infected",
          TRUE ~ Status
          
        )
        
        )



autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5,
         x = 1,
         y = 2
) 


write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")



