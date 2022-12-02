library(AnnotationHub)
library(ensembldb)
library(ggvenn)
library(ComplexHeatmap)
library(DESeq2)
library(tidyverse)

# Read Data ----

dds <- readRDS("RObjects/DESeqDataSet.interaction.rds")
res_d11 <- readRDS("RObjects/DESeqResults.interaction_d11.rds")
res_d33 <- readRDS("RObjects/DESeqResults.interaction_d33.rds")


# Annotation ----

# biomaRt
# AnnotationHub
ah <- AnnotationHub()
ah

# fetch the relevant database
mouse_ensdb <- query(ah, c("EnsDb", "Mus musculus", "102"))
mouse_ensdb <- mouse_ensdb[["AH89211"]]

# grab gene annotations and filter for existing genes
annot <- genes(mouse_ensdb, return.type = "data.frame")
annot <- annot %>% 
  dplyr::select(gene_id, gene_name, entrezid) %>% 
  dplyr::filter(gene_id %in% rownames(dds))
  
# some issues in  the annotation
# we have tidied an annotation for the practical
tidy_annot <- readRDS("RObjects/Ensembl_annotations.rds")

# join this annotation with the results table
annot_d11 <- res_d11 %>% 
  as_tibble(rownames = "GeneID") %>% 
  left_join(tidy_annot, by = "GeneID")


# Data visualisation ----

# distribution of p-values
annot_d11 %>% 
  ggplot(aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.05))


# MA plot
# average expression vs log-fold-change
plotMA(res_d11)

# shrinking the log fold change
lfc_d11 <- lfcShrink(dds, 
                     res = res_d11,
                     type = "ashr")
plotMA(lfc_d11)


d11_shrink_res <- lfc_d11 %>% 
  as_tibble(rownames = "GeneID") %>% 
  left_join(tidy_annot, by = "GeneID")

d11_shrink_res %>% 
  ggplot(aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(colour = padj < 0.05))


# EXERCISE
# Make a volcano plot for day 33 contrast
# You will have to:
# 1. Calculate the shrunken LFC values
# 2. Convert the shrunked d33 object to a tibble (retain GeneIDs)
# 3. Join it with the annotation table
# 4. Produce the volcano plot
res_d33 <- readRDS("RObjects/DESeqResults.interaction_d33.rds")
tidy_annot <- readRDS("RObjects/Ensembl_annotations.rds")

shrink_d33 <- lfcShrink(dds, 
                        res = res_d33,
                        type = "ashr")
plotMA(res_d33)
plotMA(shrink_d33)

shrink_d33 <- shrink_d33 %>% 
  as_tibble(rownames = "GeneID") %>% 
  left_join(tidy_annot, by = "GeneID")

shrink_d33 %>% 
  ggplot(aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(colour = padj < 0.05), size = 0.5) +
  labs(x = "log2(fold-change)", colour = "FDR < 5%")

shrink_d33 %>% 
  ggplot(aes(log10(baseMean), log2FoldChange)) +
  geom_point(aes(colour = padj < 0.05))


# Heatmap ----
library(ComplexHeatmap)

# differiatially expressed genes
degs <- d11_shrink_res %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>% 
  pull(GeneID)

# matrix of normalised expression
# use one of the methods in DESeq for normalising expression data
# NOT 
# TPM/FPKM/RPKM/CPM/whatever
expr_mat_d11 <- assay(vst(dds))
expr_mat_d11 <- expr_mat_d11[degs, ]

# Make heatmap
Heatmap(expr_mat_d11, 
        show_row_names = FALSE)

# Scaling expression - Zscore
# Zscore = number of standard deviations above/below the mean
t(scale(t(expr_mat_d11)))

zscore_d11 <- expr_mat_d11 %>% 
  t() %>% 
  scale() %>% 
  t()

Heatmap(zscore_d11, 
        show_row_names = FALSE)

# annotate heatmap
heat_annot1 <- HeatmapAnnotation(
  df = colData(dds)[, c("Status", "TimePoint")]
)

Heatmap(zscore_d11,
        show_row_names = FALSE,
        top_annotation = heat_annot1)
