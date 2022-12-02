
# load libraries
library( tidyverse )
library(clusterProfiler)

search_kegg_organism(  )
?search_kegg_organism
search_kegg_organism( "mouse", by = "common_name")

# KEGG enrichment analysis
shrink.d11 <- readRDS("RObjects/Shrunk_Results.d11.rds")
View(shrink.d11)
dim(shrink.d11)

# extract genes of interest
sigGenes <- shrink.d11 %>% 
  drop_na( Entrez, FDR ) %>% 
  filter( FDR < 0.05 & abs(logFC) > 1 ) %>% 
  pull( Entrez )

head(sigGenes)
length(sigGenes)

# search KEGG DB
?enrichKEGG
keggRes <- enrichKEGG( gene = sigGenes, organism = "mmu")
keggRes 
keggRes %>% 
  as_tibble() %>% 
  View()

# visuvalise pathway
browseKEGG( x = keggRes, pathID = "mmu04612")

# pathview

logFC <- shrink.d11$logFC
names(logFC) <- shrink.d11$Entrez
head(logFC)
library(pathview)
pathview( gene.data = logFC,
          pathway.id = "mmu04612",
          species = "mmu",
          limit = list( gene = 20)
          )

# exercise 1
logFC <- shrink.d11 %>% 
  drop_na( FDR, Entrez) %>% 
  filter( FDR < 0.01 ) %>% 
  pull(logFC, Entrez)
head(logFC)

pathview( gene.data = logFC,
          pathway.id = "mmu04659",
          species = "mmu"
          )

# exercise 2

sigGenes <- shrink.d11 %>% 
  drop_na( GeneID, FDR ) %>% 
  filter( FDR < 0.01 & abs(logFC) > 1 ) %>% 
  pull( GeneID )
head(sigGenes)

universe <- shrink.d11$GeneID
length(universe)

library(org.Mm.eg.db)
columns(org.Mm.eg.db)
ego <- enrichGO( gene = sigGenes,
                 universe = universe,
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 pvalueCutoff = 0.01,
                 readable = TRUE,
                 keyType = "ENSEMBL"
                  )

ego %>% 
  as_tibble() %>% 
  View()

dotplot(ego)


# GSEA

rankedGenes <- shrink.d11 %>% 
  drop_na( Entrez ) %>% 
  mutate( rank = logFC) %>% 
  arrange( desc(logFC) ) %>% 
  pull( rank, Entrez)
head(rankedGenes)
tail(rankedGenes)

library(msigdbr)
term2gene <- msigdbr( species = "Mus musculus", category = "H") %>% 
  dplyr::select( gs_name, entrez_gene)
term2name <- msigdbr( species = "Mus musculus", category = "H" ) %>% 
  dplyr::select( gs_name, gs_description) %>% 
  distinct()


# run GSEA
gseaRes <- GSEA( rankedGenes,
                 TERM2GENE = term2gene,
                 TERM2NAME = term2name,
                 pvalueCutoff = 1,
                 minGSSize = 15,
                 maxGSSize = 500
                 )

gseaRes_tab <- gseaRes %>% 
  as_tibble()

top_10_gsea_logFC <- gseaRes_tab %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, wt=-p.adjust) %>% 
  dplyr::select(-core_enrichment) %>%
  mutate(across(c("enrichmentScore", "NES"), round, digits=3)) %>% 
  mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))


# gsea plot
gseaplot( gseaRes, geneSetID  = "HALLMARK_INFLAMMATORY_RESPONSE",
          title = "HALLMARK_INFLAMMATORY_RESPONSE")


# exercise 1


rankedGenes <- shrink.d11 %>% 
  drop_na( Entrez, pvalue,  logFC) %>% 
  mutate( rank = -log10( pvalue) * sign(logFC) ) %>% 
  arrange( desc(rank) ) %>% 
  pull( rank, Entrez)

tail(rankedGenes)

gseaRes_pbval <- GSEA( rankedGenes,
                       TERM2GENE = term2gene,
                       TERM2NAME = term2name,
                       pvalueCutoff = 1,
                       minGSSize = 15,
                       maxGSSize = 500
)

gseaRes_pbval_atb <- gseaRes_pbval %>% 
  as_tibble()


gseaRes_pval_tab <- gseaRes_pbval_atb %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, wt=-p.adjust) %>% 
  dplyr::select(-core_enrichment) %>%
  mutate(across(c("enrichmentScore", "NES"), round, digits=3)) %>% 
  mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))


intersect(gseaRes_pval_tab$ID, top_10_gsea_logFC$ID)
