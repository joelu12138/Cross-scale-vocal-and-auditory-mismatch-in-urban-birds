############################################################
# GO enrichment analysis using clusterProfiler
# Author: Xu Liu
# Description:
#   This script performs Gene Ontology (GO) enrichment analysis
#   for a list of differentially expressed genes (gene symbols).
#   Gene symbols are converted to Entrez IDs prior to enrichment.
#   Enrichment is performed for BP, CC, and MF categories.
#
# Required R packages:
#   clusterProfiler
#   org.Hs.eg.db
#   DOSE
#   dplyr
#   ggplot2
#   tidyr
#
# Input file:
#   Index.csv
#     - Must contain a column named "Index"
#     - The "Index" column should contain gene SYMBOLs
#
# Output:
#   go_result.csv
#     - Full GO enrichment result table
############################################################


############################
# 1. Load required libraries
############################

library(clusterProfiler)   # GO enrichment analysis
library(org.Hs.eg.db)      # Human gene annotation database
library(dplyr)             # Data manipulation
library(DOSE)              # Enrichment result processing
library(ggplot2)           # Visualization (optional downstream use)
library(tidyr)             # Data reshaping


############################
# 2. Read input gene list
############################

# Read gene list (CSV file must be in working directory)
data <- read.csv("Index.csv", header = TRUE)

# Inspect data structure
head(data)
dim(data)

# Extract gene SYMBOL column
diff <- data$Index


############################
# 3. Convert SYMBOL to ENTREZID
############################

# Convert gene symbols to Entrez IDs
# Note: org.Hs.eg.db is used for Homo sapiens annotation
diff_entrez <- bitr(
  diff,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = org.Hs.eg.db
)

head(diff_entrez)


############################
# 4. Perform GO enrichment analysis
############################

# Perform GO enrichment
# ont = "all" includes:
#   BP (Biological Process)
#   CC (Cellular Component)
#   MF (Molecular Function)
#
# pAdjustMethod = "fdr" applies Benjamini-Hochberg correction
# pvalueCutoff and qvalueCutoff are set to 1 to retain all terms;
# filtering can be applied later if needed.

go_enrich <- clusterProfiler::enrichGO(
  gene = diff_entrez$ENTREZID,
  ont = 'all',
  keyType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  pAdjustMethod = "fdr",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)


############################
# 5. Convert ENTREZID back to SYMBOL
############################

# Replace ENTREZID with readable gene SYMBOLs
go_enrich <- DOSE::setReadable(
  go_enrich,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENTREZID'
)


############################
# 6. (Optional) Remove redundant GO terms
############################
# Uncomment if redundancy reduction is needed
# simplify() reduces semantic redundancy among enriched GO terms

# go_enrich <- simplify(
#   go_enrich,
#   cutoff = 0.7,
#   by = "p.adjust",
#   select_fun = min
# )


############################
# 7. Export enrichment results
############################

# Extract GO enrichment result table
go_result <- go_enrich@result

# Preview results
go_result

# Write results to CSV file
write.csv(go_result, "go_result.csv", row.names = FALSE)


############################################################
# End of script
############################################################