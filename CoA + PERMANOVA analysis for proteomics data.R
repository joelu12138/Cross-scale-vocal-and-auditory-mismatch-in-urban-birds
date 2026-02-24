############################################################
# PCoA + PERMANOVA analysis for proteomics data
# Author: Xu Liu
#
# Description:
#   This script performs:
#     1) Bray¨CCurtis distance calculation
#     2) Principal Coordinates Analysis (PCoA)
#     3) PERMANOVA (adonis2) to test group differences
#     4) Publication-quality PCoA visualization
#
# Input files:
#   1) OTU/Protein abundance matrix (rows = proteins, columns = samples)
#   2) Metadata file including:
#        - samples (sample ID)
#        - Group (e.g., Population)
#        - Sex (optional)
#
# Output:
#   PCoA_bray.pdf
#   PCoA_bray.png
############################################################


############################
# 1. Clean environment
############################
rm(list = ls())


############################
# 2. Load required libraries
############################
library(vegan)      # distance calculation + PERMANOVA
library(ggplot2)    # plotting
library(ggpubr)     # convenient PCoA plotting


############################
# 3. Read abundance matrix
############################
# Rows = proteins, columns = samples
otu <- read.csv("abundance_matrix.csv",
                header = TRUE,
                row.names = 1)

# Transpose so rows = samples, columns = proteins
otu <- t(otu)

# Inspect
head(otu)


############################
# 4. Calculate Bray¨CCurtis distance
############################
# Commonly used for ecological / proteomic compositional data
otu.distance <- vegdist(otu, method = "bray")


############################
# 5. Perform PCoA
############################
pcoa <- cmdscale(otu.distance, eig = TRUE)

# Extract first two axes
pc12 <- as.data.frame(pcoa$points[, 1:2])
colnames(pc12) <- c("PCoA1", "PCoA2")

# Calculate percentage variance explained
variance_percent <- round(pcoa$eig / sum(pcoa$eig) * 100, 2)

# Add sample names
pc12$samples <- rownames(pc12)


############################
# 6. Read metadata
############################
metadata <- read.csv("metadata.csv", header = TRUE)

# Ensure sample names match
df <- merge(pc12, metadata, by = "samples")

# Important: ensure order consistency for PERMANOVA
metadata <- metadata[match(rownames(otu), metadata$samples), ]


############################
# 7. PERMANOVA (adonis2)
############################
# Testing whether Group explains proteomic variation
set.seed(121314)

adonis_result <- adonis2(
  otu.distance ~ Group,
  data = metadata,
  permutations = 999
)

adonis_result

# Extract R2 and p-value
adonis_label <- paste0(
  "PERMANOVA R2 = ",
  round(adonis_result$R2[1], 4),
  ", P = ",
  adonis_result$`Pr(>F)`[1]
)


############################
# 8. Plot PCoA
############################
p <- ggplot(df,
            aes(x = PCoA1,
                y = PCoA2,
                color = Group,
                shape = Sex)) +
  theme_bw() +
  geom_point(size = 5) +
  stat_ellipse(level = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = paste0("PCoA1 (", variance_percent[1], "%)"),
    y = paste0("PCoA2 (", variance_percent[2], "%)"),
    title = adonis_label
  ) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid = element_blank()
  )

p


############################
# 9. Save figures
############################
ggsave("PCoA_bray.pdf", p, width = 8, height = 8)
ggsave("PCoA_bray.png", p, width = 8, height = 8, dpi = 300)

############################################################
# End of script

############################################################
