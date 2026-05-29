library(dplyr)
library(ggplot2)
library(dplyr)

# Load common cooccurring condition MR results
setwd("")
load("<path/to/mr/common/cooccurring/conditions", verbose=TRUE)


# FDR adjustment
df_all$fdr.pval <- NA
for (outcome in unique(df_all$id.outcome)) {
 df_all[df_all$id.outcome == outcome, "fdr.pval"] <- p.adjust(
  df_all[df_all$id.outcome == outcome, "pval"], method = "fdr")
}


# Load MS MR data
ms_mr_data <- readRDS("path/to/mr/ms")

# Filter MS data to genes in common cooccurring MR dataset
ms_df_subset <- ms_df_subset[ms_df_subset$gene %in% df_all_subset$gene, ]

# Combine datasets
df_plot <- rbind(df_all_subset, ms_df_subset)

# Add CI bounds
df_plot <- df_plot %>%
 mutate(
  ci_upper = b + 1.96 * se,
  ci_lower = b - 1.96 * se
 )

gene_order <- c(
 "FCRL3",
 "ITLN1",
 "AHSG",
 "LMAN2",
 "SCGN",
 "HSPA1L",
 "FGFBP3",
 "CD59",
 "UBASH3B",
 "KLRB1",
 "INHBC",
 "WARS",
 "MAPK3",
 "CTRB2",
 "STAT3",
 "TNFSF14",
 "ICAM5"
)

# Set factor levels in this exact order
df_plot$gene <- factor(df_plot$gene, levels = gene_order)

# Order phenotypes for display
df_plot$phenotype <- factor(df_plot$phenotype, levels = c(
 "Multiple Sclerosis",
 "Major depression",
 "Generalized anxiety disorder",
 "Hypertension",
 "Hypercholesterolaemia",
 "Asthma",
 "COPD"
))


# Add a variable to control point size
df_plot$point_size <- ifelse(df_plot$phenotype == "Multiple Sclerosis", "MS", "Other")

# Get gene positions for vertical lines
gene_positions <- seq_along(levels(df_plot$gene))[-length(levels(df_plot$gene))]

# odds ratio calculations
df_plot <- df_plot %>%
 mutate(
  odds_ratio = exp(b),
  ci_lower = exp(b - 1.96 * se),
  ci_upper = exp(b + 1.96 * se)
 )

# build a dataframe for shading
shade_df <- data.frame(
 xmin = seq(0.5, length(levels(df_plot$gene)) - 0.5, by = 2),
 xmax = seq(1.5, length(levels(df_plot$gene)) + 0.5, by = 2),
 ymin = -Inf,
 ymax = Inf
)

plot <- ggplot(df_plot, aes(x = gene, y = odds_ratio, color = phenotype)) +
 # alternating grey shading
 geom_rect(data = shade_df, inherit.aes = FALSE,
           aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
           fill = "grey85", alpha = 0.5) +
 
 # baseline
 geom_hline(yintercept = 1, linetype = "solid", color = "black", size = 0.5) +
 geom_hline(yintercept = c(2,3,4), linetype = "dashed", color = "grey80") +
 
 # error bars + points
 geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2,
               position = position_dodge(width = 0.7)) +
 geom_point(aes(size = point_size), position = position_dodge(width = 0.7)) +
 
 # manual colors and sizes
 scale_color_manual(values = c(
  "Multiple Sclerosis" = "black",
  "Major depression" = "#FFD700",
  "Generalized anxiety disorder" = "#FF8C00",
  "Hypertension" = "#6A5ACD",
  "Hypercholesterolaemia" = "#8B0000",
  "Asthma" = "#FF69B4",
  "COPD" = "#C71585"
 )) +
 scale_size_manual(values = c("MS" = 2, "Other" = 1.5), guide = "none") +
 
 theme_minimal(base_size = 13) +
 theme(
  #axis.text.x = element_text(angle = 45, hjust = 1),
  panel.grid.major.x = element_blank(),
  axis.line.y = element_line(color = "black"),
  legend.position = "none"
 ) +
 labs(x = "Gene", y = "MR Odds Ratio (95% CI)")

# flip so genes are on the y-axis
plot + coord_flip()
