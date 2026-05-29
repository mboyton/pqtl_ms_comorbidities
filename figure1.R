library(dplyr)
library(plotly)

setwd("")

# Load data
ms_mr_data <- readRDS("path/to/mr/data")
ms_coloc <- readRDS("path/to/coloc/data")

# Join colocalization results
ms_mr_data <- ms_mr_data %>%
 mutate(exposure_gene = sub(".*_", "", exposure)) %>%
 left_join(ms_coloc, by = c("exposure_gene" = "exposure"))

# Chromosome & position info
ms_mr_data <- ms_mr_data %>%
 mutate(chr = as.numeric(chr)) %>%
 arrange(chr, snp_pos_gmean)

chr_offsets <- ms_mr_data %>%
 group_by(chr) %>%
 summarise(chr_len = max(snp_pos_gmean), .groups = "drop") %>%
 arrange(chr) %>%
 mutate(offset = lag(cumsum(chr_len), default = 0))

ms_mr_data <- ms_mr_data %>%
 left_join(chr_offsets, by = "chr") %>%
 mutate(pos_cum = snp_pos_gmean + offset)

# Colour and shape encoding
ms_mr_data <- ms_mr_data %>%
 mutate(
  bg_color = ifelse(chr %% 2 == 0, "#BFBFBF", "#737373"),
  point_color = ifelse(pval.fdr < 0.05, "red", bg_color),
  marker_symbol = ifelse(PP.H4.abf > 0.8, "diamond", "circle")
 )

# Chromosome label positions
chr_labels <- ms_mr_data %>%
 group_by(chr) %>%
 summarise(label_pos = median(pos_cum), .groups = "drop") %>%
 mutate(chr = as.character(chr))

## Add OR and OR CIs
str(ms_mr_data)
ms_mr_data$OR <- exp(ms_mr_data$b)
ms_mr_data$OR_lower <- exp(ms_mr_data$b - 1.96 * ms_mr_data$se)
ms_mr_data$OR_upper <- exp(ms_mr_data$b + 1.96 * ms_mr_data$se)


# Subset FDR >= 0.05 & < 0.05
non_sig <- ms_mr_data %>% filter(pval.fdr >= 0.05)
sig <- ms_mr_data %>% filter(pval.fdr < 0.05)

# Label positioning
custom_label_positions <- c(
 FCRL3 = "bottom left",
 ITLN1 = "right",
 AHSG = "right",
 LMAN2 = "right",
 SCGN = "top center",
 HSPA1L = "right",
 FGFBP3 = "left",
 CD59 = "bottom left",
 UBASH3B = "bottom left",
 KLRB1 = "top center",
 INHBC = "right",
 WARS = "left",
 MAPK3 = "left",
 CTRB2 = "top left",
 STAT3 = "left",
 TNFSF14 = "right",
 ICAM5 = "left"
)

sig <- sig %>%
 mutate(
  OR = exp(b),
  OR_lower = exp(b - 1.96 * se),
  OR_upper = exp(b + 1.96 * se),
  label_position = ifelse(gene %in% names(custom_label_positions),
                          custom_label_positions[gene],
                          ifelse(b > 0, "top center", "bottom center")),
  sig_color = ifelse(marker_symbol == "diamond", "#F1C40F", "#E74C3C")
 )


# Split sig data for separate error bar colors
sig_coloc <- sig %>% filter(marker_symbol == "diamond")   # coloc hits
sig_coloc <- sig_coloc %>%
 mutate(
  OR       = exp(b),
  OR_lower = exp(b - 1.96 * se),
  OR_upper = exp(b + 1.96 * se),
  err_plus  = OR_upper - OR,   # above
  err_minus = OR - OR_lower    # below
 )

sig_mr_only <- sig %>% filter(marker_symbol == "circle")  # MR-only hits
sig_mr_only <- sig_mr_only %>%
 mutate(
  OR       = exp(b),
  OR_lower = exp(b - 1.96 * se),
  OR_upper = exp(b + 1.96 * se),
  err_plus  = OR_upper - OR,   # distance above the point
  err_minus = OR - OR_lower    # distance below the point
 )


# Plotting
plot <- plot_ly() %>%
 
 # Non-significant background points
 add_trace(data = non_sig,
           x = ~pos_cum,
           y = ~OR,
           type = "scatter",
           mode = "markers",
           marker = list(size = 6),
           color = ~I(point_color),
           symbol = ~marker_symbol,
           symbols = c("circle" = "circle", "diamond" = "diamond"),
           text = ~paste0("Gene: ", gene,
                          "<br>Chr: ", chr,
                          "<br>OR: ", signif(OR, 3)),
           hoverinfo = "text",
           showlegend = FALSE) %>%
 
 
 # MR-significant only (coral red) with 95% CI
 add_trace(data = sig_mr_only,
           x = ~pos_cum,
           y = ~OR,
           type = "scatter",
           mode = "markers+text",
           marker = list(size = 8, line = list(width = 1, color = "black")),
           color = I("#E74C3C"),
           symbol = ~marker_symbol,
           symbols = c("circle" = "circle"),
           text = ~paste0("<b>", gene, "</b>"),
           textposition = ~label_position,
           textfont = list(family = "Arial, sans-serif", color = "black", size = 11),
           error_y = list(
            type = "data",
            symmetric = FALSE,
            array = ~err_plus,      # above
            arrayminus = ~err_minus,# below
            color = "#E74C3C",
            width = 2
           ),
           hoverinfo = "text",
           hovertext = ~paste0(
            "Gene: ", gene,
            "<br>Chr: ", chr,
            "<br>OR (95% CI): ", signif(OR, 3), " (",
            signif(OR_lower, 3), "–", signif(OR_upper, 3), ")",
            "<br>Beta: ", signif(b, 3),
            "<br>FDR p: ", signif(pval.fdr, 3)
           ),
           showlegend = FALSE) %>%
 
 # Colocalized + MR-significant (yellow) with 95% CI
 add_trace(data = sig_coloc,
           x = ~pos_cum,
           y = ~OR,
           type = "scatter",
           mode = "markers+text",
           marker = list(size = 8, line = list(width = 1, color = "black")),
           color = I("#F1C40F"),
           symbol = ~marker_symbol,
           symbols = c("diamond" = "diamond"),
           text = ~paste0("<b>", gene, "</b>"),
           textposition = ~label_position,
           textfont = list(family = "Arial, sans-serif", color = "black", size = 11),
           error_y = list(
            type = "data",
            symmetric = FALSE,
            array = ~err_plus,      # distance above OR
            arrayminus = ~err_minus,# distance below OR
            color = "#F1C40F",
            width = 2
           ),
           hoverinfo = "text",
           hovertext = ~paste0(
            "Gene: ", gene,
            "<br>Chr: ", chr,
            "<br>OR (95% CI): ", signif(OR, 3), " (",
            signif(OR_lower, 3), "–", signif(OR_upper, 3), ")",
            "<br>Beta: ", signif(b, 3),
            "<br>FDR p: ", signif(pval.fdr, 3),
            "<br>PP.H4: ", signif(PP.H4.abf, 3)
           ),
           showlegend = FALSE) %>%
 
 # Layout
 layout(
  title = "",
  xaxis = list(
   zeroline = FALSE,
   title = "Chromosome",
   tickmode = "array",
   tickvals = chr_labels$label_pos,
   ticktext = chr_labels$chr,
   showgrid = FALSE,
   showline = TRUE,
   linecolor = "black",
   range = c(-1e7, max(ms_mr_data$pos_cum, na.rm = TRUE) + 1e7),
   tickfont = list(family = "Arial, sans-serif", size = 12, color = "black")
  ),
  yaxis = list(
   title = "Odds Ratio (95% CI)",
   type = "log",
   tickvals = c(0.25, 0.5,
                #0.67,
                1,
                #1.5,
                2, 3, 4),
   ticktext = c("0.25", "0.5",
                #"0.67",
                "1",
                #"1.5",
                "2", "3", "4"),
   showgrid = FALSE,
   showline = TRUE,
   linecolor = "black",
   tickfont = list(family = "Arial, sans-serif", size = 12, color = "black")
  ),
  shapes = list(
   list(type = "line",
        x0 = min(ms_mr_data$pos_cum, na.rm = TRUE),
        x1 = max(ms_mr_data$pos_cum, na.rm = TRUE),
        y0 = 1, y1 = 1,
        line = list(color = "lightgrey", width = 1, dash = "solid"),
        xref = "x", yref = "y")
  )
  
 )

plot %>% config(
 toImageButtonOptions = list(
  format = "svg",
  filename = "MR_beta_CI_plot"
 )
)
