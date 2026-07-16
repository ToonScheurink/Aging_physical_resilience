#Figure 4
#Assign each resilience according as four quartiles.
library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(scales)
library(purrr)
library(broom)
library(ggpubr)
library(data.table)
library(pheatmap)
library(svglite)

df <- read.csv(".../0-other-res-df.csv")

df <- df %>%
  mutate(
    physresilience_quartile = ntile(physresilience, 4)
  )

write.csv(df, ".../0-other-res-df-quartile.csv", row.names = FALSE)

df <- read.csv(".../0-other-res-df-quartile.csv")

df <- df %>%
  mutate(
    physresilience_group = ifelse(physresilience_quartile %in% c(1, 2), "Low", "High")
  )

write.csv(df, ".../0-other-res-df-group.csv", row.names = FALSE)



# Merging with library_compound_name column

# Set file path
base_path <- "..."

# Load files
metabolites <- read_csv(file.path(base_path, "0-ranked_metabolites-physres.csv"))
gnps <- read_tsv(file.path(base_path, "Library_results_with_gnps.tsv"))

# Match Metabolite with #Scan# and add Compound_Name and LibraryName as new columns
result <- metabolites %>%
  left_join(
    gnps %>% select(`#Scan#`, Compound_Name, LibraryName),
    by = c("Metabolite" = "#Scan#")
  )

# Save result
write_csv(result, file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

cat("Done! Matched rows:", sum(!is.na(result$Compound_Name)), "\n")
cat("Total rows:", nrow(result), "\n")


# Merging with classification and BA query_validation columns

# Set file path
base_path <- "..."

# Load files
metabolites <- read_csv(file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))
ba <- read_csv(file.path(base_path, "BA_annotation_results.csv"))

# Match Metabolite with #Scan# and add classification, query_validation
result <- metabolites %>%
  left_join(
    ba %>% select(`#Scan#`, classification, query_validation),
    by = c("Metabolite" = "#Scan#")
  )

# Save result (overwrite the same file)
write_csv(result, file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

cat("Done! Matched rows (classification):", sum(!is.na(result$classification)), "\n")
cat("Done! Matched rows (query_validation):", sum(!is.na(result$query_validation)), "\n")
cat("Total rows:", nrow(result), "\n")




#Figure 4a
### Wilcoxon + Bar chart: Carnitine (filter: adjusted p < 0.05; FC > 1)

# Load data
base_path <- "..."

main_df   <- read.csv(file.path(base_path, "0-other-res-df-group.csv"))
phys_meta <- read.csv(file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

# Filter carnitine features
carnitine_meta <- phys_meta %>%
  filter(grepl("carnitine", Compound_Name, ignore.case = TRUE))

cat("Number of carnitine features:", nrow(carnitine_meta), "\n")

carnitine_ids <- carnitine_meta$Metabolite

# Row-sum normalization on ALL features
all_feature_cols <- grep("^X", colnames(main_df), value = TRUE)

main_norm <- main_df %>%
  mutate(row_sum = rowSums(select(., all_of(all_feature_cols)))) %>%
  mutate(across(all_of(all_feature_cols), ~ . / row_sum)) %>%
  select(-row_sum)

# Subset carnitine features only
carnitine_cols <- paste0("X", carnitine_ids)
carnitine_cols <- intersect(carnitine_cols, colnames(main_norm))

cat("Matched carnitine columns in main_df:", length(carnitine_cols), "\n")

phys <- main_norm %>%
  filter(physresilience_group %in% c("High", "Low")) %>%
  select(physresilience_group, any_of(carnitine_cols)) %>%
  rename_with(~ sub("^X", "", .), starts_with("X"))

# Wilcoxon rank-sum test + BH correction + Log2FC
results_pval <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(Feature) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Value ~ physresilience_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# Log2FC
results_log2fc <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(physresilience_group, Feature) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = physresilience_group, values_from = mean_val) %>%
  mutate(across(c("High", "Low"), ~ ifelse(. == 0, 1e-8, .))) %>%
  mutate(Log2FC = log2(High / Low))

results <- results_log2fc %>%
  left_join(results_pval, by = "Feature")

cat("Significant features (adj_p < 0.05):", sum(results$adj_p_value < 0.05, na.rm = TRUE), "\n")

# Filter significant & add Compound_Name
label_map <- carnitine_meta %>%
  select(Metabolite, Compound_Name) %>%
  mutate(Metabolite = as.character(Metabolite))

sig_results <- results %>%
  filter(adj_p_value < 0.05 & abs(Log2FC) > 1) %>%
  mutate(Feature = as.character(Feature)) %>%
  left_join(label_map, by = c("Feature" = "Metabolite")) %>%
  group_by(Compound_Name) %>%
  mutate(label = if_else(n() > 1,
                         paste0(Compound_Name, " [", Feature, "]"),
                         Compound_Name)) %>%
  ungroup() %>%
  arrange(Log2FC)

# Bar chart
if (nrow(sig_results) == 0) {
  cat("No significant features (adj_p < 0.05) to plot.\n")
} else {
  p <- ggplot(sig_results, aes(x = reorder(label, Log2FC), y = Log2FC, fill = Log2FC > 0)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D09200")) +    coord_flip() +
    labs(
      title = "Wilcoxon Rank-Sum Test: Physical Resilience\n(Carnitine Features, BH adj_p < 0.05)",
      x     = NULL,
      y     = "Log2FC (High / Low)"
    ) +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 9),
      axis.text.y       = element_text(margin = margin(r = 5)),
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position   = "none",
      plot.margin       = margin(10, 10, 10, 120)
    )
  
  print(p)
  ggsave(file.path(base_path, "barplot_phys_carnitine_wilcoxon.svg"),
         plot = p, device = "svg", width = 20, height = max(3, nrow(sig_results) * 1.2), dpi = "retina")
  cat("Saved: barplot_phys_carnitine_wilcoxon.svg\n")
}





#Figure 4b
### Wilcoxon + Bar chart: Glutamine (filter: adjusted p < 0.05)

# Load data
base_path <- "..."

main_df   <- read.csv(file.path(base_path, "0-other-res-df-group.csv"))
phys_meta <- read.csv(file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

# Filter glutamine features
glutamine_meta <- phys_meta %>%
  filter(grepl("glutamine", Compound_Name, ignore.case = TRUE))

cat("Number of glutamine features:", nrow(glutamine_meta), "\n")

glutamine_ids <- glutamine_meta$Metabolite

# Row-sum normalization on ALL features
all_feature_cols <- grep("^X", colnames(main_df), value = TRUE)

main_norm <- main_df %>%
  mutate(row_sum = rowSums(select(., all_of(all_feature_cols)))) %>%
  mutate(across(all_of(all_feature_cols), ~ . / row_sum)) %>%
  select(-row_sum)

# Subset glutamine features only
glutamine_cols <- paste0("X", glutamine_ids)
glutamine_cols <- intersect(glutamine_cols, colnames(main_norm))

cat("Matched glutamine columns in main_df:", length(glutamine_cols), "\n")

phys <- main_norm %>%
  filter(physresilience_group %in% c("High", "Low")) %>%
  select(physresilience_group, any_of(glutamine_cols)) %>%
  rename_with(~ sub("^X", "", .), starts_with("X"))

# Wilcoxon rank-sum test + BH correction + Log2FC
results_pval <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(Feature) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Value ~ physresilience_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# Log2FC
results_log2fc <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(physresilience_group, Feature) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = physresilience_group, values_from = mean_val) %>%
  mutate(across(c("High", "Low"), ~ ifelse(. == 0, 1e-8, .))) %>%
  mutate(Log2FC = log2(High / Low))

results <- results_log2fc %>%
  left_join(results_pval, by = "Feature")

cat("Significant features (adj_p < 0.05):", sum(results$adj_p_value < 0.05, na.rm = TRUE), "\n")

# Filter significant & add Compound_Name
label_map <- glutamine_meta %>%
  select(Metabolite, Compound_Name) %>%
  mutate(Metabolite = as.character(Metabolite))

sig_results <- results %>%
  filter(adj_p_value < 0.05) %>%
  mutate(Feature = as.character(Feature)) %>%
  left_join(label_map, by = c("Feature" = "Metabolite")) %>%
  group_by(Compound_Name) %>%
  mutate(label = paste0(Compound_Name, " [", Feature, "]")) %>%
  ungroup() %>%
  arrange(Log2FC)

# Bar chart
if (nrow(sig_results) == 0) {
  cat("No significant features (adj_p < 0.05) to plot.\n")
} else {
  p <- ggplot(sig_results, aes(x = reorder(label, Log2FC), y = Log2FC, fill = Log2FC > 0)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D09200")) +
    coord_flip() +
    labs(
      title = "Wilcoxon Rank-Sum Test: Physical Resilience\n(Glutamine Features, BH adj_p < 0.05)",
      x     = NULL,
      y     = "Log2FC (High / Low)"
    ) +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 9),
      axis.text.y       = element_text(margin = margin(r = 5)),
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position   = "none",
      plot.margin       = margin(10, 10, 10, 120)
    )
  
  print(p)
  ggsave(file.path(base_path, "barplot_phys_glutamine_wilcoxon.svg"),
         plot = p, device = "svg", width = 12, height = max(3, nrow(sig_results) * 1.2), dpi = "retina")
  cat("Saved: barplot_phys_glutamine_wilcoxon.svg\n")
}





#Figure 4c
### Wilcoxon + Bar chart: Phosphocholine / PC-related (filter: adjusted p < 0.05)
# Load data
base_path <- "..."

main_df   <- read.csv(file.path(base_path, "0-other-res-df-group.csv"))
phys_meta <- read.csv(file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

# Filter phosphocholine-related features
pc_keywords <- c("phosphorylcholine", "phosphocholine", "Lyso PC", "PC\\(", "PAF", "PC-", "LPC", "SM\\(", "-PC")
pc_pattern  <- paste(pc_keywords, collapse = "|")

pc_meta <- phys_meta %>%
  filter(grepl(pc_pattern, Compound_Name, ignore.case = TRUE))

cat("Number of PC-related features:", nrow(pc_meta), "\n")

pc_ids <- pc_meta$Metabolite

# Row-sum normalization on ALL features
all_feature_cols <- grep("^X", colnames(main_df), value = TRUE)

main_norm <- main_df %>%
  mutate(row_sum = rowSums(select(., all_of(all_feature_cols)))) %>%
  mutate(across(all_of(all_feature_cols), ~ . / row_sum)) %>%
  select(-row_sum)

# Subset PC-related features only
pc_cols <- paste0("X", pc_ids)
pc_cols <- intersect(pc_cols, colnames(main_norm))

cat("Matched PC-related columns in main_df:", length(pc_cols), "\n")

phys <- main_norm %>%
  filter(physresilience_group %in% c("High", "Low")) %>%
  select(physresilience_group, any_of(pc_cols)) %>%
  rename_with(~ sub("^X", "", .), starts_with("X"))

# Wilcoxon rank-sum test + BH correction + Log2FC
results_pval <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(Feature) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Value ~ physresilience_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# Log2FC
results_log2fc <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(physresilience_group, Feature) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = physresilience_group, values_from = mean_val) %>%
  mutate(across(c("High", "Low"), ~ ifelse(. == 0, 1e-8, .))) %>%
  mutate(Log2FC = log2(High / Low))

results <- results_log2fc %>%
  left_join(results_pval, by = "Feature")

cat("Significant features (adj_p < 0.05):", sum(results$adj_p_value < 0.05, na.rm = TRUE), "\n")

# Filter significant & add Compound_Name
label_map <- pc_meta %>%
  select(Metabolite, Compound_Name) %>%
  mutate(Metabolite = as.character(Metabolite))

sig_results <- results %>%
  filter(adj_p_value < 0.05) %>%
  mutate(Feature = as.character(Feature)) %>%
  left_join(label_map, by = c("Feature" = "Metabolite")) %>%
  group_by(Compound_Name) %>%
  mutate(label = paste0(Compound_Name, " [", Feature, "]")) %>%
  ungroup() %>%
  arrange(Log2FC)


# Bar chart
if (nrow(sig_results) == 0) {
  cat("No significant features (adj_p < 0.05) to plot.\n")
} else {
  p <- ggplot(sig_results, aes(x = reorder(label, Log2FC), y = Log2FC, fill = Log2FC > 0)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D09200")) +
    coord_flip() +
    labs(
      title = "Wilcoxon Rank-Sum Test: Physical Resilience\n(PC-related Features, BH adj_p < 0.05)",
      x     = NULL,
      y     = "Log2FC (High / Low)"
    ) +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 9),
      axis.text.y       = element_text(margin = margin(r = 5)),
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position   = "none",
      plot.margin       = margin(10, 10, 10, 120)
    )
  
  print(p)
  ggsave(file.path(base_path, "barplot_phys_PC_wilcoxon.svg"),
         plot = p, device = "svg", width = 12, height = max(3, nrow(sig_results) * 1.2), dpi = "retina")
  cat("Saved: barplot_phys_PC_wilcoxon.svg\n")
}











#Figure 4d-e
#Bile acid_regression_Hurdle_Model
#Load data
base_path <- "..."
main_df   <- read.csv(file.path(base_path, "0-other-res-df-group.csv"))
phys_meta <- read.csv(file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

# Define target bile acid feature IDs
bile_acid_ids <- c("80494", "77479", "75474", "79204", "66139",
                   "65555", "71371", "79218", "71737", "72327", "78004", "62124")
bile_acid_cols <- paste0("X", bile_acid_ids)
bile_acid_cols <- intersect(bile_acid_cols, colnames(main_df))
cat("Matched bile acid columns in main_df:", length(bile_acid_cols), "\n")

# Prepare analysis dataframe with z-scored physresilience
bile_df <- main_df %>%
  select(physresilience, any_of(bile_acid_cols)) %>%
  mutate(physresilience_z = as.numeric(scale(physresilience)))

# Compound name mapping
label_map <- phys_meta %>%
  mutate(scan_number = as.character(Metabolite)) %>%
  select(scan_number, Compound_Name)

get_label <- function(col) {
  scan <- sub("^X", "", col)
  name <- label_map %>% filter(scan_number == scan) %>% pull(Compound_Name)
  if (length(name) == 0 || is.na(name[1])) return(paste0("Scan ", scan))
  return(paste0(name[1], "\n(Scan ", scan, ")"))
}

# Hurdle model function
fit_two_part_bile <- function(metab) {
  
  dat <- bile_df %>%
    dplyr::select(physresilience_z, dplyr::all_of(metab)) %>%
    dplyr::rename(abundance = dplyr::all_of(metab)) %>%
    dplyr::mutate(present = abundance > 0)
  
  out <- list()
  
  # presence/absence
  m1 <- tryCatch(
    glm(present ~ physresilience_z, data = dat, family = binomial),
    error = function(e) NULL
  )
  if (!is.null(m1)) {
    out$zero <- broom::tidy(m1) %>%
      filter(term == "physresilience_z") %>%
      mutate(component = "zero")
  }
  
  # abundance
  dat_pos <- dat %>% filter(abundance > 0)
  if (nrow(dat_pos) >= 20) {
    m2 <- tryCatch(
      lm(log(abundance) ~ physresilience_z, data = dat_pos),
      error = function(e) NULL
    )
    if (!is.null(m2)) {
      out$count <- broom::tidy(m2) %>%
        filter(term == "physresilience_z") %>%
        mutate(component = "count")
    }
  }
  
  bind_rows(out) %>% mutate(metabolite = metab)
}

# Run hurdle model across all bile acid features
bile_results <- purrr::map_dfr(bile_acid_cols, fit_two_part_bile)

# FDR correction within bile acids
bile_count <- bile_results %>%
  filter(component == "count") %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

bile_zero <- bile_results %>%
  filter(component == "zero") %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

cat("Significant in count part:", sum(bile_count$p_adj < 0.05, na.rm = TRUE), "\n")
cat("Significant in zero part:",  sum(bile_zero$p_adj  < 0.05, na.rm = TRUE), "\n")

# Select significant features
sig_bile_cols <- union(
  bile_count %>% filter(p_adj < 0.05) %>% pull(metabolite),
  bile_zero  %>% filter(p_adj < 0.05) %>% pull(metabolite)
)

cat("Total significant bile acids to plot:", length(sig_bile_cols), "\n")

# Plot function: hurdle model scatter per metabolite
plot_hurdle_bile <- function(df, metabolite_col, label = NULL, log_scale = TRUE) {
  
  # Prepare plot data
  plot_dat <- df %>%
    dplyr::select(physresilience_z, all_of(metabolite_col)) %>%
    dplyr::rename(abundance = all_of(metabolite_col)) %>%
    dplyr::mutate(
      present = abundance > 0,
      physresilience_z = as.numeric(physresilience_z)
    )
  
  # Non-zero data for count part regression
  reg_dat <- plot_dat %>% filter(abundance > 0)
  
  if (nrow(reg_dat) < 20) {
    message("Skipping ", metabolite_col, ": fewer than 20 non-zero values")
    return(NULL)
  }
  
  # Fit count part: log10-linear regression on detected samples
  m_count <- lm(log10(abundance) ~ physresilience_z, data = reg_dat)
  
  # Prediction grid for regression line and CI ribbon
  pred_dat <- data.frame(
    physresilience_z = seq(
      min(reg_dat$physresilience_z, na.rm = TRUE),
      max(reg_dat$physresilience_z, na.rm = TRUE),
      length.out = 100
    )
  )
  
  pred <- predict(m_count, newdata = pred_dat, se.fit = TRUE)
  pred_dat <- pred_dat %>%
    mutate(
      fit = 10^pred$fit,
      lwr = 10^(pred$fit - 1.96 * pred$se.fit),
      upr = 10^(pred$fit + 1.96 * pred$se.fit)
    )
  
  # Fit zero part: logistic regression for detection probability
  m_zero <- glm(present ~ physresilience_z, data = plot_dat, family = binomial)
  pred_dat$prob_present <- predict(m_zero, newdata = pred_dat, type = "response")
  
  # Pull adjusted p-values for subtitle annotation
  count_p <- bile_count %>% filter(metabolite == metabolite_col) %>% pull(p_adj)
  zero_p  <- bile_zero  %>% filter(metabolite == metabolite_col) %>% pull(p_adj)
  
  subtitle_parts <- c()
  if (length(count_p) > 0 && !is.na(count_p) && count_p < 0.05)
    subtitle_parts <- c(subtitle_parts, paste0("Count FDR = ", signif(count_p, 2)))
  if (length(zero_p) > 0 && !is.na(zero_p) && zero_p < 0.05)
    subtitle_parts <- c(subtitle_parts, paste0("Zero FDR = ", signif(zero_p, 2)))
  subtitle_text <- paste(subtitle_parts, collapse = "  |  ")
  
  # Plot title: compound name if available, otherwise scan number
  plot_title <- if (!is.null(label)) label else paste0("Scan ", sub("^X", "", metabolite_col))
  
  # Build plot
  p <- ggplot(plot_dat, aes(x = physresilience_z, y = abundance)) +
    
    # Raw data points colored by detection status
    geom_point(aes(color = present), alpha = 0.6, size = 2) +
    
    # 95% confidence ribbon from count part regression
    geom_ribbon(data = pred_dat,
                aes(x = physresilience_z, ymin = lwr, ymax = upr),
                inherit.aes = FALSE, fill = "grey70", alpha = 0.4) +
    
    # Count part regression line
    geom_line(data = pred_dat,
              aes(x = physresilience_z, y = fit),
              color = "black", linewidth = 1) +
    
    labs(
      x     = "Physical resilience (z-score)",
      y     = "Peak area",
      color = "Detected",
      title = plot_title,
      subtitle = subtitle_text
    ) +
    
    scale_color_manual(values = c(`TRUE` = "#1B9E77", `FALSE` = "#D95F02")) +
    theme_classic(base_size = 13) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))
  
  # Optional log10 y-axis
  if (log_scale) {
    p <- p + scale_y_continuous(
      trans   = "log10",
      breaks  = scales::trans_breaks("log10", function(x) 10^x),
      labels  = scales::trans_format("log10", scales::math_format(10^.x))
    )
  }
  
  return(p)
}

# Generate plots for significant bile acids only
bile_plots <- map(sig_bile_cols, function(col) {
  label <- get_label(col)
  plot_hurdle_bile(bile_df, col, label = label)
}) %>%
  compact()  # Remove NULL entries (skipped features)

names(bile_plots) <- sig_bile_cols

# Save each significant bile acid plot as individual SVG
output_dir <- file.path(base_path, "figures_bile_acids")
dir.create(output_dir, showWarnings = FALSE)

n <- length(bile_plots)

if (n == 0) {
  cat("No significant bile acids found.\n")
} else {
  iwalk(bile_plots, function(p, col) {
    filename <- paste0("bile_acid_", col, ".svg")
    filepath <- file.path(output_dir, filename)
    
    ggsave(
      filename = filepath,
      plot     = p,
      device   = "svg",
      width    = 8,
      height   = 6,
      dpi      = 300
    )
    
    cat("Saved:", filepath, "\n")
  })
  
  cat("\nDone! Total", n, "SVG files saved to:", output_dir, "\n")
}


# Figure 4f-h
# Fisher's exact test
# Load data
df <- read.csv(".../NumberOf_match_unmatch_20260412.csv",
               stringsAsFactors = FALSE)

features <- unique(df$Feature)

results_all <- data.frame()

for (feat in features) {
  cat("\n", strrep("=", 60), "\n")
  cat("Feature:", feat, "\n")
  cat(strrep("=", 60), "\n")
  
  df_feat <- df %>% filter(Feature == feat)
  healthy_row <- df_feat %>% filter(DOIDCommonName == "healthy")
  disease_rows <- df_feat %>% filter(DOIDCommonName != "healthy")
  
  if (nrow(healthy_row) == 0) {
    cat("No healthy row found, skipping\n")
    next
  }
  
  for (i in 1:nrow(disease_rows)) {
    disease <- disease_rows$DOIDCommonName[i]
    
    # Build 2x2 contingency table
    # Rows: Presence / Absence
    # Cols: Disease / Healthy
    mat <- matrix(
      c(disease_rows$Presence[i], disease_rows$Absence[i],
        healthy_row$Presence,     healthy_row$Absence),
      nrow = 2, byrow = FALSE,
      dimnames = list(
        c("Presence", "Absence"),
        c(disease, "healthy")
      )
    )
    
    cat("\n[", disease, "vs healthy ]\n")
    print(mat)
    
    # Fisher's Exact Test
    result <- fisher.test(mat)
    cat("p-value =", format(result$p.value, digits = 4),
        "| Odds Ratio =", round(result$estimate, 4), "\n")
    
    # Store result
    results_all <- rbind(results_all, data.frame(
      Feature         = feat,
      Disease         = disease,
      Disease_Presence = disease_rows$Presence[i],
      Disease_Absence  = disease_rows$Absence[i],
      Healthy_Presence = healthy_row$Presence,
      Healthy_Absence  = healthy_row$Absence,
      OddsRatio       = round(result$estimate, 4),
      P_value         = result$p.value
    ))
  }
}

# Save to CSV
output_path <- ".../Fisher_results.csv"
write.csv(results_all, output_path, row.names = FALSE)
cat("\nResults saved to:", output_path, "\n")


# Figure 4i
# Set your working directory
setwd("...")

# These packages are part of CRAN repository
#install.packages("data.table", dependencies = TRUE)
#install.packages("tidyverse", dependencies = TRUE)
#install.packages("pheatmap", dependencies = TRUE)
## Load the packages required for the analysis

# Specify the folder path - it should be the folder inside the working directory 
folder_path <- ".../FASST_ALL"

## Download/Import the ReDU metadata file - it should be in the working directory folder and NOT be in the sub-folder with the csv files from the Fast Search

# Define the filename for the ReDU metadata
processed_redu_metadata <- "all_sampleinformation.tsv"

# Check if the pre-processed metadata file exists in the working directory
if (!file.exists(file.path(getwd(), processed_redu_metadata))) {
  redu_url <- "https://redu.gnps2.org/dump"
  options(timeout = 600) # If 10 min is not enough, add more time 
  download.file(redu_url, file.path(getwd(), processed_redu_metadata), mode = "wb")
  redu_metadata <- data.table::fread(processed_redu_metadata)
} else {
  redu_metadata <- data.table::fread(processed_redu_metadata)
}

# Optional: In lieu of the previous fread() command, if memory issues occur or the program crashes, we recommend commenting out the previous command, uncommenting the line with the read_tsv() command, and following the instructions below
## Within the col_select parameter, we recommend only specifying the columns needed for a given analysis to minimize problems with memory limitations; the filename and NCBITaxonomy columns are likely to be always used
### If memory issues persist, we recommend reading in the metadata in chunks, performing the analysis desired on each chunk, and appropriately recombining the results at the end
#### redu_metadata <- readr::read_tsv("all_sampleinformation.tsv", col_select = c("filename", "NCBITaxonomy", "UBERONBodyPartName", "DOIDCommonName", "HealthStatus", "BiologicalSex", "LifeStage"), show_col_types = FALSE)

# Get the list of all .csv files in the folder
file_list <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)

# Read each .csv file and add the Compound column in each df
df_list <- lapply(file_list, function(file) {
  df <- read_csv(file)
  df$Compound <- tools::file_path_sans_ext(basename(file))
  return(df)
})

# Combine all dfs into a single df
molecules_interest <- bind_rows(df_list)

molecules_interest_filtered <- molecules_interest |> 
  dplyr::filter(`Delta Mass` >= -0.01 & `Delta Mass` <= 0.01)

#remove those with less than 2 matching peaks and cosine less than 0.8
molecules_interest_filtered2 <- molecules_interest_filtered  %>%
  dplyr::filter(`Matching Peaks` >= 2)
molecules_interest_filtered3 <- molecules_interest_filtered2  %>%
  filter(Cosine>= 0.7)

# Prepare the data tables for merging
## Create a function to extract the desired segment from the USI column

MassiveID_filename <- function(USI) {
  USI <- gsub("/", ":", USI)
  USI <- sub("\\.[^\\.]*$", "", USI)
  parts <- unlist(strsplit(USI, ":"))
  combined <- paste(parts[2], parts[length(parts)], sep = ":")
  return(combined)
}

# Apply the function to each row of the USI column in the molecules_interest
molecules_interest_filtered3$USI <- vapply(molecules_interest_filtered3$USI, MassiveID_filename, FUN.VALUE = character(1))

# Prepare the ReDU metadata USI column for merging with FASST output table
## Create a function to extract the datasetID and the last segment (filename)

ReDU_USI <- function(USI) {
  USI <- gsub("/", ":", USI)
  USI <- sub("\\.[^\\.]*$", "", USI)
  parts <- unlist(strsplit(USI, ":"))
  combined <- paste(parts[2], parts[length(parts)], sep = ":")
  return(combined)
}

# Apply the function to each row of the fxlename column in the ReDU output table
redu_metadata$USI <- vapply(redu_metadata$USI, ReDU_USI, FUN.VALUE = character(1))

# Merge the ReDU metadata table and the FASST MASST output table
ReDU_MASST <- left_join(molecules_interest_filtered3, redu_metadata, by = "USI", relationship = "many-to-many")

# Once both data tables are merged, ones can filter the table which based on the research question
## To note: not all publicly available files have associated metadata and we strongly encourage scientists to make 
### their data available with a very detailed metadata (sample information)
#### As more data are being deposited in repositories more matches will be uncovered and more results will be embedded in heatmaps 

# Standardize the body parts and Health Status
ReDU_MASST_standardize <- ReDU_MASST |> 
  dplyr::mutate(
    UBERONBodyPartName = str_replace_all(UBERONBodyPartName, 'skin of trunk|skin of pes|head or neck skin|axilla skin|skin of manus|arm skin|skin of leg', 'skin'),
    UBERONBodyPartName = str_replace_all(UBERONBodyPartName, 'blood plasma|blood serum', 'blood'),
    HealthStatus = str_replace(HealthStatus, 'Chronic Illness', 'chronic illness'),
    HealthStatus = str_replace(HealthStatus, 'Healthy', 'healthy')
  )

# Separate humans and rodents from the merged data table
df_humans <- ReDU_MASST_standardize |>  
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")

# Define a list for rodents taxonomy IDs
list_rattus_mus <- c('10088|Mus', '10090|Mus musculus', '10105|Mus minutoides', '10114|Rattus', '10116|Rattus norvegicus')

# Separate rodents
df_rodents <- ReDU_MASST_standardize |>  
  dplyr::filter(NCBITaxonomy %in% list_rattus_mus)

analyze_counts <- function(df, column_interest) {
  
  # Create a list of all unique entries in the column of interest and create a df
  df_body_parts <- df |>  distinct(across(all_of(column_interest)))
  
  # Count occurrences of each entry in the column of interest
  df_BodyPartName_counts <- df |> 
    count(across(all_of(column_interest)), name = "Counts_fastMASST")
  
  # Aggregate the number and list of unique Compounds for each entry
  compounds <- df |> 
    group_by(across(all_of(column_interest))) |> 
    summarise(Compounds = n_distinct(Compound),
              CompoundsList = toString(unique(Compound))) |> 
    ungroup()
  
  # Merge all the data into a single data frame
  combined <- df_body_parts |> 
    left_join(df_BodyPartName_counts, by = column_interest) |> 
    left_join(compounds, by = column_interest)
  
  return(combined)
}

# Get a glimpse of the number of counts per organ
body_counts_humans <- analyze_counts(df_humans, "UBERONBodyPartName")
head(body_counts_humans)

body_counts_rodents <- analyze_counts(df_rodents, "UBERONBodyPartName")
head(body_counts_rodents)

# Create a function to pivot the table for data visualization
prepare_pivot_table <- function(df, column_interest, compound) {
  
  grouped_df <- df |> 
    group_by(across(all_of(c(compound, column_interest)))) |> 
    summarise(Count = n(), .groups = 'drop')
  
  pivot_table <- grouped_df |> 
    pivot_wider(names_from = all_of(compound), values_from = Count, values_fill = list(Count = 0))
  
  return(pivot_table)
}

# Define the variables based on your research question 
## Here we are interesting in organ distribution in humans and rodents of the molecule of interest
variable <- 'UBERONBodyPartName'
pivot_table_humans <- prepare_pivot_table(df_humans, variable, 'Compound')
pivot_table_rodents <- prepare_pivot_table(df_rodents, variable, 'Compound')

# Prepare the table to be compatible with pheatmap package
humans_molecules_counts_by_bodypart <- pivot_table_humans |> 
  dplyr::arrange(UBERONBodyPartName) |> 
  tibble::column_to_rownames("UBERONBodyPartName")
# Prepare the table to be compatible with pheatmap package
rodents_molecules_counts_by_bodypart <- pivot_table_rodents |>
  dplyr::arrange(UBERONBodyPartName) |> 
  tibble::column_to_rownames("UBERONBodyPartName")

# Convert all columns to numeric for the humans df
humans_molecules_counts_by_bodypart <- humans_molecules_counts_by_bodypart |> 
  dplyr::mutate(across(everything(), as.numeric))
# Convert all columns to numeric for the rodents df
rodents_molecules_counts_by_bodypart <- rodents_molecules_counts_by_bodypart |> 
  dplyr::mutate(across(everything(), as.numeric))

# Define your chosen colors
colors_version <- c("#FFFFFF", "#C7D6F0", "#EBB0A6")
# Creating the gradient function
color_gradient <- colorRampPalette(colors_version)
# Generate 30 discrete colors from this gradient
gradient_colors <- color_gradient(30)

# The users can log scale or not the data
log_humans_molecules_counts_by_bodypart <- log10(1 + humans_molecules_counts_by_bodypart)
#write.csv(log_humans_molecules_counts_by_bodypart, 
#          file = "log_humans_molecules_counts_by_bodypart.csv", 
#          row.names = TRUE)
log_rodents_molecules_counts_by_bodypart <- log10(1 + rodents_molecules_counts_by_bodypart)
#write.csv(log_rodents_molecules_counts_by_bodypart, 
#         file = "log_rodents_molecules_counts_by_bodypart.csv", 
#         row.names = TRUE)
# Organ distribution in humans
## Use heatmap for data visualization or organ distribution - humans
### If one MS/MS spectrum is used in reverse metabolomics, the cluster_rows and cluster_cols should be set to FALSE
Organ_humans <- pheatmap(log_humans_molecules_counts_by_bodypart,
                         color = gradient_colors,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         angle_col = 90,
                         main = "Organ distribution in humans",
                         fontsize = 10,
                         cellwidth = 15,
                         cellheight = 15,
                         treeheight_row = 100,
                         fontsize_row = 12,
                         fontsize_col = 12,
                         legend_fontsize = 10,
                         border_color = NA)
Organ_humans
ggsave("Organ_distribution_in_humans_All.svg", plot = Organ_humans, width = 14, height = 10, dpi = 900)
getwd()