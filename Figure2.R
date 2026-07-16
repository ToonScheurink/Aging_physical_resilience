# Figure 2b
# A Stacked bar chart for un-annotated features using a SIRIUS result

library(dplyr)

# Define file paths
canopus_path <- ".../canopus_formula_summary.csv"

metabolites_path <- ".../0-ranked_metabolites-physres.csv"

output_path <- ".../canopus_formula_summary_merged.csv"

# Load input files
canopus_df     <- read.csv(canopus_path,     check.names = FALSE, stringsAsFactors = FALSE)
metabolites_df <- read.csv(metabolites_path, check.names = FALSE, stringsAsFactors = FALSE)

cat("canopus_formula_summary rows     :", nrow(canopus_df), "\n")
cat("0-ranked_metabolites-physres rows:", nrow(metabolites_df), "\n")

# Select only the columns needed from the metabolites file
metabolites_subset <- metabolites_df %>%
  select(Metabolite_physical, Direction_physical, AbsWeight_physical)

merged_df <- canopus_df %>%
  left_join(
    metabolites_subset,
    by = c("mappingFeatureId" = "Metabolite_physical")
  )

n_matched <- sum(!is.na(merged_df$Direction_physical))
cat("Matched rows :", n_matched, "/", nrow(canopus_df), "\n")

# Export the merged data frame as a new CSV file
write.csv(merged_df, output_path, row.names = FALSE)
cat("Output saved to:", output_path, "\n")


library(dplyr)
library(ggplot2)

# Define file paths
input_path <- ".../canopus_formula_summary_merged.csv"

output_path <- ".../stacked_bar_npc_superclass.svg"

# Load data
df <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)

cat("Total rows loaded:", nrow(df), "\n")

# Filter: keep only rows where NPC#superclass Probability > 0.7
df_filtered <- df %>%
  filter(`NPC#superclass Probability` > 0.7)

cat("Rows after probability filter (> 0.7):", nrow(df_filtered), "\n")

# Filter: keep only rows where AbsWeight_physical >= 0.02
df_filtered <- df_filtered %>%
  filter(!is.na(AbsWeight_physical), AbsWeight_physical >= 0.02)

cat("Rows after AbsWeight_physical filter (>= 0.02):", nrow(df_filtered), "\n")

# Remove rows with missing NPC#superclass or Direction_physical
df_filtered <- df_filtered %>%
  filter(
    !is.na(`NPC#superclass`),
    `NPC#superclass` != "",
    !is.na(Direction_physical),
    Direction_physical != ""
  )

cat("Rows after removing NA/empty NPC#superclass & Direction_physical:", nrow(df_filtered), "\n")

direction_levels <- sort(unique(df_filtered$Direction_physical))
df_filtered$Direction_physical <- factor(df_filtered$Direction_physical,
                                         levels = direction_levels)

# Count features per NPC#superclass × Direction_physical
df_counts <- df_filtered %>%
  group_by(`NPC#superclass`, Direction_physical) %>%
  summarise(Count = n(), .groups = "drop")

# Order superclasses by total count (descending) for readability
superclass_order <- df_counts %>%
  group_by(`NPC#superclass`) %>%
  summarise(Total = sum(Count)) %>%
  arrange(desc(Total)) %>%
  pull(`NPC#superclass`)

df_counts$`NPC#superclass` <- factor(df_counts$`NPC#superclass`,
                                     levels = superclass_order)

# Define fill colors for Direction_physical
direction_colors <- c(
  "Negative" = "#CF9102",
  "Positive" = "#2263A6"
)


# Build the stacked bar chart
p <- ggplot(df_counts,
            aes(x     = `NPC#superclass`,
                y     = Count,
                fill  = Direction_physical)) +
  
  geom_bar(stat = "identity", width = 0.7) +
  
  scale_fill_manual(
    values = direction_colors,
    name   = "Direction\n(physical)"
  ) +
  
  labs(
    title    = "NPC Superclass Distribution by Direction (Physical Resilience)",
    subtitle = "NPC#superclass Probability > 0.7; AbsWeight_physical >= 0.02; Unannotated features",
    y        = "Count"
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 11, color = "gray40"),
    axis.text.y        = element_text(size = 10),
    axis.text.x        = element_text(size = 10, angle = 45, hjust = 1),
    legend.position    = "right",
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "black"),
    axis.line.y        = element_line(color = "black"),
    panel.border       = element_blank()
  )

print(p)
message("Preview displayed.")
invisible(readline())

# Export as SVG
n_superclasses <- length(unique(df_counts$`NPC#superclass`))
svg_width      <- max(5, n_superclasses * 0.25 + 2) 

svg(output_path, width = svg_width, height = 5)
print(p)
dev.off()

cat("SVG saved to:", output_path, "\n")