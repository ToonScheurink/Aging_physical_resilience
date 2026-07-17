### Packages
library(data.table)
library(dplyr)
library(tibble)
library(ggpubr)
library(stringr)
library(tidyr)
library(caret)
library(mixOmics)
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(ppcor)
library(lme4)
library(lmerTest)       
library(broom.mixed)   
library(future.apply) 
library(patchwork)
library(broom)
library(forcats)
library(flextable)
library(gtsummary)
library(gt)      
library(zCompositions)
library(compositions) 
library(pls)
library(pheatmap)
library(MASS)
library(viridis)
library(ggvegan)
library(ggfortify)
library(miaViz)
library(ggrepel)
library(gridExtra)
library(kableExtra)

df <- readRDS("clean_data_wellcome.rds")

df$host_age<-as.numeric(df$host_age)
df$host_body_mass_index<-as.numeric(df$host_body_mass_index)
df$physresilience_pctl<-as.numeric(df$physresilience_pctl)
df$physresilience<-as.numeric(df$physresilience)



df <- df %>%
  group_by(record_ID) %>%
  slice_min(order_by = host_age, with_ties = FALSE) %>%
  ungroup()

#antihypertensives
cols_to_check <- c("10712", "38586", "55898", "70755")

df <- df %>%
  filter(if_all(all_of(cols_to_check), ~ . == 0))#Selction of metabolites only
df_metab <- df %>%
  column_to_rownames(var = "filename") %>%  
  dplyr::select(3:22663) %>%                
  filter(rowSums(.) > 0)           

#Counting metabolites
n_present_metabolites <- sum(colSums(df_metab != 0) > 0)
n_present_metabolites

#RCLR conversion
df_rclr <- decostand(df_metab, method = "rclr")
df_rclr_matrix <- as.matrix(df_rclr)


#Removing zero variance metabolites
var_resid <- apply(df_rclr_matrix, 2, var)
metab_resid_filtered <- df_rclr_matrix[, var_resid > 1e-6]


X <- metab_resid_filtered
Y <- df$physresilience

var_resid <- apply(metab_resid_filtered, 2, var)
top_idx <- order(var_resid, decreasing = TRUE)[1:5000]
X_sub <- metab_resid_filtered[, top_idx]

spls_final <- spls(X_sub, Y, ncomp = 1, keepX = 2810) 

# Get component scores
comp1_scores <- spls_final$variates$X[,1]

# Fit linear model: CR ~ sPLS Component 1
model <- lm(CR ~ Comp1, data = data.frame(Comp1 = comp1_scores, CR = Y))
summary_model <- summary(model)

summary_model
# Extract statistics
beta <- round(summary_model$coefficients[2, 1], 2)   # slope
p_val <- signif(summary_model$coefficients[2, 4], 3) # p-value
r2 <- round(summary_model$r.squared, 2)             # R-squared

# Set position for annotation
x_pos <- max(comp1_scores) - 0.05 * diff(range(comp1_scores))
y_pos <- min(Y) + 0.25 * diff(range(Y))

# Create plot with annotation
sPLS_physres <- ggplot(
  data.frame(Comp1 = comp1_scores, CR = Y),
  aes(x = Comp1, y = CR)
) +
  geom_point(alpha = 0.8, color = "#2C3E50", size = 3) +
  geom_smooth(method = "lm", color = "#0072B2", fill = "#0072B2",
              alpha = 0.2, se = TRUE, size = 1.2) +
  labs(
    x = "sPLS Component 1 Score",
    y = NULL
  ) +
  annotate(
    "label",
    x = x_pos, y = y_pos,
    label = paste0("R² = ", round(r2, 2),
                   "\nβ = ", round(beta, 2),
                   "\np = ", format.pval(p_val, 2)),
    hjust = 1, vjust = 1,
    size = 8,
    label.size = 0.3,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray40"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),       
    axis.ticks.length = unit(3, "pt"),               
    plot.margin = margin(15, 15, 15, 15)
  )


sPLS_physres

#ggsave("sPLS_physres_no_antihypertensives.svg", plot = sPLS_physres, dpi =300, width = 10, height = 7)


set.seed(500)

# Define cross-validation setup 
folds <- createFolds(Y, k = 10, list = TRUE)

cv_spls <- function(X_sub, Y, folds, ncomp = 1, keepX) {
  results <- data.frame(R2 = NA, Spearman = NA, RMSE = NA)
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(Y), test_idx)
    
    X_train <- X_sub[train_idx, , drop = FALSE]
    y_train <- Y[train_idx]
    X_test  <- X_sub[test_idx, , drop = FALSE]
    y_test  <- Y[test_idx]
    
    # Fit sPLS on training set
    model <- spls(X_train, y_train, ncomp = ncomp, keepX = keepX)
    
    # Predict Y on test samples
    preds <- predict(model, X_test)
    y_pred <- preds$predict[, 1, 1]   # first Y, first component
    
    # Evaluate metrics in test data
    lm_fit <- lm(y_test ~ y_pred)
    R2  <- summary(lm_fit)$r.squared
    rho <- suppressWarnings(cor(y_test, y_pred, method = "spearman"))
    RMSE <- sqrt(mean((y_test - y_pred)^2))
    
    results[i, ] <- c(R2, rho, RMSE)
  }
  
  return(results)
}

cv_selected <- cv_spls(X_sub, Y, folds, ncomp = 1, keepX = rep(2569, 1))

summary(cv_selected)
mean(cv_selected$R2)
mean(cv_selected$Spearman)
mean(cv_selected$RMSE)


### Permutation testing
set.seed(500)

n_perm <- 100
perm_results <- numeric(n_perm)

for (p in 1:n_perm) {
  # Shuffle the outcome
  Y_perm <- sample(Y)
  
  # Run the same CV with permuted Y
  cv_perm <- cv_spls(X_sub, Y_perm, folds, ncomp = 1, keepX = rep(2569, 1))
  
  # Store mean R²
  perm_results[p] <- mean(cv_perm$R2)
  
  cat("Permutation", p, "done\n")
}

### Feature stability
set.seed(500)

folds <- createFolds(Y, k = 10, list = TRUE)

# Prepare matrix to record selections
var_selection <- matrix(0, nrow = ncol(X_sub), ncol = length(folds))
rownames(var_selection) <- colnames(X_sub)

keepX_value <- 50  

for (i in seq_along(folds)) {
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_along(Y), test_idx)
  
  X_train <- X_sub[train_idx, , drop = FALSE]
  y_train <- Y[train_idx]
  
  model <- spls(X_train, y_train, ncomp = 1, keepX = rep(keepX_value, 1))
  
  selected_vars <- selectVar(model, comp = 1)$X$name
  
  # Mark selected variables
  var_selection[selected_vars, i] <- 1
}

# Compute selection frequencies
var_freq <- rowSums(var_selection) / length(folds)

var_selection <- matrix(0, nrow = ncol(X_sub), ncol = length(folds))
rownames(var_selection) <- colnames(X_sub)

for (i in seq_along(folds)) {
  train_idx <- setdiff(seq_along(Y), folds[[i]])
  model <- spls(X_sub[train_idx, , drop = FALSE], Y[train_idx],
                ncomp = 1, keepX = rep(2569, 1))
  sel <- selectVar(model, comp = 1)$X$name
  if (!is.null(sel)) var_selection[sel, i] <- 1
}

# Compute frequency across folds
var_freq <- rowSums(var_selection) / length(folds)

# Create a dataframe
df_stability <- data.frame(
  Metabolite = names(var_freq),
  Stability = var_freq
)

### Extracting top X, rank and weight
# Loadings for the first component
loadings <- spls_final$loadings$X[, 1]

# Keep only nonzero loadings
nonzero_idx <- loadings != 0

# Creating a df for extraction
df_loadings <- data.frame(
  Metabolite = names(loadings)[nonzero_idx],
  Weight = loadings[nonzero_idx]
)

#Adding direction to the df
df_loadings$Direction <- ifelse(df_loadings$Weight > 0, "Positive", "Negative")

# Sort df by absolute weight and add rank
df_loadings <- df_loadings %>%
  arrange(desc(abs(Weight))) %>%
  mutate(Rank = row_number())

#Adding stability
df_loadings_full <- df_loadings %>%
  dplyr::left_join(df_stability, by = "Metabolite") %>%
  arrange(desc(abs(Weight)))


