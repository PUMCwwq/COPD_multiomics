##############################################################################
# Script information                                                      
# Title: adonis_analysis
# Author: Minyu Zhou
# Date: 2024-10-09
# Description: None
##############################################################################

library(dplyr)
library(vegan)
library(readxl)
library(tidyr)

# ===== Load data =====
cat("Loading data...\n")
SS <- read_excel("bacteria_rel.xlsx", sheet = "SS") %>% as.data.frame()
RS <- read_excel("bacteria_rel.xlsx", sheet = "RS") %>% as.data.frame()
CS <- read_excel("bacteria_rel.xlsx", sheet = "CS") %>% as.data.frame()
geo_name <- read_excel("meta.xlsx", sheet = 'geo_name') %>% as.data.frame()
num <- read_excel("meta.xlsx", sheet = 'num') %>% as.data.frame()
cat <- read_excel("meta.xlsx", sheet = 'cat') %>% as.data.frame()

# ===== Preprocess metadata =====
colnames(geo_name) <- gsub("\\W", "_", colnames(geo_name))
colnames(num) <- gsub("\\W", "_", colnames(num))
colnames(cat) <- gsub("\\W", "_", colnames(cat))

meta <- merge(geo_name, merge(num, cat, by = "SampleID"), by = "SampleID")
cat("Metadata merged:", nrow(meta), "rows,", ncol(meta), "columns\n")

# Clean sample names and set row names
meta <- meta %>% mutate(SampleID = trimws(as.character(SampleID)))
rownames(SS) <- SS$SampleID; SS <- SS[, -1]
rownames(RS) <- RS$SampleID; RS <- RS[, -1]
rownames(CS) <- CS$SampleID; CS <- CS[, -1]

# Group list
group_list <- list("SS" = SS, "RS" = RS, "CS" = CS)

# ===== Missing value check =====
cat("\n===== Missing values in metadata =====\n")
na_stats <- sapply(meta, function(x) sum(is.na(x)))
print(na_stats)

# ===== Define variables for adonis2 analysis =====
all_vars <- c("District", "Age", "BMI", "FEV1/FVC", "FEV1/pred", 
              "CAT score", "CS pack year", "Altitude", "EOS", 
              "NEUT", "LYMPH", "MONO", "Sex", 
              "Dust exposure", "GOLD", "AE", "TTFC")
vars2analyze <- all_vars[all_vars %in% colnames(meta)]

cat_vars <- c("District", "Sex", "Dust exposure", "GOLD", "AE", "TTFC")
num_vars <- setdiff(vars2analyze, cat_vars)

# Initialize results
Adonis_res <- NULL

# ===== Main loop for adonis2 analysis =====
for (df_name in names(group_list)) {
  cat("\n===== Processing dataset:", df_name, "=====\n")
  mf_df <- group_list[[df_name]] %>% as.data.frame()
  cat("Microbiome data dimensions:", dim(mf_df), "\n")
  
  for (v in vars2analyze) {
    cat("\n--- Variable:", v, "---\n")
    meta_sub <- meta %>% select(SampleID, all_of(v))
    
    if (v %in% cat_vars) {
      meta_sub <- meta_sub %>%
        mutate(!!sym(v) := as.character(!!sym(v))) %>%
        mutate(!!sym(v) := trimws(!!sym(v))) %>%
        filter(!!sym(v) != "" & !is.na(!!sym(v)))
      meta_sub[[v]] <- as.factor(meta_sub[[v]])
      
      lvl_counts <- table(meta_sub[[v]])
      min_samples <- 3
      if (any(lvl_counts < min_samples)) {
        small_levels <- names(lvl_counts[lvl_counts < min_samples])
        meta_sub <- meta_sub %>% filter(!(!!sym(v) %in% small_levels))
      }
      if(length(unique(meta_sub[[v]])) <= 1) next
      
    } else {
      meta_sub <- meta_sub %>%
        mutate(!!sym(v) := as.character(!!sym(v))) %>%
        mutate(!!sym(v) := trimws(!!sym(v))) %>%
        mutate(!!sym(v) := ifelse(!!sym(v) %in% c("NA", "", " ", "-", "?", "N/A", "#N/A"), NA, !!sym(v)))
      meta_sub[[v]] <- suppressWarnings(as.numeric(meta_sub[[v]]))
      meta_sub <- meta_sub %>% filter(!is.na(!!sym(v)))
      if(nrow(meta_sub) < 5) next
    }
    
    common_samples <- intersect(rownames(mf_df), meta_sub$SampleID)
    if (length(common_samples) < 5) next
    
    mf_df_sub <- mf_df[common_samples, ]
    meta_sub <- meta_sub %>% filter(SampleID %in% common_samples) %>%
      arrange(match(SampleID, common_samples))
    rownames(meta_sub) <- meta_sub$SampleID
    mf_df_sub <- mf_df_sub[match(meta_sub$SampleID, rownames(mf_df_sub)), ]
    
    zero_cols <- which(colSums(mf_df_sub) == 0)
    if(length(zero_cols) > 0) {
      mf_df_sub <- mf_df_sub[, -zero_cols]
    }
    if (any(is.na(mf_df_sub))) {
      mf_df_sub[is.na(mf_df_sub)] <- 0
    }
    
    fml <- as.formula(paste("mf_df_sub ~", v))
    
    adonis_result <- try(
      adonis2(
        fml, 
        data = meta_sub, 
        permutations = 999, 
        method = "bray", 
        by = "margin"
      ),
      silent = FALSE
    )
    
    if (inherits(adonis_result, "try-error")) {
      R2 <- NA; p <- NA; F_value <- NA
    } else {
      R2 <- adonis_result$R2[1]
      p <- adonis_result$`Pr(>F)`[1]
      F_value <- adonis_result$F[1]
    }
    
    res <- data.frame(
      "data" = df_name,
      "variable" = v,
      "adonis.r2" = R2,
      "adonis.F" = F_value,
      "adonis.p" = p,
      "sample_size" = length(common_samples)
    )
    
    Adonis_res <- bind_rows(Adonis_res, res)
  }
}

# ===== Output results =====
print(Adonis_res)

