library(tidyverse)
library(stringr)

# Set working directory as the study folder
setwd("/sc/arion/projects/breen_lab/AD_editing/MSBB/")

# Load data ---------------------------------------------------------------
raw_edits <- read_tsv("edit_matrix_clean.tsv") # Table of edit sites and samples
summarized_edits <- read_tsv("summarized_edits.tsv") # Edit regions and elements
AEI <- read_csv("../../AEI/total_AEI_clean.csv") # AEI scores
meta <- read_tsv("../../MSBB_meta") # Metadata
aggregate_data <- inner_join(meta, summarized_edits, by = "filename") %>%
  left_join(AEI, by = "filename") # Combine dataframes into mega table
write.csv(aggregate_data, file = "/sc/arion/projects/breen_lab/AD_editing/MSBB/MSBB_aggregate_data.csv", row.names = FALSE)

data <- read_csv("MSBB_aggregate_data.csv")

# Create data frames for each region --------------------------------------
# Select column names from 'meta' for each BrodmannArea
BM10_samples <- meta$filename[meta$BrodmannArea == 10]
BM22_samples <- meta$filename[meta$BrodmannArea == 22]
BM36_samples <- meta$filename[meta$BrodmannArea == 36]
BM44_samples <- meta$filename[meta$BrodmannArea == 44]

# Subset meta
BM10_aggregate <- aggregate_data %>%
  filter(BrodmannArea == 10)
BM22_aggregate <- aggregate_data %>%
  filter(BrodmannArea == 22)
BM36_aggregate <- aggregate_data %>%
  filter(BrodmannArea == 36)
BM44_aggregate <- aggregate_data %>%
  filter(BrodmannArea == 44)

# Subset raw edits into Brodmann Areas
BM10_raw <- raw_edits %>%
  select(site, all_of(BM10_samples)) %>%
  select(site, sort(names(.))) %>%
  column_to_rownames("site")
BM22_raw <- raw_edits %>%
  select(site, all_of(BM22_samples)) %>%
  select(site, sort(names(.))) %>%
  column_to_rownames("site")
BM36_raw <- raw_edits %>%
  select(site, all_of(BM36_samples)) %>%
  select(site, sort(names(.))) %>%
  column_to_rownames("site")
BM44_raw <- raw_edits %>%
  select(site, all_of(BM44_samples)) %>%
  select(site, sort(names(.))) %>%
  column_to_rownames("site")

# AEI analysis ------------------------------------------------------------
library(ggplot2)

# Visualize means of AEI 
BM10_AEI <- ggplot(BM10_aggregate, aes(x = CERAD, y = A2GEditingIndex, group = CERAD)) +
  geom_violin() +
  labs(title = "A-to-I Indices")

# Visualize genic region




# Differential editing analysis -------------------------------------------
# Load in txt file that contains metadata for all editing sites
editing_meta <- read_tsv("/sc/arion/projects/breen_lab/OCD_RNAseq/CNS_REDi_combined.txt") %>%
  mutate(site = paste(chromosome, position, sep = "_")) %>%
  relocate(site, .before = 1) %>%
  select(site, Repeat, Region, Gene) %>%
  distinct(site, .keep_all = TRUE) 

library(limma)
library(edgeR)

CDR = as.factor(BM10_aggregate$CDR)
Braak = as.factor(BM10_aggregate$Braak)
CERAD = as.factor(BM10_aggregate$CERAD)
Ancestry = as.factor(BM10_aggregate$ethnicity)
Sex = as.factor(BM10_aggregate$sex)
RIN = (BM10_aggregate$RIN)

design <- model.matrix(~0+CDR+Sex+Ancestry+RIN) #Create design matrix
fit <- lmFit(BM10_raw,design) # Fit linear model
cm <-makeContrasts(DevEffect = (CDR5 - CDR0),levels=design)
fit2 <- contrasts.fit(fit, cm)
fitDupCor <- eBayes(fit2)
DE_sites<- topTable(fitDupCor, coef="DevEffect", n=nrow(BM10_raw))

topTable(fitDupCor, coef="DevEffect")
DE_sites<- topTable(fitDupCor, coef="DevEffect", n=nrow(GeneExprs))
write.table(DE_sites, "DEG_BrainVar_sites.txt", sep="\t")
