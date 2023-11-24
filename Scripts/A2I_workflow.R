library(tidyverse)
library(stringr)

# Set working directory as the study folder
setwd("/sc/arion/projects/breen_lab/AD_editing/MSBB/")

# Load data ---------------------------------------------------------------
data <- read_csv("MSBB_aggregate_data.csv") %>%
  mutate(AD = factor(ifelse(CERAD %in% c(1, 2), "pos", ifelse(CERAD == 4, "neg", NA)))) 
data$BrodmannArea <- as.factor(data$BrodmannArea)
raw_edits <- read_tsv("./editing_sites/known_sites/edit_matrix_clean.tsv") # Table of edit sites and samples


# Create region-specific data frames
BM10 <- data %>%
  filter(BrodmannArea == 10)

# Create data frames for each region --------------------------------------
# Select column names from 'meta' for each BrodmannArea
BM10_samples <- data$filename[data$BrodmannArea == 10]

# Subset meta
BM10_aggregate <- aggregate_data %>%
  filter(BrodmannArea == 10)

# Subset raw edits into Brodmann Areas
BM10_raw <- raw_edits %>%
  select(site, all_of(BM10_samples)) %>%
  select(site, sort(names(.))) %>%
  column_to_rownames("site")
BM10_raw <- BM10_raw[, names(BM10_raw) %in% BM10$filename[!is.na(BM10$AD)]]


# AEI analysis ------------------------------------------------------------
library(ggplot2)

# AEI scores of all regions AD vs non-AD
sub_data <- data %>%
  filter(!is.na(AD))
ggplot(sub_data, aes(x = BrodmannArea, y = A2GEditingIndex, fill = AD)) +
  scale_fill_manual(values = c("#d80c8c", "#00aeef", "pink")) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Brodmann Area", y = "AEI Index", title = "AEIs of AD vs Controls") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Check statistical significance
# Separate the data into two groups: AD and Control
# Split the data into subsets based on 'BrodmannArea'
data_split <- split(sub_data, sub_data$BrodmannArea)

# Perform Mann-Whitney U test for each subset
test_results <- lapply(data_split, function(subset) {
  data_ad <- subset[subset$AD == "pos", ]
  data_control <- subset[subset$AD == "neg", ]
  result <- wilcox.test(data_ad$A2GEditingIndex, data_control$A2GEditingIndex)
  return(result)
})

# Print the test results for each subset
for (i in 1:length(test_results)) {
  cat("Brodmann Area:", names(test_results)[i], "\n")
  print(test_results[[i]])
  cat("\n")
}

# Combine the test results into a dataframe
results_df <- do.call(rbind, lapply(names(test_results), function(area) {
  result <- test_results[[area]]
  data.frame(
    BrodmannArea = area,
    p_value = result$p.value,
    test_statistic = result$statistic
  )
}))
results_df$BrodmannArea <- as.factor(results_df$BrodmannArea)

ggplot(sub_data, aes(x = BrodmannArea, y = A2GEditingIndex, fill = AD)) +
  scale_fill_manual(values = c("#d80c8c", "#00aeef")) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Brodmann Area", y = "AEI Index", title = "AEIs of AD vs Controls") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Add text labels for p-values
  geom_text(data = results_df, aes(x = BrodmannArea, label = paste("p =", round(p_value, 3))),
            vjust = -0.5, hjust = 0.5, size = 3) +
  # Add test statistic as a subtitle
  labs(subtitle = "Mann-Whitney U Test Statistic") +
  geom_text(data = results_df, aes(label = sprintf("U = %.2f", test_statistic)),
            vjust = -0.8, hjust = 0.5, size = 3)

# Visualize means of AEI 
BM10_AEI <- ggplot(BM10, aes(x = CERAD, y = A2GEditingIndex, group = CERAD)) +
  geom_violin() +
  labs(title = "A-to-I Indices")

# Visualize distribution of edited genic regions using grouped bar plot
genic_regions_data <- data %>%
  select("filename","AD","UTR3","UTR5","Downstream","Exonic","Intergenic","Intronic","ncRNA_exonic",
         "ncRNA_intronic","Splicing","Upstream","Upstream;Downstream") %>%
  pivot_longer(cols = c("UTR3","UTR5","Downstream","Exonic","Intergenic","Intronic","ncRNA_exonic",
                        "ncRNA_intronic","Splicing","Upstream","Upstream;Downstream"),
               names_to = "region",
               values_to = "value")

genic_regions_mean_se <- genic_regions_data %>%
  filter(!is.na(AD)) %>%
  group_by(AD, region) %>%
  summarize(mean_value = mean(value),
            sd_value = sd(value))

ggplot(genic_regions_mean_se, aes(x = region, y = mean_value, fill = AD)) +
  scale_fill_manual(values = c("#d80c8c", "#00aeef", "pink")) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                position = position_dodge(width = 0.8), width = 0.25) +
  theme_classic() +
  labs(x = "Genic Region", y = "Proportion of Edits", title = "Proportion of Edits in Various Genic Regions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Differential editing analysis -------------------------------------------
# Load in txt file that contains metadata for all editing sites
editing_meta <- read_tsv("/sc/arion/projects/breen_lab/OCD_RNAseq/CNS_REDi_combined.txt") %>%
  mutate(site = paste(chromosome, position, sep = "_")) %>%
  relocate(site, .before = 1) %>%
  select(site, Repeat, Region, Gene) %>%
  distinct(site, .keep_all = TRUE) 

library(limma)
library(edgeR)

CDR = as.factor(BM10$CDR)
Braak = as.factor(BM10$Braak)
CERAD = as.factor(BM10$CERAD)
AD = as.factor(BM10$AD)
Ancestry = as.factor(BM10$ethnicity)
Sex = as.factor(BM10$sex)
RIN = (BM10$RIN)

design <- model.matrix(~0+AD+Sex+Ancestry+RIN) #Create design matrix
fit <- lmFit(BM10_raw,design) # Fit linear model
cm <-makeContrasts(DevEffect = (ADpos - ADneg),levels=design)
fit2 <- contrasts.fit(fit, cm)
fitDupCor <- eBayes(fit2)
DE_sites<- topTable(fitDupCor, coef="DevEffect", n=nrow(BM10_raw)) 
DE_sites <- DE_sites %>% mutate(neglog10.adj.P.Val = -log10(adj.P.Val))

ggplot(DE_sites, aes(x = logFC, y = neglog10.adj.P.Val)) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC) > 0.15, "red", "black")), size = 2) +
  scale_color_identity() +
  labs(x = "logFC", y = "adj.P.Val") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.15, 0.15), linetype = "dashed", color = "blue") +
  theme_minimal()

# Look at AEI in AD vs non-AD
# Visualize means of AEI 
BM10_AEI <- ggplot(BM10, aes(x = CERAD, y = A2GEditingIndex, group = CERAD)) +
  geom_violin() +
  labs(title = "A-to-I Indices")
