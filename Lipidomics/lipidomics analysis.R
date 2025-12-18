# ==============================================================================
# Lipidomics Analysis 
# Author: Feargal Ryan
# Last Updated: 18/12/2025


# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)
library(readxl)
library(vsn)      # For variance-stabilizing normalization
library(ape)      # For principal coordinate analysis (PCoA)
library(limma)    # For linear modeling
library(dplyr)    # For data manipulation
library(ggrepel)

##########################################################################

# Load lipidomics data and metadata
lipid_counts <- as.data.frame(read_xlsx("Lipidomics data.xlsx", sheet = "Counts"))
rownames(lipid_counts) <- lipid_counts$FeatureID
lipid_counts <- lipid_counts[,-1]
lipid_counts <- lipid_counts[apply(is.na(lipid_counts), 1, sum) < 7, ]
lipid_counts <- as.data.frame(impute.knn(as.matrix(lipid_counts), k = 4)$data)


## Sample meta data
meta_data = as.data.frame(read_xlsx("Lipidomics data.xlsx",sheet = "Meta",)) ; rownames(meta_data) = meta_data$ParticipantID
meta_data = meta_data[meta_data$Cognitive_impairment != "Missing",]

## This is a table which contains which lipid family each species is from. 
taxa = as.data.frame(read_xlsx("Lipidomics data.xlsx",sheet = "Taxa",)) ; rownames(taxa) = taxa$FeatureID


#Ensure all the tables overlap
features = intersect(rownames(lipid_counts),rownames(taxa))
taxa = taxa[features,]
lipid_counts = lipid_counts[features,]


samples = intersect(colnames(lipid_counts),rownames(meta_data))
meta_data = meta_data[samples,]
lipid_counts = lipid_counts[,samples]

# ------------------------------------------------------------------------------
# 3. Normalize Data and Perform PCoA
# ------------------------------------------------------------------------------

# Normalize lipid counts using variance-stabilizing normalization (VSN)
normalized_counts <- normalizeVSN(lipid_counts)

# Plot mean-SD before and after normalization
vsn::meanSdPlot(as.matrix(lipid_counts))
vsn::meanSdPlot(normalized_counts)

# Perform MDS and PCoA for visualization
mds_out <- plotMDS(normalized_counts)
pcoa_out <- pcoa(as.dist(mds_out$distance.matrix))

# Create a dataframe for plotting
normy_data.df <- data.frame(meta_data, PC1 = mds_out$x, PC2 = mds_out$y)

# Plot MDS/PCoA results
ggplot(normy_data.df, aes(x = PC1, y = PC2, color = Cognitive_impairment)) +
  geom_point(size = 4) +
  theme_classic() +
  xlab(paste("PC1:", round(mds_out$var.explained[1] * 100, 2), "% variation")) +
  ylab(paste("PC2:", round(mds_out$var.explained[2] * 100, 2), "% variation")) +
  scale_color_manual(values = c("deepskyblue3", "firebrick3"))




# Load limma library (already loaded)
library(limma)

# Ensure your data is properly formatted (features as rows, samples as columns)
filtered_counts <- lipid_counts  # You can apply variance filtering here if needed

# Create the design matrix for your model
design <- model.matrix(~ Cognitive_impairment, data = meta_data)
# Fit the linear model for each feature
fit <- lmFit( normalized_counts, design)
# Empirical Bayes moderation
fit <- eBayes(fit)
# Extract p-values and fold changes (log2FC)
results <- topTable(fit, coef = "Cognitive_impairmentYes", number = Inf, adjust.method = "none")
results = data.frame(results,taxa[rownames(results),])

# Filter significant results (adjust this threshold as needed)
significant_results <- results[results$adj.P.Val < 0.1, ]

