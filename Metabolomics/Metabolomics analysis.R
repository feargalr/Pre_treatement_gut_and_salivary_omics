# ==============================================================================
# Metabolomics Analysis 
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

##########################################################################

## metabolite counts
metabolite_counts = as.data.frame(read_xlsx("Metabolomics data.xlsx",sheet = "Counts")) ; rownames(metabolite_counts) = metabolite_counts$FeatureID ; metabolite_counts = metabolite_counts[,-1]


## Sample meta data
meta_data = as.data.frame(read_xlsx("Metabolomics data.xlsx",sheet = "Meta",)) ; rownames(meta_data) = meta_data$ParticipantID
meta_data = meta_data[meta_data$Cognitive_impairment != "Missing",]

## This is a table which contains which metabolite family each species is from. 
taxa = as.data.frame(read_xlsx("Metabolomics data.xlsx",sheet = "Taxa",)) ; rownames(taxa) = taxa$FeatureID
taxa = taxa[!grepl("Metabolite_",taxa$Name),]

# Ensure that the filtered counts match with the sample metadata
features <- intersect(rownames(metabolite_counts), rownames(taxa))
taxa <- taxa[features, ]
metabolite_counts <- metabolite_counts[features, ]

# Ensure that sample names still overlap between counts and metadata
samples <- intersect(colnames(metabolite_counts), rownames(meta_data))
meta_data <- meta_data[samples, ]
metabolite_counts <- metabolite_counts[, samples]


#Ensure all the tables overlap
features = intersect(rownames(metabolite_counts),rownames(taxa))
taxa = taxa[features,]
metabolite_counts = metabolite_counts[features,]


samples = intersect(colnames(metabolite_counts),rownames(meta_data))
meta_data = meta_data[samples,]
metabolite_counts = metabolite_counts[,samples]

# ------------------------------------------------------------------------------
# 3. Normalize Data and Perform PCoA
# ------------------------------------------------------------------------------

# Normalize metabolite counts using variance-stabilizing normalization (VSN)
normalized_counts <- normalizeVSN(metabolite_counts)

# Plot mean-SD before and after normalization
vsn::meanSdPlot(as.matrix(metabolite_counts))
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


# Load necessary libraries
library(umap)
library(ggplot2)

# Perform UMAP on normalized counts
umap_out <- umap(t(normalized_counts)) # Transpose to have samples as rows, features as columns

# Create a dataframe for plotting
normy_data.df <- data.frame(meta_data, UMAP1 = umap_out$layout[, 1], UMAP2 = umap_out$layout[, 2])

# Plot UMAP results
ggplot(normy_data.df, aes(x = UMAP1, y = UMAP2, color = Cognitive_impairment)) +
  geom_point(size = 4) +
  theme_classic() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_color_manual(values = c("deepskyblue3", "firebrick3"))

# Prepare data for differential analysis
all_data.df <- data.frame(meta_data, t(normalized_counts))
all_data.df$Cognitive_impairment = 
  factor(all_data.df$Cognitive_impairment,levels=c("No","Yes"))
p_values <- c()
fold_changes <- c()


# Load limma library (already loaded)
library(limma)

# Ensure your data is properly formatted (features as rows, samples as columns)
metabolite_counts <- metabolite_counts  # You can apply variance filtering here if needed

# Create the design matrix for your model
design <- model.matrix(~ Cognitive_impairment, data = meta_data)
# Fit the linear model for each feature
fit <- lmFit(normalized_counts, design)
# Empirical Bayes moderation
fit <- eBayes(fit)
# Extract p-values and fold changes (log2FC)
results_cognitive_impairment <- topTable(fit, coef = "Cognitive_impairmentYes", number = Inf, adjust.method = "none")
results_cognitive_impairment = data.frame(results_cognitive_impairment,taxa[rownames(results_cognitive_impairment),])
# Filter significant results (adjust this threshold as needed)
significant_results_cognitive_impairment <- results_cognitive_impairment[results_cognitive_impairment$adj.P.Val < 0.1, ]

data2plot = data.frame(meta_data,t(normalized_counts))

ggplot(data2plot,aes(x=Cognitive_impairment,y=Metabolite_38,color=Cognitive_impairment))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.1)+
  stat_compare_means(label="p.signif",comparisons = list(c("No","Yes")))+
  theme_classic()+
  scale_color_manual(values=c("No"="deepskyblue3","Yes"="firebrick3"))+
  ylab("6-Hydroxynicotinic acid")

grid.arrange(
ggplot(data2plot,aes(x=Cognitive_impairment,y=Metabolite_38,color=Cognitive_impairment))+
  geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.06)+theme_classic()+
  xlab("Cognitive impairment")+scale_color_manual(values=c("No"="deepskyblue4","Yes"="firebrick4"))+
  ylab(taxa["Metabolite_38",2])+guides(color="none")+
  stat_compare_means(comparisons = list(c("No","Yes"))),
ggplot(data2plot,aes(x=Cognitive_impairment,y=Metabolite_232,color=Cognitive_impairment))+
  geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.06)+theme_classic()+
  xlab("Cognitive impairment")+scale_color_manual(values=c("No"="deepskyblue4","Yes"="firebrick4"))+
  ylab(taxa["Metabolite_232",2])+guides(color="none")+
  stat_compare_means(comparisons = list(c("No","Yes"))),
ggplot(data2plot,aes(x=Cognitive_impairment,y=Metabolite_105,color=Cognitive_impairment))+
  geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.06)+theme_classic()+
  xlab("Cognitive impairment")+scale_color_manual(values=c("No"="deepskyblue4","Yes"="firebrick4"))+
  ylab(taxa["Metabolite_105",2])+guides(color="none")+
  stat_compare_means(comparisons = list(c("No","Yes"))),
ggplot(data2plot,aes(x=Cognitive_impairment,y=Metabolite_36,color=Cognitive_impairment))+
  geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.06)+theme_classic()+
  xlab("Cognitive impairment")+scale_color_manual(values=c("No"="deepskyblue4","Yes"="firebrick4"))+
  ylab(taxa["Metabolite_36",2])+guides(color="none")+
  stat_compare_means(comparisons = list(c("No","Yes"))),nrow=1)