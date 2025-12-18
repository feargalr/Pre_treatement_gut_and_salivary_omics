# ==============================================================================
# 16S Analysis 
# Author: Feargal Ryan
# Last Updated: 18/12/2025

# Load necessary libraries

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(phyloseq)
library(reshape)
library(ape)
library(vegan)
library(plyr)
library(scales)
library(readxl)



# ==============================================================================


## Reading in data
raw_counts = read.table("RSV_count_table.txt",sep="\t",header=TRUE,row.names = 1,check.names = FALSE) 
tree = read.tree("tree.nwk")
meta_data = as.data.frame(read_xlsx("meta_data.xlsx",sheet = 1)) ; 
taxa = read.delim("SILVA138_taxonomy.tsv",sep="\t",header=TRUE,row.names = 1) 

taxa[taxa==""] = "Unclassified" ; taxa[taxa=="uncultured_bacterium"] = "Unclassified" ; taxa[taxa=="uncultured"] = "Unclassified" ; 
taxa[taxa=="Incertae_Sedis"] = "Unclassified" ; taxa[grepl("_sp.",taxa$Species),] = "Unclassified" ; taxa[grepl("uncultured_",taxa$Species),] = "Unclassified" ;
taxa[grepl("mouse_gut",taxa$Species),] = "Unclassified" ; taxa[grepl("unidentified",taxa$Species),] = "Unclassified" ; taxa$Genus[grepl("GCA-",taxa$Genus)] = "Unclassified" ;

meta_data = meta_data[meta_data$Type=="Gut",]
meta_data = meta_data[meta_data$Cognitive_impairment!="Missing",]
meta_data = meta_data[meta_data$Diarrhoea!="Missing",]

samples2keep = intersect(meta_data$LibraryID,colnames(raw_counts))

raw_counts = raw_counts[,samples2keep]
meta_data = meta_data[meta_data$LibraryID %in% samples2keep,]
rownames(meta_data) = meta_data$SampleID
raw_counts = raw_counts[,meta_data$LibraryID]
colnames(raw_counts) = rownames(meta_data)


## Tidying up the taxonomic data 
taxa$Genus = as.character(taxa$Genus)
taxa$Genus = as.character(taxa$Genus)
taxa$Order = as.character(taxa$Order)


#Sets anything that is unclassified at the Genus level as it's Family, and if still unclassified to the order. 
taxa[taxa$Genus=="Unclassified","Genus"] = taxa[taxa$Genus=="Unclassified","Family"]
taxa[taxa$Genus=="Unclassified","Genus"] = taxa[taxa$Genus=="Unclassified","Order"]



raw_counts = raw_counts[rownames(taxa),]
taxa = taxa[rownames(raw_counts),]

genus_counts = aggregate(raw_counts, by=list(taxa$Genus), FUN=sum)
rownames(genus_counts) = genus_counts$Group.1
genus_counts = genus_counts[,-1]

meta_data = meta_data[colnames(genus_counts),]


genus_counts =genus_counts[apply(genus_counts>1000,1,sum)>2,]

#Convert relative abundance to a phyloseq object
OTU = otu_table(genus_counts,taxa_are_rows =TRUE)
map_data = sample_data(meta_data)
physeq = phyloseq(OTU,map_data)
ps_clr <- microbiome::transform(physeq, "clr")
Genus_clr = data.frame(otu_table(ps_clr))

#Differential abundance
CIPs = sort(apply(Genus_clr,1,function(x){kruskal.test(x,meta_data$Cognitive_impairment)$p.value}))
CIPs = p.adjust(CIPs,method="fdr")

data2plot = data.frame(meta_data,t(Genus_clr))

ggplot(data2plot,aes(x = Cognitive_impairment, y = Sutterella,fill=Cognitive_impairment))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.05)+theme_classic()+
  scale_fill_manual(values=c("deepskyblue3","firebrick4"))+guides(fill="none")+
  stat_compare_means(comparisons = list(c("No","Yes")))+
  ylab("Sutterella CLR")+
  theme(axis.title = element_text(size=10),axis.text = element_text(size=8,color="black"))

