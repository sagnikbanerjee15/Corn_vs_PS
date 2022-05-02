#! /usr/bin/env Rscript

############################################################################################################################################################################################################
# Program finds differentially expressed genes in a timecourse study

# Usage: find_differentially_expressed_genes_deseq2_timecourse_expt <metadatafilename> <countsfilename> <outputfileanme_prefix>
############################################################################################################################################################################################################

############################################################################################################################################################################################################
# Load libraries
############################################################################################################################################################################################################
library("DESeq2")
library("ggplot2")
library("dplyr")
library("reshape2")
library("data.table")

args = commandArgs(trailingOnly=TRUE)
print(args)
############################################################################################################################################################################################################

############################################################################################################################################################################################################
# Read metadata
############################################################################################################################################################################################################

if (length(args)!=3) 
{
  stop("You must provide exactly 3 arguments", call.=FALSE)
} 

metadatafilename = args[1]
metadatafilename = "/Users/sagnik/work/Corn_vs_PS/metadata.csv"
metadata = read.csv2(metadatafilename, sep = ",")
metadata$Replicate = NULL
metadata$Location = NULL
metadata$Ended = NULL
metadata$Cultivar = as.factor(metadata$Cultivar)
metadata$Condition = as.factor(metadata$Condition)
metadata$Timepoints = as.factor(metadata$Timepoints)

count_data_filename = args[2]
count_data_filename = "/Users/sagnik/work/Corn_vs_PS_data/analysis_with_argentinian_isolate/counts_across_all_samples.csv"
count_data = read.csv2(count_data_filename, sep = "\t", check.names = FALSE, nrows=10000000, stringsAsFactors=FALSE)
row.names(count_data) = count_data$Name
count_data$Name = NULL
count_data[, c(1,dim(count_data)[2])] <- apply(count_data[, c(1,dim(count_data)[2])],2,
                                                function(x) as.numeric(as.character(x))
)
column_names = colnames(count_data)
row_names = rownames(count_data)
count_data = data.frame(lapply(count_data, as.numeric))
colnames(count_data) = column_names
rownames(count_data) = row_names
count_data = round(count_data)

all <- apply(count_data, 1, function(x) all(x==0) )
count_data <- count_data[!all,]

output_filename_prefix = args[3]
output_filename_prefix = "/Users/sagnik/work/Corn_vs_PS_data/analysis_with_argentinian_isolate/deseq2_results/timecourse"
############################################################################################################################################################################################################

cultivars = unique(metadata$Cultivar)
timepoints = unique(metadata$Timepoints)

############################################################################################################################################################################################################
# PCA plots to demonstrate replicate correlation
############################################################################################################################################################################################################
gm_mean = function(x, na.rm=TRUE)
{
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

generatePCAPlotsToDetermineBiologicalReplicateSimilarity <- function(counts, prefix, cultivars)
{
  flag = 0
  for(cultivar in cultivars)
  {
    if(file.exists(paste0(prefix,"_",cultivar,".pdf"))==FALSE)
    {
      flag = 1
      break
    }
  }
  if(flag == 1)
  {
    #counts = count_data
    for(cultivar in cultivars)
    {
      #cultivar = "H95"
      metadata_cultivar = subset(metadata, Cultivar == cultivar)
      counts_cultivar = counts %>% select(metadata_cultivar$SampleName)
      
      counts_cultivar_pca = prcomp( counts_cultivar, center = TRUE, scale = TRUE)
      counts_cultivar_pca = data.frame(counts_cultivar_pca$rotation[,1:2])
      
      counts_cultivar_pca$SampleName = rownames(counts_cultivar_pca)
      counts_cultivar_pca = merge(counts_cultivar_pca, metadata, by = "SampleName")
      print(prefix)
      print(cultivar)
      print(counts_cultivar_pca)
      ggplot(counts_cultivar_pca, aes(x = PC1, y = PC2, fill = Timepoints, shape = Condition)) +
        geom_point(size = 2, stroke = 1, aes(color=Timepoints))
      
      output_filename = paste0(prefix,"_",cultivar,".pdf")
      ggsave(output_filename, device = "pdf")
    }
  }
}

# Generate PCA plots for raw counts
generatePCAPlotsToDetermineBiologicalReplicateSimilarity(counts = count_data, prefix = paste0(output_filename_prefix,"/","pca_raw_counts"), cultivars = cultivars)

# Generate PCA plots for DESeq2 normalized counts

dds_entire_dataset <- DESeqDataSetFromMatrix(countData = round(count_data),
                              colData = metadata,
                              design = ~ Condition )
dds_entire_dataset = estimateSizeFactors(dds_entire_dataset)
sizeFactors(dds_entire_dataset)
sizeFactors(dds_entire_dataset) = colSums(count_data)/gm_mean(colSums(count_data))
sizeFactors(dds_entire_dataset)
normalized_counts = as.data.frame(counts(dds_entire_dataset, normalized = TRUE))
generatePCAPlotsToDetermineBiologicalReplicateSimilarity(counts = normalized_counts, prefix = paste0(output_filename_prefix,"/","pca_deseq2_normalized_counts"), cultivars = cultivars)

cor(count_data$`1`,count_data$`3`, method = "pearson")
cor(normalized_counts$`1`,normalized_counts$`3`, method = "pearson")
############################################################################################################################################################################################################

############################################################################################################################################################################################################
# Correlation heatmaps to demonstrate replicate correlation
############################################################################################################################################################################################################

get_lower_tri<-function(cormat)
{
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat)
{
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

generateHeatmapsToDetermineBiologicalReplicateSimilarity <- function(counts, prefix, cultivars)
{
  flag = 0
  for(cultivar in cultivars)
  {
    if(file.exists(paste0(prefix,"_",cultivar,".pdf"))==FALSE)
    {
      flag = 1
      break
    }
  }
  if(flag == 1)
  {
    for(cultivar in cultivars)
    {
      metadata_cultivar = subset(metadata, Cultivar == cultivar)
      count_data_cultivar = data.matrix(select(counts, metadata_cultivar$SampleName))
      counts_corr = cor(data.matrix(count_data_cultivar))
      #counts_corr = get_upper_tri(counts_corr)
      melted_counts_corr = melt(counts_corr)
      print(prefix)
      print(cultivar)
      print(head(melted_counts_corr))
      melted_counts_corr$Var1 = as.factor(melted_counts_corr$Var1)
      melted_counts_corr$Var2 = as.factor(melted_counts_corr$Var2)
      melted_counts_corr = merge(melted_counts_corr, metadata_cultivar, by.x = c("Var1"), by.y = c("SampleName"))
      melted_counts_corr$labels1 = paste0(melted_counts_corr$Condition,"_" , melted_counts_corr$Timepoints, "_", melted_counts_corr$Var1)
      melted_counts_corr$Cultivar = NULL
      melted_counts_corr$Condition = NULL
      melted_counts_corr$Timepoints = NULL
      melted_counts_corr = merge(melted_counts_corr, metadata_cultivar, by.x = c("Var2"), by.y = c("SampleName"))
      melted_counts_corr$labels2 = paste0(melted_counts_corr$Condition,"_" , melted_counts_corr$Timepoints, "_", melted_counts_corr$Var2)
      melted_counts_corr$Cultivar = NULL
      melted_counts_corr$Condition = NULL
      melted_counts_corr$Timepoints = NULL
      
      ggplot(data = melted_counts_corr, aes(x=labels1, y=labels2, fill=value)) + 
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                              midpoint = 0.5, limit = c(0,1), space = "Lab", 
                              )+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 90))
      
      output_filename = paste0(prefix,"_",cultivar,".pdf")
      ggsave(output_filename, device = "pdf")
    }
  }
}

generateHeatmapsToDetermineBiologicalReplicateSimilarity(counts = count_data, prefix = paste0(output_filename_prefix,"/","correlation_raw_counts"), cultivars = cultivars)
generateHeatmapsToDetermineBiologicalReplicateSimilarity(counts = as.data.frame(normalized_counts), prefix = paste0(output_filename_prefix,"/","correlation_deseq2_normalized_counts"), cultivars = cultivars)

############################################################################################################################################################################################################
# DESeq2 analysis - 
# 1. Compare between inoculated and control for each timepoint and each cultivar
############################################################################################################################################################################################################

for(cultivar in cultivars)
{
  for(timepoint in timepoints)
  {
    metadata_cultivar_timepoint = subset(metadata, Cultivar == cultivar & Timepoints == timepoint)
    count_data_cultivar_timepoint = count_data[,metadata_cultivar_timepoint$SampleName]
    #print(colnames(count_data_cultivar_timepoint))
    #print(count_data_cultivar_timepoint)
    #print(data.matrix(count_data_cultivar_timepoint))
    dds = DESeqDataSetFromMatrix(countData = data.matrix(count_data_cultivar_timepoint),
                                  colData = metadata_cultivar_timepoint,
                                  design = ~ Condition)
    dds = DESeq(dds)
    print(resultsNames(dds))
    output_filename = paste0(output_filename_prefix,"/","each_cultivar_each_timepoint_",cultivar,"_",timepoint,".csv")
    write.csv(results(dds),output_filename)
    break
  }
  break
}

############################################################################################################################################################################################################
# DESeq2 analysis - 
# 2. Compare average between inoculated and control for each cultivar (all timepoints averaged)
############################################################################################################################################################################################################