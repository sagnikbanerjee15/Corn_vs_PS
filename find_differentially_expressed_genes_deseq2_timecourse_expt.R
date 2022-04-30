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
metadata = read.csv2(metadatafilename, sep = ",")
metadata$Replicate = NULL
metadata$Location = NULL
metadata$Ended = NULL

count_data_filename = args[2]
count_data = read.csv2(count_data_filename, sep = "\t", check.names = FALSE, nrows=1000000, stringsAsFactors=FALSE)
row.names(count_data) = count_data$Name
count_data$Name = NULL
df <- data.frame(apply(count_data, 2, function(x) as.numeric(as.character(x))))
count_data[, c(1,dim(count_data)[2])] <- apply(count_data[, c(1,dim(count_data)[2])],2,
                                                function(x) as.numeric(as.character(x))
)

count_data[, c(1,dim(count_data)[2])] <- apply(count_data[, c(1,dim(count_data)[2])],2,
                                                 function(x) trimws(x)
)

all <- apply(count_data, 1, function(x) all(x==0) )
count_data <- count_data[!all,]

output_filename_prefix = args[3]
############################################################################################################################################################################################################

cultivars = unique(metadata$Cultivar)
timepoints = unique(metadata$Timepoints)

############################################################################################################################################################################################################
# PCA plot to demonstrate replicate correlation
############################################################################################################################################################################################################

raw_counts = data.matrix(count_data)
raw_counts_pca = prcomp( raw_counts, center = TRUE, scale = TRUE)
raw_counts_pca = data.frame(raw_counts_pca$rotation[,1:2])
raw_counts_pca$SampleName = rownames(raw_counts_pca)
raw_counts_pca = merge(raw_counts_pca, metadata, by = 'SampleName')
#raw_counts_pca$Condition <- as.factor(as.numeric(as.factor(raw_counts_pca$Condition)))
raw_counts_pca
for(cultivar in cultivars)
{
  raw_counts_pca_cultivar = subset(raw_counts_pca, Cultivar == cultivar)
  ggplot(raw_counts_pca_cultivar, aes(x = PC1, y = PC2, fill = Timepoints, shape = Condition)) +
    geom_point(size = 2, stroke = 1, aes(color=Timepoints))
  
  output_filename = paste0(output_filename_prefix,"/","pca_raw_counts_",cultivar,".pdf")
  ggsave(output_filename, device = "pdf")
}

dds <- DESeqDataSetFromMatrix(countData = data.matrix(count_data),
                              colData = metadata,
                              design = ~ Condition)
dds = DESeq(dds)
normalized_counts = counts(dds, normalized = TRUE)
normalized_counts_pca = prcomp( normalized_counts, center = TRUE, scale = TRUE)
normalized_counts_pca = data.frame(normalized_counts_pca$rotation[,1:2])
normalized_counts_pca$SampleName = rownames(normalized_counts_pca)
normalized_counts_pca = merge(normalized_counts_pca, metadata, by = 'SampleName')
#raw_counts_pca$Condition <- as.factor(as.numeric(as.factor(raw_counts_pca$Condition)))
raw_counts_pca
for(cultivar in cultivars)
{
  normalized_counts_pca_cultivar = subset(normalized_counts_pca, Cultivar == cultivar)
  ggplot(normalized_counts_pca_cultivar, aes(x = PC1, y = PC2, fill = Timepoints, shape = Condition)) +
    geom_point(size = 2, stroke = 1, aes(color=Timepoints))
  
  output_filename = paste0(output_filename_prefix,"/","pca_normalized_counts_",cultivar,".pdf")
  ggsave(output_filename, device = "pdf")
}
############################################################################################################################################################################################################

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