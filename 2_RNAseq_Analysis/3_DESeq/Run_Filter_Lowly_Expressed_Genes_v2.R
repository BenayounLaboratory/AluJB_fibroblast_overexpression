# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Process_RNASeq_functions.R")

# Load libraries
library(DESeq2)


    
# Set 1: T24 IMR-90 AluJb OE ------------------------------------------------------------------------------------------------------------------------
 

  
# Define the output directory
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/'

# Define the directory with the counts tables and collect individual file paths in a list.
count_Table.dir <- c("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/2_Read_counting/counts_T24_AluJb_IMR90_Normal_Media/")
cntTable_list <- list.files(count_Table.dir, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables into one df.
gene_TE_counts_raw <- aggregate_counts(cntTable_complete.list = cntTable_list)

# Remove useless info. from sample names
colnames(gene_TE_counts_raw) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub("*fastp_", "", colnames(gene_TE_counts_raw))

# Reorder columns so controls are first
gene_TE_counts_raw <- gene_TE_counts_raw[, c(1, 6:9, 2:5)]

# Check, but Remove plasmid gene features to exclude from downstream analysis
  
    # Define the gene features
    plasmid_gene_features_to_exclude <- c('U6_AluJb')
    
    # Obtain gene feature indices in the matrix
    plasmid_gene_feature_indices <- which(gene_TE_counts_raw$gene.TE %in% plasmid_gene_features_to_exclude)
    
    # Check the counts for these features
    gene_TE_counts_raw[plasmid_gene_feature_indices, ]
    
    # Remove these features from the counts table
    gene_TE_counts_raw <- gene_TE_counts_raw[-plasmid_gene_feature_indices, ]
      
# Run function to remove gene/transcript version info, combine transcripts from the same gene, and filter lowly expressed genes *IMPORTANT: PSEUDOAUTOSOMAL GENES (PAR_Y) ARE COMBINED
counts_filtered <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 4/8)

    # Save counts files.
    write.table(counts_filtered, file = paste(dir.output, "Filtered_counts_T24_AluJb_OE", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)




    
# Set 2: T72 IMR-90 AluJb OE for Proteomics ------------------------------------------------------------------------------------------------------------------------
 

  
# Define the output directory
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/'

# Define the directory with the counts tables and collect individual file paths in a list.
count_Table.dir <- c("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/2_Read_counting/counts_T72_AluJb_IMR90_for_Proteomics/")
cntTable_list <- list.files(count_Table.dir, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables into one df.
gene_TE_counts_raw <- aggregate_counts(cntTable_complete.list = cntTable_list)

# Remove useless info. from sample names
colnames(gene_TE_counts_raw) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub("*fastp_", "", colnames(gene_TE_counts_raw))

# Reorder columns so controls are first
gene_TE_counts_raw <- gene_TE_counts_raw[, c(1, 6:9, 2:5)]

# Check, but Remove plasmid gene features to exclude from downstream analysis
  
    # Define the gene features
    plasmid_gene_features_to_exclude <- c('U6_AluJb')
    
    # Obtain gene feature indices in the matrix
    plasmid_gene_feature_indices <- which(gene_TE_counts_raw$gene.TE %in% plasmid_gene_features_to_exclude)
    
    # Check the counts for these features
    gene_TE_counts_raw[plasmid_gene_feature_indices, ]
    
    # Remove these features from the counts table
    gene_TE_counts_raw <- gene_TE_counts_raw[-plasmid_gene_feature_indices, ]
      
# Run function to remove gene/transcript version info, combine transcripts from the same gene, and filter lowly expressed genes *IMPORTANT: PSEUDOAUTOSOMAL GENES (PAR_Y) ARE COMBINED
counts_filtered <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 4/8)

    # Save counts files.
    write.table(counts_filtered, file = paste(dir.output, "Filtered_counts_T72_Proteomics", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)




    

# Set 3: T72 IMR-90 Lamivudine ------------------------------------------------------------------------------------------------------------------------
 

  
# Define the output directory
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/'

# Define the directory with the counts tables and collect individual file paths in a list.
count_Table.dir <- c("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/2_Read_counting/counts_T72_AluJb_IMR90_Normal_Media_3TC/")
cntTable_list <- list.files(count_Table.dir, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables into one df.
gene_TE_counts_raw <- aggregate_counts(cntTable_complete.list = cntTable_list)

# Remove useless info. from sample names
colnames(gene_TE_counts_raw) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub("*fastp_", "", colnames(gene_TE_counts_raw))

# Reorder columns so controls are first
gene_TE_counts_raw <- gene_TE_counts_raw[, c(1, 13:16, 10:12, 6:9, 2:5)]

# Check, but Remove plasmid gene features to exclude from downstream analysis
  
    # Define the gene features
    plasmid_gene_features_to_exclude <- c('U6_AluJb')
    
    # Obtain gene feature indices in the matrix
    plasmid_gene_feature_indices <- which(gene_TE_counts_raw$gene.TE %in% plasmid_gene_features_to_exclude)
    
    # Check the counts for these features
    gene_TE_counts_raw[plasmid_gene_feature_indices, ]
    
    # Remove these features from the counts table
    gene_TE_counts_raw <- gene_TE_counts_raw[-plasmid_gene_feature_indices, ]
      
# Run function to remove gene/transcript version info, combine transcripts from the same gene, and filter lowly expressed genes *IMPORTANT: PSEUDOAUTOSOMAL GENES (PAR_Y) ARE COMBINED
counts_filtered <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 3/15)

    # Save counts files.
    write.table(counts_filtered, file = paste(dir.output, "Filtered_counts_T72_Lamivudine", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)




    

# Set 4: T72 WI-38 AluJb OE ------------------------------------------------------------------------------------------------------------------------
 

  
# Define the output directory
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/'

# Define the directory with the counts tables and collect individual file paths in a list.
count_Table.dir <- c("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/2_Read_counting/counts_T72_AluJb_WI38_Normal_Media/")
cntTable_list <- list.files(count_Table.dir, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables into one df.
gene_TE_counts_raw <- aggregate_counts(cntTable_complete.list = cntTable_list)

# Remove useless info. from sample names
colnames(gene_TE_counts_raw) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub("*fastp_", "", colnames(gene_TE_counts_raw))

# Reorder columns (To get U6_1...U6_4, Alu_1...Alu_4)
gene_TE_counts_raw <- gene_TE_counts_raw[, c('gene.TE', 'JB08', 'JB07', 'JB05', 'JB03', 'JB04', 'JB01', 'JB02', 'JB06')]

# Update colnames
colnames(gene_TE_counts_raw) <- c('gene.TE', 'U6_1', 'U6_2', 'U6_3', 'U6_4', 'Alu_1', 'Alu_2', 'Alu_3', 'Alu_4')

# Check, but Remove plasmid gene features to exclude from downstream analysis
  
    # Define the gene features
    plasmid_gene_features_to_exclude <- c('U6_AluJb')
    
    # Obtain gene feature indices in the matrix
    plasmid_gene_feature_indices <- which(gene_TE_counts_raw$gene.TE %in% plasmid_gene_features_to_exclude)
    
    # Check the counts for these features
    gene_TE_counts_raw[plasmid_gene_feature_indices, ]
    
    # Remove these features from the counts table
    gene_TE_counts_raw <- gene_TE_counts_raw[-plasmid_gene_feature_indices, ]
      
# Run function to remove gene/transcript version info, combine transcripts from the same gene, and filter lowly expressed genes *IMPORTANT: PSEUDOAUTOSOMAL GENES (PAR_Y) ARE COMBINED
counts_filtered <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 4/8)

    # Save counts files.
    write.table(counts_filtered, file = paste(dir.output, "Filtered_counts_T72_WI38", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)







# Environment cleanup ------------------------------------------------------------------------------------------------------------------------
    

# Clean the environment
rm(list=ls())
    
    