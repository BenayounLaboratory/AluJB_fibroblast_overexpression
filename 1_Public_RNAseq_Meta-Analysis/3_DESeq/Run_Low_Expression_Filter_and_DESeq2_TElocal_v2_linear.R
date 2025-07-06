# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/Process_RNASeq_functions.R")

# Load libraries
library(limma) # Correcting for batch effects
library(DESeq2) # For differential expression
library(BiocParallel) # Used by DESeq2 for parallelization
  register(MulticoreParam(6)) # Use six cores
library(biomaRt) # To map Ensembl IDs to gene symbols
library(DEGreport) # only needed for LRT
library(pheatmap)
library(genefilter)
library(splitstackshape)

# Define the output directories
counts.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/Processed_counts/'
DESeq.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/DESeq_Results/'
MDS.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/MDS_PCA/'
Heatmaps.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/Heatmaps/'

# Load sample metadata
my.meta_data <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/0_Metadata_and_External_Resources/GSE113957_Metadata_Aging_Subset2_White-Caucasian_ArmSkin_CLEAN.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)


    
# FILTER LOW EXPRESSION GENES ------------------------------------------------------------------------------------------------------------------------
  

# Define the directory with the counts tables and collect individual file paths in a list.
count_Table.dir <- c("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/2_Read_counting/TElocal_counts_GSE113957_Human_Dermal_Fibroblasts/")
cntTable_list <- list.files(count_Table.dir, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables into one df.
gene_TE_counts_raw <- aggregate_counts(cntTable_complete.list = cntTable_list)

# Remove genes/TEs with all zeros
gene_TE_counts_raw <- gene_TE_counts_raw[which(rowSums(gene_TE_counts_raw[, -c(1)]) > 0), ]

# Remove useless info. from sample names
colnames(gene_TE_counts_raw) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub(".sra*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub("*fastp_", "", colnames(gene_TE_counts_raw))

# Keep only the final set of samples that will be analyzed (i.e. that are in the cleaned up meta-data file)
gene_TE_counts_raw <- gene_TE_counts_raw[, c('gene.TE', rownames(my.meta_data))]

# Load repeat subsets
repeats.intergenic.nearby <- read.csv("/Users/juanb/Library/CloudStorage/Dropbox/Research/2_Shared_Bioinformatic_Resources/TETranscripts_TElocal_annotation_files/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE44/Repeats_intergenic_within_5kb.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
repeats.intergenic.distal <- read.csv("/Users/juanb/Library/CloudStorage/Dropbox/Research/2_Shared_Bioinformatic_Resources/TETranscripts_TElocal_annotation_files/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE44/Repeats_intergenic_farther_than_5kb.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
repeats.intronic <- read.csv("/Users/juanb/Library/CloudStorage/Dropbox/Research/2_Shared_Bioinformatic_Resources/TETranscripts_TElocal_annotation_files/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE44/Repeats_intronic.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
repeats.exonic <- read.csv("/Users/juanb/Library/CloudStorage/Dropbox/Research/2_Shared_Bioinformatic_Resources/TETranscripts_TElocal_annotation_files/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE44/Repeats_exonic.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)

    # Assign repeat_name to rownames
    rownames(repeats.intergenic.nearby) <- repeats.intergenic.nearby$V9
    rownames(repeats.intergenic.distal) <- repeats.intergenic.distal$V9
    rownames(repeats.intronic) <- repeats.intronic$V9
    rownames(repeats.exonic) <- repeats.exonic$V9

# Combine TE loci by subfamily, but also splitting depending on location relative to genes (intergenic/exonic/intronic)
# Function gives a warning message that does not impact results:
# In type.convert.default(X[[i]], ...) :
# 'as.is' should be specified by the caller; using TRUE
gene_TE_counts_raw.aggregate <- prefilter_TE_aggregation(gene_TE_counts = gene_TE_counts_raw, 
                                                         aggregate_all = FALSE, 
                                                         intergenic_nearby_repeats = repeats.intergenic.nearby, 
                                                         intergenic_distal_repeats = repeats.intergenic.distal, 
                                                         intronic_repeats = repeats.intronic, 
                                                         exonic_repeats = repeats.exonic)

# Calculate the fraction of samples in each age group (use this info to set the minimum expression threshold)
fraction_group_1 <- length(which(my.meta_data$Age_Group2 == 'A_Young'))/nrow(my.meta_data)
fraction_group_2 <- length(which(my.meta_data$Age_Group2 == 'B_Middle_Aged'))/nrow(my.meta_data)
fraction_group_3 <- length(which(my.meta_data$Age_Group2 == 'C_Old'))/nrow(my.meta_data)
fraction_group_4 <- length(which(my.meta_data$Age_Group2 == 'D_Geriatric'))/nrow(my.meta_data)

# Run function to remove gene/transcript version info, combine transcripts from the same gene, and filter lowly expressed genes
counts_filtered <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw.aggregate, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 0.20)

    # Save the final counts file
    write.table(counts_filtered, file = paste(counts.dir, "TElocal_GSE113957_Aging_Fibroblasts_filtered_counts", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)

# Convert counts to matrix format
counts_filtered <- as.matrix(counts_filtered)

# Remove unneeded files
rm(gene_TE_counts_raw, gene_TE_counts_raw.aggregate)




   
# COVARIATE AND BATCH CORRECTIONS ------------------------------------------------------------------------------------------------------------------------
    
    
# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = rownames(my.meta_data),
                         Age = as.numeric(my.meta_data$Age),
                         Sex = as.factor(my.meta_data$sex),
                         Instrument = as.factor(my.meta_data$Instrument)
                         )

# Set the full models (adjustment + variables of interest).
mod_full = model.matrix(~ Age + Sex + Instrument, data = SampleInfo)

# Remove batch effects and covariates with Limma (everything except age )
my.corrected.data <- removeBatchEffect(x = log2(counts_filtered + 1.0),
                                       batch = NULL,
                                       covariates = mod_full[, c(3:4)],
                                       design = mod_full[, c(1:2)])

# Delog and round data.
gene_TE_counts <- ceiling(2^my.corrected.data-1.0)

    # Save batch corrected counts (to visualize downstream results)
    write.table(gene_TE_counts, file = paste(counts.dir, 'TElocal_GSE113957_Aging_Fibroblasts_filtered_counts_BATCH-CORRECTED', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
  


    
  
# DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS (WALT TEST/LINEAR MODEL ACROSS AGE) ------------------------------------------------------------------------------------------------------------------------


# GENERAL PARAMETERS

# Define alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
padj_limit <- 0.05 
    
# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = rownames(my.meta_data),
                         Age = as.numeric(my.meta_data$Age))

# Create DESeq2 object
my.dds <- DESeqDataSetFromMatrix(countData = gene_TE_counts,
                                  colData = SampleInfo,
                                  design = ~ Age) 

# run DESeq2
my.dds <- DESeq(my.dds, parallel = TRUE) 

# Get VST data (includes size factor normalization)
normalized_counts <- assay(varianceStabilizingTransformation(my.dds, blind = FALSE))

    # Save VST data
    write.table(normalized_counts, file = paste(counts.dir, 'TElocal_GSE113957_Aging_Fibroblasts_filtered_counts_VST', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
 
# Extract results
my_results <- results(my.dds, name = c("Age"), alpha = padj_limit, independentFiltering = TRUE)

# DESeq Stats
summary(my_results, alpha = padj_limit)

# Extract significant results and save.
my_sig <- Extract_DESeq_stats(DESeq_results = my_results, padj_limit = padj_limit, organism = 'hs', 
                                          output.dir = DESeq.dir,
                                          output_file_prefix_all = 'TElocal_GSE113957_Linear_Across_Age_All_Genes',
                                          output_file_prefix_sig = 'TElocal_GSE113957_Linear_Across_Age_FDR5')



    
    
# PLOT ALUJB EXPRESSION SCATTERPLOTS ---------------------------------------------------------------------------------------------------------------------------------


# Collect covariates in a dataframe
my.Alu.exp <- data.frame(row.names = rownames(my.meta_data),
                         Age = my.meta_data$Age,
                         AluJb_exonic = normalized_counts['exonic:AluJb:Alu:SINE', rownames(my.meta_data)],
                         AluJb_inter_near = normalized_counts['intergenic_near:AluJb:Alu:SINE', rownames(my.meta_data)]
                         )

# Extract FDR
exonic.FDR <- my_sig['exonic:AluJb:Alu:SINE', 'padj']
intergenic_near.FDR <- my_sig['intergenic_near:AluJb:Alu:SINE', 'padj']

# Find median expression in young group

    # Define young samples
    my.young <- rownames(my.meta_data[which(my.meta_data$Age_Group2 == 'A_Young'), ])
    
    # Calculate media expression
    young.median.exon.AluJb <- median(my.Alu.exp[my.young, 'AluJb_exonic'])
    young.median.nearby.AluJb <- median(my.Alu.exp[my.young, 'AluJb_inter_near'])

# Generate boxplots (exon-overlapping AluJb)
pdf(paste(DESeq.dir, "GSE113957_Scatterplot_Expression_AluJb_Exon-overlapping", ".pdf", sep=""), width = 6, height = 6)

        # make boxplots 
        plot(AluJb_exonic ~ Age, 
             data = my.Alu.exp,
             pch = 16,
             ylim = c(10.2, 11.2),
             main = paste('Exon-overlapping AluJb Expression', sep = ''),
             ylab = paste('Normalized log2(counts)', sep = ''),
             xlab = c('Age'))
        
        # Draw dashed line at median AluJb expression
        abline(h = young.median.exon.AluJb,lty = 3)
        
        # Add FDR
        text(60, 11.2, paste('FDR = ', signif(exonic.FDR, 3), sep = ''), cex = 1)

# End pdf
dev.off()

# Generate boxplots (intergenic nearby AluJb)
pdf(paste(DESeq.dir, "GSE113957_Scatterplot_Expression_AluJb_intergenic-nearby", ".pdf", sep=""), width = 6, height = 6)

        # make boxplots 
        plot(AluJb_inter_near ~ Age, 
             data = my.Alu.exp,
             pch = 16,
             ylim = c(8.25, 10.25),
             main = paste('Nearby Intergenic AluJb Expression', sep = ''),
             ylab = paste('Normalized log2(counts)', sep = ''),
             xlab = c('Age'))
        
        # Draw dashed line at median AluJb expression
        abline(h = young.median.nearby.AluJb,lty = 3)
        
        # Add FDR
        text(60, 10.25, paste('FDR = ', signif(intergenic_near.FDR, 3), sep = ''), cex = 1)

# End pdf
dev.off()



    

# MDS/PCA PLOTS ------------------------------------------------------------------------------------------------------------------------
      
# Define sample groups   
SampleInfo <- data.frame(row.names = rownames(my.meta_data),
                         Age = my.meta_data$Age,
                         Age_Group = my.meta_data$Age_Group2
                         )

# Generate df to hold plot point colors and shapes
point.parameters <- SampleInfo[, c('Age_Group'), drop = FALSE]
point.parameters$Age_Group <- as.character(point.parameters$Age_Group)
    
# Assign a color to each age group
point.parameters[which(point.parameters$Age_Group == 'A_Young'), "Age_Group"] <- '#FFB000'
point.parameters[which(point.parameters$Age_Group == 'B_Middle_Aged'), "Age_Group"] <- '#FE6100'
point.parameters[which(point.parameters$Age_Group == 'C_Old'), "Age_Group"] <- '#DC267F'
point.parameters[which(point.parameters$Age_Group == 'D_Geriatric'), "Age_Group"] <- '#785EF0'


# Define point shapes and colors
my.shapes <- 16
my.colors <- point.parameters$Age_Group

# Run MDS
Run_MDS(VST_expression = normalized_counts,
        treatment_label = 'GSE113957',
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = paste(MDS.dir, 'TElocal_', sep = ''))

# Run PCA
Run_PCA(VST_expression = normalized_counts,
        treatment_label = 'GSE113957',
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = paste(MDS.dir, 'TElocal_', sep = ''))
    
    

    

# HEATMAPS ---------------------------------------------------------------------------------------------------------------------------------

  
# Load significant DEGs
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/DESeq_Results/TElocal_GSE113957_Linear_Across_Age_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant gene/TE names
my.gene.names <- rownames(my.DEGs[grepl('ENSG', rownames(my.DEGs)),])
my.TE.names <- rownames(my.DEGs[!grepl('ENSG', rownames(my.DEGs)),])

# Extract DESeq results for significant genes/TEs
DESeq.genes <- my.DEGs[my.gene.names, ]
DESeq.TEs <- my.DEGs[my.TE.names, ]

    # Count the number of up/down genes/TEs in each analysis
    Number.genes.up <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange > 0), ])
    Number.genes.down <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange < 0), ])
    
    Number.TEs.up <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange > 0), ])
    Number.TEs.down <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange < 0), ])

# Extract expression values
expression.genes <- normalized_counts[my.gene.names, ]
expression.TEs <- normalized_counts[my.TE.names, ]

# Add column to metadata to hold heatmap labels
my.meta_data$Plot_Group <- my.meta_data$Age_Group1

# Update the column as a factor
my.meta_data$Plot_Group <- as.factor(my.meta_data$Plot_Group)
  
# Define heatmap column group labels and group colors
annot_column <- data.frame(row.names = rownames(my.meta_data), 
                           Plot_Group = my.meta_data$Plot_Group)

annot_color <- list(Plot_Group = c(A_20s = '#F0EAD6',
                                   B_30s = '#CEC9B7',
                                   C_40s = '#ABA799',
                                   D_50s = '#89867A',
                                   E_60s = '#67645C',
                                   F_70s = '#45433D',
                                   G_80s = '#22211F',
                                   H_90s = '#000000')
                    )

# Define row label groups and colors (FOR REPEATS)

    # Make df to hold TE family info
    TE_family_info <- data.frame(row.names = rownames(expression.TEs), 
                                 TEs = rownames(expression.TEs))
    
    # Split the TE label
    TE_family_info <- as.data.frame(splitstackshape::cSplit(TE_family_info, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = FALSE, stripWhite = FALSE, makeEqual = TRUE)) # warning about specifying 'as.is' by the caller can be ignored
    rownames(TE_family_info) <- TE_family_info$TEs
    colnames(TE_family_info) <- c('TEs', 'Location', 'Subfamily', 'Family', 'Class')
  
    # Extract the info for the row labels
    annot_row <- TE_family_info[, c('Class', 'Location'), drop = FALSE]
    
    # Consolidate class names for TEs with unclear or undefined class designations
    annot_row[which(annot_row$Class == 'Retroposon'), 'Class'] <- 'Other'
    annot_row[which(annot_row$Class == 'Satellite'), 'Class'] <- 'Other'
    annot_row[which(annot_row$Class == 'RC'), 'Class'] <- 'Other'
    annot_row[which(annot_row$Class == 'RNA'), 'Class'] <- 'Other'
    
    annot_row[which(annot_row$Class == 'LTR?'), 'Class'] <- 'Other'
    annot_row[which(annot_row$Class == 'RC?'), 'Class'] <- 'Other'
    annot_row[which(annot_row$Class == 'DNA?'), 'Class'] <- 'Other'
    
    annot_row[which(annot_row$Class == 'Unknown'), 'Class'] <- 'Other'
    annot_row[which(is.na(annot_row$Class)), 'Class'] <- 'Other'
        
    # Make the final annotation df, ordering by location then by repeat class 
    annot_row <- annot_row[order(annot_row$Location, annot_row$Class), , drop = FALSE]
  
# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)

# Save heatmap (TEs)
pdf(paste(Heatmaps.dir, "TElocal_GSE113957_Linear_Age_FDR5_REPEATS_Heatmap", ".pdf", sep=""), height = 8, width = 10)

    pheatmap(expression.TEs, 
             cluster_rows = FALSE, 
             show_rownames = FALSE,
             fontsize_row = 4,
             cluster_cols = FALSE, 
             show_colnames = FALSE,
             scale = 'row',
             border_color = NA,
             annotation_col = annot_column,
             annotation_colors = annot_color,
             annotation_row = annot_row,
             color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
             breaks = breaksList,
             main = paste('GSE113957 Aging Human Dermal Fibroblasts \n FDR 5% Repeats Across Age: ', Number.TEs.up, ' Up, ', Number.TEs.down, ' Down', sep = ''))
    
dev.off()





# Session Info ------------------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_DESeq2_MDS_PCA_Heatmaps_TElocal.txt", sep =""))
sessionInfo()
sink()      
    


    
# Clean the environment
rm(list=ls())


