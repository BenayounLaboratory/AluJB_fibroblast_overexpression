# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(pheatmap) # For heatmaps
library(gridExtra) # to plot pheatmaps side by side

# Define the output directories
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/4_Heatmaps_and_Comparisons/Plots/'



# Section 1: T24 AluJb ---------------------------------------------------------------------------------------------------------------------------------

  
# Load VST data
my.expression <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/VST_Expression_T24_AluJb.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Load significant DEGs
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T24_AluJb_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant gene/TE names
my.gene.names <- rownames(my.DEGs[grepl('ENSG', rownames(my.DEGs)),])
my.TE.names <- rownames(my.DEGs[!grepl('ENSG', rownames(my.DEGs)),])

# Extract expression values
expression.genes <- my.expression[my.gene.names, ]
expression.TEs <- my.expression[my.TE.names, ]
    
# Extract DESeq results for significant genes/TEs
DESeq.genes <- my.DEGs[my.gene.names, ]
DESeq.TEs <- my.DEGs[my.TE.names, ]

# Count the number of up/down genes/TEs in each analysis
Number.genes.up <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange > 0), ])
Number.genes.down <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange < 0), ])

Number.TEs.up <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange > 0), ])
Number.TEs.down <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange < 0), ])
  
# Define column label groups and group colors
annot_column <- data.frame(row.names = colnames(my.expression), 
                           Treatment = c(rep('Empty', 4), rep('AluJb_OE', 4)))

annot_color <- list(Treatment = c(Empty = 'dodgerblue', AluJb_OE = 'dodgerblue4'))

# Define row label groups and colors

    # Make df to hold TE family info
    TE_family_info <- data.frame(row.names = rownames(expression.TEs), 
                                 TEs = rownames(expression.TEs))
    
    # Split the TE label
    # The warning message that " 'as.is' should be specified by the caller" does not affect results and can be ignored
    TE_family_info <- as.data.frame(splitstackshape::cSplit(TE_family_info, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = FALSE, stripWhite = FALSE, makeEqual = TRUE)) # warning about specifying 'as.is' by the caller can be ignored
    rownames(TE_family_info) <- TE_family_info$TEs
    colnames(TE_family_info) <- c('TEs', 'Subfamily', 'Family', 'Class')
  
    # Make the final annotation df, ordering by subfamily name
    annot_row <- TE_family_info[, 'Class', drop = FALSE]
    annot_row <- annot_row[order(annot_row$Class), , drop = FALSE]
    
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
  
# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)
    
# Generate heatmaps

    # Genes 
    p1 <- pheatmap(expression.genes, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   fontsize_row = 1,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T24 RNA-Seq \n FDR 5% Genes: ', Number.genes.up, ' Up, ', Number.genes.down, ' Down', sep = ''))

    # TEs 
    p2 <- pheatmap(expression.TEs[rownames(annot_row), ], 
                   cluster_rows = FALSE, 
                   show_rownames = FALSE,
                   fontsize_row = 2.0,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   annotation_row = annot_row,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T24 RNA-Seq \n FDR 5% TEs: ', Number.TEs.up, ' Up, ', Number.TEs.down, ' Down', sep = ''))

    
# Save heatmaps (in one file, side by side)
pdf(paste(dir.output, "Plot_Heatmap_DEGs_FDR5_T24_Alu_IMR90", ".pdf", sep=""), height = 8, width = 10)

  # Modify plot list to plot heatmaps side by side
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
    
  # Arrange the plot grid
  grid.arrange(grobs=plot_list, ncol=2)
    
dev.off()
dev.off()


  
  

# Section 2: Proteomics ---------------------------------------------------------------------------------------------------------------------------------


# Load VST data
my.expression <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/VST_Expression_T72_Proteomics.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Load significant DEGs
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant gene/TE names
my.gene.names <- rownames(my.DEGs[grepl('ENSG', rownames(my.DEGs)),])
my.TE.names <- rownames(my.DEGs[!grepl('ENSG', rownames(my.DEGs)),])

# Extract expression values
expression.genes <- my.expression[my.gene.names, ]
expression.TEs <- my.expression[my.TE.names, ]
    
# Extract DESeq results for significant genes/TEs
DESeq.genes <- my.DEGs[my.gene.names, ]
DESeq.TEs <- my.DEGs[my.TE.names, ]

# Count the number of up/down genes/TEs in each analysis
Number.genes.up <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange > 0), ])
Number.genes.down <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange < 0), ])

Number.TEs.up <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange > 0), ])
Number.TEs.down <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange < 0), ])
  
# Define column label groups and group colors
annot_column <- data.frame(row.names = colnames(my.expression), 
                           Treatment = c(rep('Empty', 4), rep('AluJb_OE', 4)))

annot_color <- list(Treatment = c(Empty = 'dodgerblue', AluJb_OE = 'dodgerblue4'))

# Define row label groups and colors

    # Make df to hold TE family info
    TE_family_info <- data.frame(row.names = rownames(expression.TEs), 
                                 TEs = rownames(expression.TEs))
    
    # Split the TE label
    # The warning message that " 'as.is' should be specified by the caller" does not affect results and can be ignored
    TE_family_info <- as.data.frame(splitstackshape::cSplit(TE_family_info, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = FALSE, stripWhite = FALSE, makeEqual = TRUE)) # warning about specifying 'as.is' by the caller can be ignored
    rownames(TE_family_info) <- TE_family_info$TEs
    colnames(TE_family_info) <- c('TEs', 'Subfamily', 'Family', 'Class')
  
    # Make the final annotation df, ordering by subfamily name
    annot_row <- TE_family_info[, 'Class', drop = FALSE]
    annot_row <- annot_row[order(annot_row$Class), , drop = FALSE]
    
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
  
# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)
    
# Generate heatmaps

    # Genes 
    p1 <- pheatmap(expression.genes, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   fontsize_row = 1,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE Proteomics RNA-Seq \n FDR 5% Genes: ', Number.genes.up, ' Up, ', Number.genes.down, ' Down', sep = ''))

    # TEs 
    p2 <- pheatmap(expression.TEs[rownames(annot_row), ], 
                   cluster_rows = FALSE, 
                   show_rownames = FALSE,
                   fontsize_row = 2.0,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   annotation_row = annot_row,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE Proteomics RNA-Seq \n FDR 5% TEs: ', Number.TEs.up, ' Up, ', Number.TEs.down, ' Down', sep = ''))

    
# Save heatmaps (in one file, side by side)
pdf(paste(dir.output, "Plot_Heatmap_DEGs_FDR5_Proteomics_Alu_IMR90", ".pdf", sep=""), height = 8, width = 10)

  # Modify plot list to plot heatmaps side by side
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
    
  # Arrange the plot grid
  grid.arrange(grobs=plot_list, ncol=2)
    
dev.off()
dev.off()


  
  


# Section 3: Lamivudine ---------------------------------------------------------------------------------------------------------------------------------


# DEFINE GENERAL PARAMETERS 

# Load VST data
my.expression <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/VST_Expression_T72_Lamivudine.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Define sample indices in the expression matrix
index.AluEff.Veh <- c(1:4, 8:11)
index.AluEff.3TC <- c(5:7, 12:15)
index.LamEff.Empty <- c(1:4, 5:7)
index.LamEff.Alu <- c(8:11, 12:15)

# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)





# ALU EFFECT, VEHICLE BACKGROUND

# Load significant DEGs
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_LamExperiment_Alu_vs_Empty_Veh_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant gene/TE names
my.gene.names <- rownames(my.DEGs[grepl('ENSG', rownames(my.DEGs)),])
my.TE.names <- rownames(my.DEGs[!grepl('ENSG', rownames(my.DEGs)),])

# Extract expression values
expression.genes <- my.expression[my.gene.names, index.AluEff.Veh]
expression.TEs <- my.expression[my.TE.names, index.AluEff.Veh]
    
# Extract DESeq results for significant genes/TEs
DESeq.genes <- my.DEGs[my.gene.names, ]
DESeq.TEs <- my.DEGs[my.TE.names, ]

# Count the number of up/down genes/TEs in each analysis
Number.genes.up <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange > 0), ])
Number.genes.down <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange < 0), ])

Number.TEs.up <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange > 0), ])
Number.TEs.down <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange < 0), ])
  
# Define column label groups and group colors
annot_column <- data.frame(row.names = colnames(expression.genes), 
                           Treatment = c(rep('Empty', 4), rep('AluJb_OE', 4)))

annot_color <- list(Treatment = c(Empty = 'dodgerblue', AluJb_OE = 'dodgerblue4'))

# Define row label groups and colors

    # Make df to hold TE family info
    TE_family_info <- data.frame(row.names = rownames(expression.TEs),
                                 TEs = rownames(expression.TEs))

    # Split the TE label
    # The warning message that " 'as.is' should be specified by the caller" does not affect results and can be ignored
    TE_family_info <- as.data.frame(splitstackshape::cSplit(TE_family_info, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = FALSE, stripWhite = FALSE, makeEqual = TRUE)) # warning about specifying 'as.is' by the caller can be ignored
    rownames(TE_family_info) <- TE_family_info$TEs
    colnames(TE_family_info) <- c('TEs', 'Subfamily', 'Family', 'Class')

    # Make the final annotation df, ordering by subfamily name
    annot_row <- TE_family_info[, 'Class', drop = FALSE]
    annot_row <- annot_row[order(annot_row$Class), , drop = FALSE]

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
    
# Generate heatmaps

    # Genes 
    p1 <- pheatmap(expression.genes, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   fontsize_row = 1,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T72 Alu Effect (Vehicle Background) \n FDR 5% Genes: ', Number.genes.up, ' Up, ', Number.genes.down, ' Down', sep = ''))

    # TEs
    p2 <- pheatmap(expression.TEs[rownames(annot_row), ],
                   cluster_rows = TRUE,
                   show_rownames = FALSE,
                   fontsize_row = 2.0,
                   cluster_cols = FALSE,
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   annotation_row = annot_row,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T72 Alu Effect (Vehicle Background) \n FDR 5% TEs: ', Number.TEs.up, ' Up, ', Number.TEs.down, ' Down', sep = ''))

    
# Save heatmaps (in one file, side by side)
pdf(paste(dir.output, "Plot_Heatmap_DEGs_FDR5_LamExperiment_IMR90_Alu_Effect_Vehicle", ".pdf", sep=""), height = 8, width = 10)

  # Modify plot list to plot heatmaps side by side
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
    
  # Arrange the plot grid
  grid.arrange(grobs=plot_list, ncol=2)
    
dev.off()
dev.off()





# ALU EFFECT, 3TC BACKGROUND

# Load significant DEGs
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_LamExperiment_Alu_vs_Empty_3TC_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant gene/TE names
my.gene.names <- rownames(my.DEGs[grepl('ENSG', rownames(my.DEGs)),])
my.TE.names <- rownames(my.DEGs[!grepl('ENSG', rownames(my.DEGs)),])

# Extract expression values
expression.genes <- my.expression[my.gene.names, index.AluEff.3TC]
expression.TEs <- my.expression[my.TE.names, index.AluEff.3TC]
    
# Extract DESeq results for significant genes/TEs
DESeq.genes <- my.DEGs[my.gene.names, ]
DESeq.TEs <- my.DEGs[my.TE.names, ]

# Count the number of up/down genes/TEs in each analysis
Number.genes.up <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange > 0), ])
Number.genes.down <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange < 0), ])

Number.TEs.up <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange > 0), ])
Number.TEs.down <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange < 0), ])
  
# Define column label groups and group colors
annot_column <- data.frame(row.names = colnames(expression.genes), 
                           Treatment = c(rep('Empty', 3), rep('AluJb_OE', 4)))

annot_color <- list(Treatment = c(Empty = 'dodgerblue', AluJb_OE = 'dodgerblue4'))

# # Define row label groups and colors
# 
#     # Make df to hold TE family info
#     TE_family_info <- data.frame(row.names = rownames(expression.TEs), 
#                                  TEs = rownames(expression.TEs))
#     
#     # Split the TE label
#     TE_family_info <- as.data.frame(splitstackshape::cSplit(TE_family_info, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = FALSE, stripWhite = FALSE, makeEqual = TRUE)) # warning about specifying 'as.is' by the caller can be ignored
#     rownames(TE_family_info) <- TE_family_info$TEs
#     colnames(TE_family_info) <- c('TEs', 'Subfamily', 'Family', 'Class')
#   
#     # Make the final annotation df, ordering by subfamily name
#     annot_row <- TE_family_info[, 'Class', drop = FALSE]
#     annot_row <- annot_row[order(annot_row$Class), , drop = FALSE]
#     
#     # Consolidate class names for TEs with unclear or undefined class designations
#     annot_row[which(annot_row$Class == 'Retroposon'), 'Class'] <- 'Other'
#     annot_row[which(annot_row$Class == 'Satellite'), 'Class'] <- 'Other'
#     annot_row[which(annot_row$Class == 'RC'), 'Class'] <- 'Other'
#     annot_row[which(annot_row$Class == 'RNA'), 'Class'] <- 'Other'
#     
#     annot_row[which(annot_row$Class == 'LTR?'), 'Class'] <- 'Other'
#     annot_row[which(annot_row$Class == 'RC?'), 'Class'] <- 'Other'
#     annot_row[which(annot_row$Class == 'DNA?'), 'Class'] <- 'Other'
#     
#     annot_row[which(annot_row$Class == 'Unknown'), 'Class'] <- 'Other'
#     annot_row[which(is.na(annot_row$Class)), 'Class'] <- 'Other'
    
# Generate heatmaps

    # Genes 
    p1 <- pheatmap(expression.genes, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   fontsize_row = 1,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T72 Alu Effect (3TC Background) \n FDR 5% Genes: ', Number.genes.up, ' Up, ', Number.genes.down, ' Down', sep = ''))

    # # TEs 
    # p2 <- pheatmap(expression.TEs[rownames(annot_row), ], 
    #                cluster_rows = FALSE, 
    #                show_rownames = FALSE,
    #                fontsize_row = 2.0,
    #                cluster_cols = FALSE, 
    #                show_colnames = FALSE,
    #                scale = 'row',
    #                border_color = NA,
    #                annotation_col = annot_column,
    #                annotation_colors = annot_color,
    #                annotation_row = annot_row,
    #                color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
    #                breaks = breaksList,
    #                main = paste('AluJb OE T72 Alu Effect (3TC Background) \n FDR 5% TEs: ', Number.TEs.up, ' Up, ', Number.TEs.down, ' Down', sep = ''))

    
# Save heatmaps (in one file, side by side)
pdf(paste(dir.output, "Plot_Heatmap_DEGs_FDR5_LamExperiment_IMR90_Alu_Effect_3TC", ".pdf", sep=""), height = 8, width = 10)

  # Modify plot list to plot heatmaps side by side
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  #plot_list[['p2']]=p2[[4]]
    
  # Arrange the plot grid
  grid.arrange(grobs=plot_list, ncol=2)
    
dev.off()
dev.off()





# 3TC EFFECT, EMPTY PLASMID BACKGROUND

# Load significant DEGs
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_LamExperiment_3TC_vs_Veh_Empty_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant gene/TE names
my.gene.names <- rownames(my.DEGs[grepl('ENSG', rownames(my.DEGs)),])
my.TE.names <- rownames(my.DEGs[!grepl('ENSG', rownames(my.DEGs)),])

# Extract expression values
expression.genes <- my.expression[my.gene.names, index.LamEff.Empty]
expression.TEs <- my.expression[my.TE.names, index.LamEff.Empty]
    
# Extract DESeq results for significant genes/TEs
DESeq.genes <- my.DEGs[my.gene.names, ]
DESeq.TEs <- my.DEGs[my.TE.names, ]

# Count the number of up/down genes/TEs in each analysis
Number.genes.up <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange > 0), ])
Number.genes.down <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange < 0), ])

Number.TEs.up <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange > 0), ])
Number.TEs.down <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange < 0), ])
  
# Define column label groups and group colors
annot_column <- data.frame(row.names = colnames(expression.genes), 
                           Treatment = c(rep('Vehicle', 4), rep('Lamivudine', 3)))

annot_color <- list(Treatment = c(Vehicle = 'mediumpurple', Lamivudine = 'mediumpurple4'))

# Define row label groups and colors

    # Make df to hold TE family info
    TE_family_info <- data.frame(row.names = rownames(expression.TEs),
                                 TEs = rownames(expression.TEs))

    # Split the TE label
    TE_family_info <- as.data.frame(splitstackshape::cSplit(TE_family_info, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = FALSE, stripWhite = FALSE, makeEqual = TRUE)) # warning about specifying 'as.is' by the caller can be ignored
    rownames(TE_family_info) <- TE_family_info$TEs
    colnames(TE_family_info) <- c('TEs', 'Subfamily', 'Family', 'Class')

    # Make the final annotation df, ordering by subfamily name
    annot_row <- TE_family_info[, 'Class', drop = FALSE]
    annot_row <- annot_row[order(annot_row$Class), , drop = FALSE]

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
    
# Generate heatmaps

    # Genes 
    p1 <- pheatmap(expression.genes, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   fontsize_row = 1,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T72 3TC Effect (Empty Background) \n FDR 5% Genes: ', Number.genes.up, ' Up, ', Number.genes.down, ' Down', sep = ''))

    # TEs
    p2 <- pheatmap(expression.TEs[rownames(annot_row), ],
                   cluster_rows = FALSE,
                   show_rownames = FALSE,
                   fontsize_row = 2.0,
                   cluster_cols = FALSE,
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   annotation_row = annot_row,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T72 3TC Effect (Empty Background) \n FDR 5% TEs: ', Number.TEs.up, ' Up, ', Number.TEs.down, ' Down', sep = ''))

    
# Save heatmaps (in one file, side by side)
pdf(paste(dir.output, "Plot_Heatmap_DEGs_FDR5_LamExperiment_IMR90_3TC_Effect_Empty", ".pdf", sep=""), height = 8, width = 10)

  # Modify plot list to plot heatmaps side by side
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
    
  # Arrange the plot grid
  grid.arrange(grobs=plot_list, ncol=2)
    
dev.off()
dev.off()






# 3TC EFFECT, ALU PLASMID BACKGROUND

# Load significant DEGs
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_LamExperiment_3TC_vs_Veh_Alu_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant gene/TE names
my.gene.names <- rownames(my.DEGs[grepl('ENSG', rownames(my.DEGs)),])
my.TE.names <- rownames(my.DEGs[!grepl('ENSG', rownames(my.DEGs)),])

# Extract expression values
expression.genes <- my.expression[my.gene.names, index.LamEff.Alu]
expression.TEs <- my.expression[my.TE.names, index.LamEff.Alu]
    
# Extract DESeq results for significant genes/TEs
DESeq.genes <- my.DEGs[my.gene.names, ]
DESeq.TEs <- my.DEGs[my.TE.names, ]

# Count the number of up/down genes/TEs in each analysis
Number.genes.up <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange > 0), ])
Number.genes.down <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange < 0), ])

Number.TEs.up <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange > 0), ])
Number.TEs.down <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange < 0), ])
  
# Define column label groups and group colors
annot_column <- data.frame(row.names = colnames(expression.genes), 
                           Treatment = c(rep('Vehicle', 4), rep('Lamivudine', 4)))

annot_color <- list(Treatment = c(Vehicle = 'mediumpurple', Lamivudine = 'mediumpurple4'))

# Define row label groups and colors

    # Make df to hold TE family info
    TE_family_info <- data.frame(row.names = rownames(expression.TEs),
                                 TEs = rownames(expression.TEs))

    # Split the TE label
    TE_family_info <- as.data.frame(splitstackshape::cSplit(TE_family_info, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = FALSE, stripWhite = FALSE, makeEqual = TRUE)) # warning about specifying 'as.is' by the caller can be ignored
    rownames(TE_family_info) <- TE_family_info$TEs
    colnames(TE_family_info) <- c('TEs', 'Subfamily', 'Family', 'Class')

    # Make the final annotation df, ordering by subfamily name
    annot_row <- TE_family_info[, 'Class', drop = FALSE]
    annot_row <- annot_row[order(annot_row$Class), , drop = FALSE]

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
    
# Generate heatmaps

    # Genes 
    p1 <- pheatmap(expression.genes, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   fontsize_row = 1,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T72 3TC Effect (Alu Background) \n FDR 5% Genes: ', Number.genes.up, ' Up, ', Number.genes.down, ' Down', sep = ''))

    # TEs
    p2 <- pheatmap(expression.TEs[rownames(annot_row), ],
                   cluster_rows = TRUE,
                   show_rownames = FALSE,
                   fontsize_row = 2.0,
                   cluster_cols = FALSE,
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   annotation_row = annot_row,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T72 3TC Effect (Alu Background) \n FDR 5% TEs: ', Number.TEs.up, ' Up, ', Number.TEs.down, ' Down', sep = ''))

    
# Save heatmaps (in one file, side by side)
pdf(paste(dir.output, "Plot_Heatmap_DEGs_FDR5_LamExperiment_IMR90_3TC_Effect_Alu", ".pdf", sep=""), height = 8, width = 10)

  # Modify plot list to plot heatmaps side by side
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
    
  # Arrange the plot grid
  grid.arrange(grobs=plot_list, ncol=2)
    
dev.off()
dev.off()  
  




# Section 4: T72 WI38 ---------------------------------------------------------------------------------------------------------------------------------

  


# Load VST data
my.expression <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/VST_Expression_T72_WI38.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Load significant DEGs
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_WI38_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant gene/TE names
my.gene.names <- rownames(my.DEGs[grepl('ENSG', rownames(my.DEGs)),])
my.TE.names <- rownames(my.DEGs[!grepl('ENSG', rownames(my.DEGs)),])

# Extract expression values
expression.genes <- my.expression[my.gene.names, ]
expression.TEs <- my.expression[my.TE.names, ]
    
# Extract DESeq results for significant genes/TEs
DESeq.genes <- my.DEGs[my.gene.names, ]
DESeq.TEs <- my.DEGs[my.TE.names, ]

# Count the number of up/down genes/TEs in each analysis
Number.genes.up <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange > 0), ])
Number.genes.down <- nrow(DESeq.genes[which(DESeq.genes$log2FoldChange < 0), ])

Number.TEs.up <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange > 0), ])
Number.TEs.down <- nrow(DESeq.TEs[which(DESeq.TEs$log2FoldChange < 0), ])
  
# Define column label groups and group colors
annot_column <- data.frame(row.names = colnames(my.expression), 
                           Treatment = c(rep('Empty', 4), rep('AluJb_OE', 4)))

annot_color <- list(Treatment = c(Empty = 'dodgerblue', AluJb_OE = 'dodgerblue4'))

# Define row label groups and colors

    # Make df to hold TE family info
    TE_family_info <- data.frame(row.names = rownames(expression.TEs), 
                                 TEs = rownames(expression.TEs))
    
    # Split the TE label
    # The warning message that " 'as.is' should be specified by the caller" does not affect results and can be ignored
    TE_family_info <- as.data.frame(splitstackshape::cSplit(TE_family_info, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = FALSE, stripWhite = FALSE, makeEqual = TRUE)) # warning about specifying 'as.is' by the caller can be ignored
    rownames(TE_family_info) <- TE_family_info$TEs
    colnames(TE_family_info) <- c('TEs', 'Subfamily', 'Family', 'Class')
  
    # Make the final annotation df, ordering by subfamily name
    annot_row <- TE_family_info[, 'Class', drop = FALSE]
    annot_row <- annot_row[order(annot_row$Class), , drop = FALSE]
    
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
  
# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)
    
# Generate heatmaps

    # Genes 
    p1 <- pheatmap(expression.genes, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   fontsize_row = 1,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T72 WI-38 RNA-Seq \n FDR 5% Genes: ', Number.genes.up, ' Up, ', Number.genes.down, ' Down', sep = ''))

    # TEs 
    p2 <- pheatmap(expression.TEs[rownames(annot_row), ], 
                   cluster_rows = FALSE, 
                   show_rownames = FALSE,
                   fontsize_row = 2.0,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   annotation_row = annot_row,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('AluJb OE T72 WI-38 RNA-Seq \n FDR 5% TEs: ', Number.TEs.up, ' Up, ', Number.TEs.down, ' Down', sep = ''))

    
# Save heatmaps (in one file, side by side)
pdf(paste(dir.output, "Plot_Heatmap_DEGs_FDR5_T72_Alu_WI38", ".pdf", sep=""), height = 8, width = 10)

  # Modify plot list to plot heatmaps side by side
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
    
  # Arrange the plot grid
  grid.arrange(grobs=plot_list, ncol=2)
    
dev.off()
dev.off()


  
  

# Session Info ---------------------------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/4_Heatmaps_and_Comparisons/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_Heatmaps.txt", sep =""))
sessionInfo()
sink()      

    
# Clean the environment
rm(list=ls())

