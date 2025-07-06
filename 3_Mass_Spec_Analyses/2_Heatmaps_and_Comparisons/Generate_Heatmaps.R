# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(pheatmap) # For heatmaps
library(gridExtra) # to plot pheatmaps side by side

# Define the output directories
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/2_Heatmaps_and_Comparisons/Plots/'


# CELLULAR PROTEOMICS ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Load protein abundance data (log10 transformed)
my.expression <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Cellular_Proteome_Abundances_log10.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Only keep quantification
    my.expression <- my.expression[, -c(1:2)]

# Load differential abundance results
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Cellular_Proteome_Differential_Abundance_Results.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant protein names
my.gene.names <- rownames(my.DEGs[which(my.DEGs$Qvalue < 0.05 & abs(my.DEGs$AVG.Log2.Ratio) > 0.20), ])

# Extract expression values
expression.genes <- my.expression[my.gene.names, ]

# Extract differential abundance results for significant proteins
diff.abundance.proteins <- my.DEGs[my.gene.names, ]

# Count the number of up/down proteins 
Number.genes.up <- nrow(diff.abundance.proteins[which(diff.abundance.proteins$AVG.Log2.Ratio > 0), ])
Number.genes.down <- nrow(diff.abundance.proteins[which(diff.abundance.proteins$AVG.Log2.Ratio < 0), ])
  
# Define column label groups and group colors
annot_column <- data.frame(row.names = colnames(my.expression), 
                           Treatment = c(rep('Empty', 5), rep('AluJb_OE', 5)))

annot_color <- list(Treatment = c(Empty = 'dodgerblue', AluJb_OE = 'dodgerblue4'))

# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)
    
# Generate heatmaps

    # Genes 
    p1 <- pheatmap(expression.genes, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   fontsize_row = 4,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('Cellular Proteomics \n Q < 5% and |log2 FC| > 0.20 Genes: ', Number.genes.up, ' Up, ', Number.genes.down, ' Down', sep = ''))

# Save heatmaps (in one file, side by side)
pdf(paste(dir.output, "Plot_Heatmap_Significant_Differential_Abundance_Cellular_Proteomics", ".pdf", sep=""), height = 8, width = 10)

  # Modify plot list to plot heatmaps side by side
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  #plot_list[['p2']]=p2[[4]]
    
  # Arrange the plot grid
  grid.arrange(grobs=plot_list, ncol=2)
    
dev.off()
dev.off()


  
# SECRETED PROTEOMICS ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Load protein abundance data (log10 transformed)
my.expression <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Secreted_Proteome_Abundances_log10.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Only keep quantification
    my.expression <- my.expression[, -c(1:7)]

# Load differential abundance results
my.DEGs <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Secreted_Proteome_Differential_Abundance_Results.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Extract significant protein names
my.gene.names <- rownames(my.DEGs[which(my.DEGs$Qvalue < 0.05 & abs(my.DEGs$AVG.Log2.Ratio) > 0.58), ])

# Extract expression values
expression.genes <- my.expression[my.gene.names, ]

# Extract differential abundance results for significant proteins
diff.abundance.proteins <- my.DEGs[my.gene.names, ]

# Count the number of up/down proteins 
Number.genes.up <- nrow(diff.abundance.proteins[which(diff.abundance.proteins$AVG.Log2.Ratio > 0), ])
Number.genes.down <- nrow(diff.abundance.proteins[which(diff.abundance.proteins$AVG.Log2.Ratio < 0), ])
  
# Define column label groups and group colors
annot_column <- data.frame(row.names = colnames(my.expression), 
                           Treatment = c(rep('Empty', 4), rep('AluJb_OE', 4)))

annot_color <- list(Treatment = c(Empty = 'dodgerblue', AluJb_OE = 'dodgerblue4'))

# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)
    
# Generate heatmaps

    # Genes 
    p1 <- pheatmap(expression.genes, 
                   cluster_rows = TRUE, 
                   show_rownames = FALSE,
                   fontsize_row = 4,
                   cluster_cols = FALSE, 
                   show_colnames = FALSE,
                   scale = 'row',
                   border_color = NA,
                   annotation_col = annot_column,
                   annotation_colors = annot_color,
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList,
                   main = paste('Secreted Proteomics \n Q < 5% and |log2 FC| > 0.58 Genes: ', Number.genes.up, ' Up, ', Number.genes.down, ' Down', sep = ''))

# Save heatmaps (in one file, side by side)
pdf(paste(dir.output, "Plot_Heatmap_Significant_Differential_Abundance_Secreted_Proteomics", ".pdf", sep=""), height = 8, width = 10)

  # Modify plot list to plot heatmaps side by side
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  #plot_list[['p2']]=p2[[4]]
    
  # Arrange the plot grid
  grid.arrange(grobs=plot_list, ncol=2)
    
dev.off()
dev.off()  



# SESSION INFO ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/2_Heatmaps_and_Comparisons/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_Heatmaps.txt", sep =""))
sessionInfo()
sink()      
    

# Clean the environment
rm(list=ls())