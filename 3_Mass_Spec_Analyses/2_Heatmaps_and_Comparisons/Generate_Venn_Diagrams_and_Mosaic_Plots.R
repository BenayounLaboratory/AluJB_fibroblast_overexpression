# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(readxl) # To read excel
library(eulerr) # for making venn diagrams

# Define the output directories
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/2_Heatmaps_and_Comparisons/Plots/'


# AluJb OE (Cell Proteins) vs Core SASP  ----------------------------------------------------------------------------------------


# Load Basisty Core SASP table
Basisty.core <- as.data.frame(read_excel("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Aging_and_Senescence_Gene_Lists/Basisty_2020_Core_SASP_CLEANED_UP.xlsx", sheet = "IR_RAS_ATV fibroblast SASP", skip = 2))

# Load differential expression results
DEGs.all <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Cellular_Proteome_Differential_Abundance_Results.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Extract significant protein names
    sig.protein.names <- DEGs.all[which(DEGs.all$Qvalue < 0.05 & abs(DEGs.all$AVG.Log2.Ratio) > 0.20), "UniProtIds"]
    
    # Extract sig protein results
    DEGs.sig <- DEGs.all[sig.protein.names, ]
    
    # make dataframe with non-sig proteins only
    DEGs.nonsig <- DEGs.all[setdiff(rownames(DEGs.all), rownames(DEGs.sig)), ]
    
# Set a seed
set.seed(90280)

# Calculate the size of each set *in the plot
size.overlap <- sum(rownames(DEGs.sig) %in% Basisty.core$Uniprot)
size.SASP <- nrow(Basisty.core) - size.overlap
size.Alu <- nrow(DEGs.sig) - size.overlap
size.universe <- length(union(rownames(DEGs.all), Basisty.core$Uniprot)) # universe = all expressed proteins + proteins in senescence/aging list

# Define the element to plot
my.overlaps <- c('Core_SASP' = size.SASP, 'AluOE' = size.Alu, 'Core_SASP&AluOE' = size.overlap)

# Run stats on the overlap

  # Generate a contingency table 
  my.contingency <- matrix(c(size.universe - sum(size.overlap, size.SASP, size.Alu), size.Alu, size.SASP, size.overlap), nrow=2)
  
  # Run Fischer's Exact Test for enrichment of the overlap
  my.Fischer.res <- fisher.test(my.contingency, alternative = "greater")

  
# Generate the plot
pdf(file = paste(dir.output, "Venn_Diagram_Alu_Cellular_Proteome_vs_CORE_SASP", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(euler(my.overlaps), quantities = TRUE, main = paste('p = ', signif(my.Fischer.res$p.value, 3), sep = ''))

dev.off()


# MAKE MOSAIC PLOT

# calculate missing values
nonsig_in_list <- sum(rownames(DEGs.nonsig) %in% Basisty.core$Uniprot)
nonsig_not_in_list <- nrow(DEGs.nonsig) - nonsig_in_list

# make matrix to hold the data
my.mosaic <- matrix(c(size.overlap, nonsig_in_list, size.Alu, nonsig_not_in_list), 2, 2)
rownames(my.mosaic) <- c('Sig', 'Nonsignificant')
colnames(my.mosaic) <- c('Yes', 'No')

# Run Fischer's Exact Test comparing frequencies
my.Fischer.mosaic <- fisher.test(my.mosaic)

# Generate the plot
pdf(file = paste(dir.output, "Mosaic_plot_Alu_Cellular_Proteome_vs_CORE_SASP", '.pdf', sep=""), width = 5, height = 5)

    mosaicplot(my.mosaic,
               color = TRUE,
               xlab = "Cellular Proteome",
               ylab = "Core SASP",
               main = paste('Fishers exact p-value = ', signif(my.Fischer.mosaic$p.value, 3), sep = '')
               )

dev.off()



# AluJb OE (Secreted Proteins) vs Core SASP  ----------------------------------------------------------------------------------------


# Load Basisty Core SASP table
Basisty.core <- as.data.frame(read_excel("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Aging_and_Senescence_Gene_Lists/Basisty_2020_Core_SASP_CLEANED_UP.xlsx", sheet = "IR_RAS_ATV fibroblast SASP", skip = 2))

# Load differential expression results
DEGs.all <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Secreted_Proteome_Differential_Abundance_Results.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Extract significant protein names
    sig.protein.names <- DEGs.all[which(DEGs.all$Qvalue < 0.05 & abs(DEGs.all$AVG.Log2.Ratio) > 0.58), "UniProtIds"]
    
    # Extract sig protein results
    DEGs.sig <- DEGs.all[sig.protein.names, ]
    
    # make dataframe with non-sig proteins only
    DEGs.nonsig <- DEGs.all[setdiff(rownames(DEGs.all), rownames(DEGs.sig)), ]
    
# Set a seed
set.seed(90280)

# Calculate the size of each set *in the plot
size.overlap <- sum(rownames(DEGs.sig) %in% Basisty.core$Uniprot)
size.SASP <- nrow(Basisty.core) - size.overlap
size.Alu <- nrow(DEGs.sig) - size.overlap
size.universe <- length(union(rownames(DEGs.all), Basisty.core$Uniprot)) # universe = all secreted proteins + proteins in senescence/aging list

# Define the element to plot
my.overlaps <- c('Core_SASP' = size.SASP, 'AluOE' = size.Alu, 'Core_SASP&AluOE' = size.overlap)

# Run stats on the overlap

  # Generate a contingency table 
  my.contingency <- matrix(c(size.universe - sum(size.overlap, size.SASP, size.Alu), size.Alu, size.SASP, size.overlap), nrow=2)
  
  # Run Fischer's Exact Test for enrichment of the overlap
  my.Fischer.res <- fisher.test(my.contingency, alternative = "greater")

  
# Generate the plot
pdf(file = paste(dir.output, "Venn_Diagram_Alu_Secreted_Proteome_vs_CORE_SASP", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(euler(my.overlaps), quantities = TRUE, main = paste('p = ', signif(my.Fischer.res$p.value, 3), sep = ''))

dev.off()


# MAKE MOSAIC PLOT

# calculate missing values
nonsig_in_list <- sum(rownames(DEGs.nonsig) %in% Basisty.core$Uniprot)
nonsig_not_in_list <- nrow(DEGs.nonsig) - nonsig_in_list

# make matrix to hold the data
my.mosaic <- matrix(c(size.overlap, nonsig_in_list, size.Alu, nonsig_not_in_list), 2, 2)
rownames(my.mosaic) <- c('Sig', 'Nonsignificant')
colnames(my.mosaic) <- c('Yes', 'No')

# Run Fischer's Exact Test comparing frequencies
my.Fischer.mosaic <- fisher.test(my.mosaic)

# Generate the plot
pdf(file = paste(dir.output, "Mosaic_plot_Alu_Secreted_Proteome_vs_CORE_SASP", '.pdf', sep=""), width = 5, height = 5)

    mosaicplot(my.mosaic,
               color = TRUE,
               xlab = "Secreted Proteome",
               ylab = "Core SASP",
               main = paste('Fishers exact p-value = ', signif(my.Fischer.mosaic$p.value, 3), sep = '')
               )

dev.off()



# Common proteins across transcriptome, cell proteins, and secreted proteins  ----------------------------------------------------------------------------------------


# Load significant RNA DESeq results
my.transcriptome.sig <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
        
        # Only keep genes (remove TEs)
        my.transcriptome.sig <- my.transcriptome.sig[grepl('ENSG', rownames(my.transcriptome.sig)), ]
    
        # Assign EnsembleID to entries with an NA or Empty cell for HGNC Symbol
        my.transcriptome.sig[is.na(my.transcriptome.sig$hgnc_symbol), 'hgnc_symbol'] <- rownames(my.transcriptome.sig[is.na(my.transcriptome.sig$hgnc_symbol), ])
        my.transcriptome.sig[which(my.transcriptome.sig$hgnc_symbol == ''), 'hgnc_symbol'] <- rownames(my.transcriptome.sig[which(my.transcriptome.sig$hgnc_symbol == ''), ])
        
        # Remove duplicate gene names
        my.transcriptome.sig <- my.transcriptome.sig[!duplicated(my.transcriptome.sig$hgnc_symbol), ]
        
        # Only keep symbol 
        my.transcriptome.sig <- my.transcriptome.sig$hgnc_symbol
        
# Load differential cell protein results
cell.proteins <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Cellular_Proteome_Differential_Abundance_Results.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Split concatenated gene names
    # The warning message that " 'as.is' should be specified by the caller" does not affect results and can be ignored
    cell.proteins <- as.data.frame(splitstackshape::cSplit(cell.proteins, splitCols = "Genes", sep = ";", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE)) # warning about specifying 'as.is' by the caller can be ignored

    # Remove duplicate gene names
    cell.proteins <- cell.proteins[!duplicated(cell.proteins$Genes), ]
    
    # Extract gene names for significant proteins
    cell.proteins.sig <- cell.proteins[which(cell.proteins$Qvalue < 0.05 & abs(cell.proteins$AVG.Log2.Ratio) > 0.20), "Genes"]
    
# Load differential secreted protein results
secreted.proteins <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Secreted_Proteome_Differential_Abundance_Results.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Split concatenated gene names
    # The warning message that " 'as.is' should be specified by the caller" does not affect results and can be ignored
    secreted.proteins <- as.data.frame(splitstackshape::cSplit(secreted.proteins, splitCols = "Genes", sep = ";", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE)) # warning about specifying 'as.is' by the caller can be ignored

    # Remove duplicate gene names
    secreted.proteins <- secreted.proteins[!duplicated(secreted.proteins$Genes), ]
    
    # Extract gene names for significant proteins
    secreted.proteins.sig <- secreted.proteins[which(secreted.proteins$Qvalue < 0.05 & abs(secreted.proteins$AVG.Log2.Ratio) > 0.58), "Genes"]

# Generate Venn diagram object
my.venn <- list(Genes = my.transcriptome.sig,
                Cell_Proteins = cell.proteins.sig,
                Secreted_Proteins = secreted.proteins.sig
                )

# Generate plot and save
pdf(file = paste(dir.output, "Venn_Diagram_Overlap_Transcriptome_Cell-Proteins_Secretome", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(venn(my.venn))
    #plot(euler(my.venn, shape = "ellipse"), quantities = TRUE)

dev.off()

# Find genes/proteins in common across three conditions (can be listed in a figure)
shared.all <- intersect(my.transcriptome.sig, intersect(cell.proteins.sig, secreted.proteins.sig))



# SESSION INFO  ----------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/2_Heatmaps_and_Comparisons/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_Venn_Diagrams_and_Mosaic_Plots.txt", sep =""))
sessionInfo()
sink()      

    
# Clean the environment
rm(list=ls())

