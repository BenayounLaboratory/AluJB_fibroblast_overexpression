# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(readxl) # To read excel
library(eulerr) # for making venn diagrams

# Define the output directories
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/4_Heatmaps_and_Comparisons/Plots/'


# AluJb OE (Genes from Multiomic Samples) vs Core SASP  ----------------------------------------------------------------------------------------


# Load Basisty Core SASP table
Basisty.core <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Aging_and_Senescence_Gene_Lists/Basisty_2020_Core_SASP_CLEANED_UP_gProfiler.csv", header = TRUE, sep = ',')

# Load differential expression results
DEGs.all <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_All.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
DEGs.sig <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Only keep genes (remove TEs)
    DEGs.all <- DEGs.all[grepl('ENSG', rownames(DEGs.all)), ]
    DEGs.sig <- DEGs.sig[grepl('ENSG', rownames(DEGs.sig)), ]
    
    # make dataframe with non-sig genes/proteins only
    DEGs.nonsig <- DEGs.all[setdiff(rownames(DEGs.all), rownames(DEGs.sig)), ]
    
# Set a seed
set.seed(90280)

# Calculate the size of each set
size.overlap <- sum(rownames(DEGs.sig) %in% Basisty.core$converted_alias)
size.SASP <- length(Basisty.core$converted_alias) - size.overlap
size.Alu <- nrow(DEGs.sig) - size.overlap
size.universe <- length(union(rownames(DEGs.all), Basisty.core$converted_alias)) # universe = all expressed genes + genes in senescence/aging list


# Run stats on the overlap

  # Generate a contingency table 
  my.contingency <- matrix(c(size.universe - sum(size.overlap, size.SASP, size.Alu), size.Alu, size.SASP, size.overlap), nrow=2)
  
  # Run Fischer's Exact Test for enrichment of overlap
  my.Fischer.res <- fisher.test(my.contingency, alternative = "greater")

  
# MAKE VENN DIAGRAMS 
  
# Define the element for the plot
my.overlaps <- c('Core_SASP' = size.SASP, 'AluOE' = size.Alu, 'Core_SASP&AluOE' = size.overlap)

# Generate the plot
pdf(file = paste(dir.output, "Venn_Diagram_Alu_Transcriptomics_vs_CORE_SASP", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(euler(my.overlaps), quantities = TRUE, main = paste('p = ', signif(my.Fischer.res$p.value, 3), sep = ''))

dev.off()


# MAKE MOSAIC PLOT

# calculate missing values
nonsig_in_list <- sum(rownames(DEGs.nonsig) %in% Basisty.core$converted_alias)
nonsig_not_in_list <- nrow(DEGs.nonsig) - nonsig_in_list

# make matrix to hold the data
my.mosaic <- matrix(c(size.overlap, nonsig_in_list, size.Alu, nonsig_not_in_list), 2, 2)
rownames(my.mosaic) <- c('Sig', 'Nonsignificant')
colnames(my.mosaic) <- c('Yes', 'No')

# Run Fischer's Exact Test comparing frequencies
my.Fischer.mosaic <- fisher.test(my.mosaic)

# Generate the plot
pdf(file = paste(dir.output, "Mosaic_plot_Alu_Transcriptomics_vs_CORE_SASP", '.pdf', sep=""), width = 5, height = 5)

    mosaicplot(my.mosaic,
               color = TRUE,
               xlab = "Transcriptome",
               ylab = "Core SASP",
               main = paste('Fishers exact p-value = ', signif(my.Fischer.mosaic$p.value, 3), sep = '')
               )

dev.off()



# Aging Fibroblast DEGs vs AluJb OE Genes (Multiomic samples)  ----------------------------------------------------------------------------------------


# Load AluJb OE DEGs
Alu.sig <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Separate up and down genes
    Alu.sig.up <- Alu.sig[which(Alu.sig$log2FoldChange > 0), ]
    Alu.sig.down <- Alu.sig[which(Alu.sig$log2FoldChange < 0), ]
    
    # Separate genes vs repeats
    Alu.repeat.sig.up <- Alu.sig.up[!grepl('ENSG', rownames(Alu.sig.up)), ]
    Alu.repeat.sig.down <- Alu.sig.down[!grepl('ENSG', rownames(Alu.sig.down)), ]
    Alu.DEGs.sig.up <- Alu.sig.up[grepl('ENSG', rownames(Alu.sig.up)), ]
    Alu.DEGs.sig.down <- Alu.sig.down[grepl('ENSG', rownames(Alu.sig.down)), ]
    
    # Remove variables that aren't needed
    rm(Alu.sig, Alu.sig.up, Alu.sig.down)
    
# Load differential expression results
aging.all <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/DESeq_Results/TETranscripts_GSE113957_Linear_Across_Age_All_Genes.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
aging.sig <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/DESeq_Results/TETranscripts_GSE113957_Linear_Across_Age_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Split all expressed genes and repeats (to use as universes)
    aging.all.repeats <- aging.all[!grepl('ENSG', rownames(aging.all)), ]
    aging.all.genes <- aging.all[grepl('ENSG', rownames(aging.all)), ]
    
    # Separate up and down genes/repeats
    aging.all.repeats.up <- aging.all.repeats[which(aging.all.repeats$log2FoldChange > 0), ]
    aging.all.repeats.down <- aging.all.repeats[which(aging.all.repeats$log2FoldChange < 0), ]
    aging.all.genes.up <- aging.all.genes[which(aging.all.genes$log2FoldChange > 0), ]
    aging.all.genes.down <- aging.all.genes[which(aging.all.genes$log2FoldChange < 0), ]

    # Split all sig genes and repeats 
    aging.sig.repeats <- aging.sig[!grepl('ENSG', rownames(aging.sig)), ]
    aging.sig.genes <- aging.sig[grepl('ENSG', rownames(aging.sig)), ]
    
    # Only keep sig genes/repeats split by up/down change
    aging.sig.repeats.up <- aging.sig.repeats[which(aging.sig.repeats$log2FoldChange > 0), ]
    aging.sig.repeats.down <- aging.sig.repeats[which(aging.sig.repeats$log2FoldChange < 0), ]
    aging.sig.genes.up <- aging.sig.genes[which(aging.sig.genes$log2FoldChange > 0), ]
    aging.sig.genes.down <- aging.sig.genes[which(aging.sig.genes$log2FoldChange < 0), ]
    
    # make dataframe with non-sig genes/repeats only
    aging.nonsig.repeats.up <- aging.all.repeats.up[setdiff(rownames(aging.all.repeats.up), rownames(aging.sig.repeats.up)), ]
    aging.nonsig.repeats.down <- aging.all.repeats.down[setdiff(rownames(aging.all.repeats.down), rownames(aging.sig.repeats.down)), ]
    aging.nonsig.genes.up <- aging.all.genes.up[setdiff(rownames(aging.all.genes.up), rownames(aging.sig.genes.up)), ]
    aging.nonsig.genes.down <- aging.all.genes.down[setdiff(rownames(aging.all.genes.down), rownames(aging.sig.genes.down)), ]
    
    # Remove variables that aren't needed
    rm(aging.all, aging.sig, aging.all.repeats, aging.all.genes, aging.sig.repeats, aging.sig.genes)
    
    
    
    
    
# COMPARE REPEATS THAT ARE UP IN BOTH AGING AND ALU OE
    
# Set a seed
set.seed(90280)

# Calculate the size of each set *in the plot
size.overlap <- sum(rownames(aging.sig.repeats.up) %in% rownames(Alu.repeat.sig.up))
size.Alu <- nrow(Alu.repeat.sig.up) - size.overlap
size.Aging <- nrow(aging.sig.repeats.up) - size.overlap
size.universe <- length(union(rownames(aging.all.repeats.up), rownames(Alu.repeat.sig.up))) # universe = all aging upregulated repeats + repeats in Alu list

# Define the element to plot
my.overlaps <- c('Alu_up_TE' = size.Alu, 'Aging_up_TE' = size.Aging, 'Alu_up_TE&Aging_up_TE' = size.overlap)

# Run stats on the overlap

  # Generate a contingency table 
  my.contingency <- matrix(c(size.universe - sum(size.overlap, size.Alu, size.Aging), size.Aging, size.Alu, size.overlap), nrow=2)
  
  # Run Fischer's Exact Test for enrichment of overlap
  my.Fischer.res <- fisher.test(my.contingency, alternative = "greater")

  
# Generate the plot
pdf(file = paste(dir.output, "Venn_Diagram_Aging_Transcriptomic_Repeats_Up_vs_Alu_Up", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(euler(my.overlaps), quantities = TRUE, main = paste('p = ', signif(my.Fischer.res$p.value, 3), sep = ''))

dev.off()


# MAKE MOSAIC PLOT

# calculate missing values
nonsig_in_list <- sum(rownames(aging.nonsig.repeats.up) %in% rownames(Alu.repeat.sig.up))
nonsig_not_in_list <- nrow(aging.nonsig.repeats.up) - nonsig_in_list

# make matrix to hold the data
my.mosaic <- matrix(c(size.overlap, nonsig_in_list, size.Aging, nonsig_not_in_list), 2, 2)
rownames(my.mosaic) <- c('Sig', 'Nonsignificant')
colnames(my.mosaic) <- c('Yes', 'No')

# Run Fischer's Exact Test comparing frequencies
my.Fischer.mosaic <- fisher.test(my.mosaic)

# Generate the plot
pdf(file = paste(dir.output, "Mosaic_plot_Aging_Transcriptome_Repeats_Up_vs_AluOE_Repeats_Up", '.pdf', sep=""), width = 5, height = 5)

    mosaicplot(my.mosaic,
               color = TRUE,
               xlab = "Aging Upregulated Repeats",
               ylab = "Alu Upregulated Repeats",
               main = paste('Fishers exact p-value = ', signif(my.Fischer.mosaic$p.value, 3), sep = '')
               )

dev.off()





# COMPARE REPEATS THAT ARE DOWN IN BOTH AGING AND ALU OE
    
# Set a seed
set.seed(90280)

# Calculate the size of each set *in the plot
size.overlap <- sum(rownames(aging.sig.repeats.down) %in% rownames(Alu.repeat.sig.down))
size.Alu <- nrow(Alu.repeat.sig.down) - size.overlap
size.Aging <- nrow(aging.sig.repeats.down) - size.overlap
size.universe <- length(union(rownames(aging.all.repeats.down), rownames(Alu.repeat.sig.down))) # universe = all aging downregulated repeats + repeats in Alu list

# Define the element to plot
my.overlaps <- c('Alu_down_TE' = size.Alu, 'Aging_down_TE' = size.Aging, 'Alu_down_TE&Aging_down_TE' = size.overlap)

# Run stats on the overlap

  # Generate a contingency table 
  my.contingency <- matrix(c(size.universe - sum(size.overlap, size.Alu, size.Aging), size.Aging, size.Alu, size.overlap), nrow=2)
  
  # Run Fischer's Exact Test for enrichment of overlap
  my.Fischer.res <- fisher.test(my.contingency, alternative = "greater")

  
# Generate the plot
pdf(file = paste(dir.output, "Venn_Diagram_Aging_Transcriptomic_Repeats_Down_vs_Alu_Down", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(euler(my.overlaps), quantities = TRUE, main = paste('p = ', signif(my.Fischer.res$p.value, 3), sep = ''))

dev.off()


# MAKE MOSAIC PLOT

# calculate missing values
nonsig_in_list <- sum(rownames(aging.nonsig.repeats.down) %in% rownames(Alu.repeat.sig.down))
nonsig_not_in_list <- nrow(aging.nonsig.repeats.down) - nonsig_in_list

# make matrix to hold the data
my.mosaic <- matrix(c(size.overlap, nonsig_in_list, size.Aging, nonsig_not_in_list), 2, 2)
rownames(my.mosaic) <- c('Sig', 'Nonsignificant')
colnames(my.mosaic) <- c('Yes', 'No')

# Run Fischer's Exact Test comparing frequencies
my.Fischer.mosaic <- fisher.test(my.mosaic)

# Generate the plot
pdf(file = paste(dir.output, "Mosaic_plot_Aging_Transcriptome_Repeats_Down_vs_AluOE_Repeats_Down", '.pdf', sep=""), width = 5, height = 5)

    mosaicplot(my.mosaic,
               color = TRUE,
               xlab = "Aging Downregulated Repeats",
               ylab = "Alu Downregulated Repeats",
               main = paste('Fishers exact p-value = ', signif(my.Fischer.mosaic$p.value, 3), sep = '')
               )

dev.off()





# COMPARE GENES THAT ARE UP IN BOTH AGING AND ALU OE
    
# Set a seed
set.seed(90280)

# Calculate the size of each set *in the plot
size.overlap <- sum(rownames(aging.sig.genes.up) %in% rownames(Alu.DEGs.sig.up))
size.Alu <- nrow(Alu.DEGs.sig.up) - size.overlap
size.Aging <- nrow(aging.sig.genes.up) - size.overlap
size.universe <- length(union(rownames(aging.all.genes.up), rownames(Alu.DEGs.sig.up))) # universe = all aging upregulated genes + genes in Alu list

# Define the element to plot
my.overlaps <- c('Alu_up_DEGs' = size.Alu, 'Aging_up_DEGs' = size.Aging, 'Alu_up_DEGs&Aging_up_DEGs' = size.overlap)

# Run stats on the overlap

  # Generate a contingency table 
  my.contingency <- matrix(c(size.universe - sum(size.overlap, size.Alu, size.Aging), size.Aging, size.Alu, size.overlap), nrow=2)
  
  # Run Fischer's Exact Test for enrichment of overlap
  my.Fischer.res <- fisher.test(my.contingency, alternative = "greater")

  
# Generate the plot
pdf(file = paste(dir.output, "Venn_Diagram_Aging_Transcriptomic_Genes_Up_vs_Alu_Up", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(euler(my.overlaps), quantities = TRUE, main = paste('p = ', signif(my.Fischer.res$p.value, 3), sep = ''))

dev.off()


# MAKE MOSAIC PLOT

# calculate missing values
nonsig_in_list <- sum(rownames(aging.nonsig.genes.up) %in% rownames(Alu.DEGs.sig.up))
nonsig_not_in_list <- nrow(aging.nonsig.genes.up) - nonsig_in_list

# make matrix to hold the data
my.mosaic <- matrix(c(size.overlap, nonsig_in_list, size.Aging, nonsig_not_in_list), 2, 2)
rownames(my.mosaic) <- c('Sig', 'Nonsignificant')
colnames(my.mosaic) <- c('Yes', 'No')

# Run Fischer's Exact Test comparing frequencies
my.Fischer.mosaic <- fisher.test(my.mosaic)

# Generate the plot
pdf(file = paste(dir.output, "Mosaic_plot_Aging_Transcriptome_Genes_Up_vs_AluOE_Genes_Up", '.pdf', sep=""), width = 5, height = 5)

    mosaicplot(my.mosaic,
               color = TRUE,
               xlab = "Aging Upregulated Genes",
               ylab = "Alu Upregulated Genes",
               main = paste('Fishers exact p-value = ', signif(my.Fischer.mosaic$p.value, 3), sep = '')
               )

dev.off()





# COMPARE GENES THAT ARE DOWN IN BOTH AGING AND ALU OE
    
# Set a seed
set.seed(90280)

# Calculate the size of each set *in the plot
size.overlap <- sum(rownames(aging.sig.genes.down) %in% rownames(Alu.DEGs.sig.down))
size.Alu <- nrow(Alu.DEGs.sig.down) - size.overlap
size.Aging <- nrow(aging.sig.genes.down) - size.overlap
size.universe <- length(union(rownames(aging.all.genes.down), rownames(Alu.DEGs.sig.down))) # universe = all aging downregulated genes + genes in Alu list

# Define the element to plot
my.overlaps <- c('Alu_down_DEGs' = size.Alu, 'Aging_down_DEGs' = size.Aging, 'Alu_down_DEGs&Aging_down_DEGs' = size.overlap)

# Run stats on the overlap

  # Generate a contingency table 
  my.contingency <- matrix(c(size.universe - sum(size.overlap, size.Alu, size.Aging), size.Aging, size.Alu, size.overlap), nrow=2)
  
  # Run Fischer's Exact Test for enrichment of overlap
  my.Fischer.res <- fisher.test(my.contingency, alternative = "greater")

  
# Generate the plot
pdf(file = paste(dir.output, "Venn_Diagram_Aging_Transcriptomic_Genes_Down_vs_Alu_Down", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(euler(my.overlaps), quantities = TRUE, main = paste('p = ', signif(my.Fischer.res$p.value, 3), sep = ''))

dev.off()


# MAKE MOSAIC PLOT

# calculate missing values
nonsig_in_list <- sum(rownames(aging.nonsig.genes.down) %in% rownames(Alu.DEGs.sig.down))
nonsig_not_in_list <- nrow(aging.nonsig.genes.down) - nonsig_in_list

# make matrix to hold the data
my.mosaic <- matrix(c(size.overlap, nonsig_in_list, size.Aging, nonsig_not_in_list), 2, 2)
rownames(my.mosaic) <- c('Sig', 'Nonsignificant')
colnames(my.mosaic) <- c('Yes', 'No')

# Run Fischer's Exact Test comparing frequencies
my.Fischer.mosaic <- fisher.test(my.mosaic)

# Generate the plot
pdf(file = paste(dir.output, "Mosaic_plot_Aging_Transcriptome_Genes_Down_vs_AluOE_Genes_Down", '.pdf', sep=""), width = 5, height = 5)

    mosaicplot(my.mosaic,
               color = TRUE,
               xlab = "Aging Downregulated Genes",
               ylab = "Alu Downregulated Genes",
               main = paste('Fishers exact p-value = ', signif(my.Fischer.mosaic$p.value, 3), sep = '')
               )

dev.off()



# AluS DEGs (Cantarella 2019) vs AluJb OE Genes (Multiomic samples)  ----------------------------------------------------------------------------------------


# Load AluS OE DEGs
AluSx.sig <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Cantarella_2019_AluS_OE_Gene_Lists/Cantarella_2019_AluSx_DEGs.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
AluSq.sig <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Cantarella_2019_AluS_OE_Gene_Lists/Cantarella_2019_AluSq2_DEGs.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Aggregate into 1 list
    Alu.sig <- rbind(AluSx.sig, AluSq.sig)
    
    # Remove transcript or ENSEMBL version info. Keeping this can raise issues with gene name mapping.
    Alu.sig$Ensembl.gene.ID <- sub("\\..*", "", Alu.sig$Ensembl.gene.ID, fixed=FALSE)
  
    # Separate up and down genes
    Alu.DEGs.sig.up <- Alu.sig[which(Alu.sig$log2FoldChange > 0), ]
    Alu.DEGs.sig.down <- Alu.sig[which(Alu.sig$log2FoldChange < 0), ]
    
    # Remove duplicates
    Alu.DEGs.sig.up <- Alu.DEGs.sig.up[which(!duplicated(Alu.DEGs.sig.up$Ensembl.gene.ID)), ]
    Alu.DEGs.sig.down <- Alu.DEGs.sig.down[which(!duplicated(Alu.DEGs.sig.down$Ensembl.gene.ID)), ]

    # Remove variables that aren't needed
    rm(AluSx.sig, AluSq.sig, Alu.sig)
    
# Load differential expression results
my.all <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_All.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
my.sig <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Split all expressed genes and repeats (to use as universes)
    my.all.genes <- my.all[grepl('ENSG', rownames(my.all)), ]
    
    # Separate up and down genes/repeats
    my.all.genes.up <- my.all.genes[which(my.all.genes$log2FoldChange > 0), ]
    my.all.genes.down <- my.all.genes[which(my.all.genes$log2FoldChange < 0), ]

    # Split all sig genes and repeats 
    my.sig.genes <- my.sig[grepl('ENSG', rownames(my.sig)), ]
    
    # Only keep sig genes/repeats split by up/down change
    my.sig.genes.up <- my.sig.genes[which(my.sig.genes$log2FoldChange > 0), ]
    my.sig.genes.down <- my.sig.genes[which(my.sig.genes$log2FoldChange < 0), ]
    
    # make dataframe with non-sig genes/repeats only
    my.nonsig.genes.up <- my.all.genes.up[setdiff(rownames(my.all.genes.up), rownames(my.sig.genes.up)), ]
    my.nonsig.genes.down <- my.all.genes.down[setdiff(rownames(my.all.genes.down), rownames(my.sig.genes.down)), ]
    
    # Remove variables that aren't needed
    rm(my.all, my.sig, my.all.genes, my.sig.genes)
    





# COMPARE GENES THAT ARE UP IN BOTH MY AND CANTARELLA'S ALU 2019
    
# Set a seed
set.seed(90280)

# Calculate the size of each set *in the plot
size.overlap <- sum(rownames(my.sig.genes.up) %in% Alu.DEGs.sig.up$Ensembl.gene.ID)
size.Alu <- nrow(Alu.DEGs.sig.up) - size.overlap
size.my <- nrow(my.sig.genes.up) - size.overlap
size.universe <- length(union(rownames(my.all.genes.up), Alu.DEGs.sig.up$Ensembl.gene.ID)) # universe = all my upregulated genes + genes in Alu list

# Define the element to plot
my.overlaps <- c('Alu_up_DEGs' = size.Alu, 'my_up_DEGs' = size.my, 'Alu_up_DEGs&my_up_DEGs' = size.overlap)

# Run stats on the overlap

  # Generate a contingency table 
  my.contingency <- matrix(c(size.universe - sum(size.overlap, size.Alu, size.my), size.my, size.Alu, size.overlap), nrow=2)
  
  # Run Fischer's Exact Test for enrichment of overlap
  my.Fischer.res <- fisher.test(my.contingency, alternative = "greater")

  
# Generate the plot
pdf(file = paste(dir.output, "Venn_Diagram_AluJb_Transcriptomic_Genes_Up_vs_Cantarella_2019_Genes_Up", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(euler(my.overlaps), quantities = TRUE, main = paste('p = ', signif(my.Fischer.res$p.value, 3), sep = ''))

dev.off()


# MAKE MOSAIC PLOT

# calculate missing values
nonsig_in_list <- sum(rownames(my.nonsig.genes.up) %in% Alu.DEGs.sig.up$Ensembl.gene.ID)
nonsig_not_in_list <- nrow(my.nonsig.genes.up) - nonsig_in_list

# make matrix to hold the data
my.mosaic <- matrix(c(size.overlap, nonsig_in_list, size.my, nonsig_not_in_list), 2, 2)
rownames(my.mosaic) <- c('Sig', 'Nonsignificant')
colnames(my.mosaic) <- c('Yes', 'No')

# Run Fischer's Exact Test comparing frequencies
my.Fischer.mosaic <- fisher.test(my.mosaic)

# Generate the plot
pdf(file = paste(dir.output, "Mosaic_plot_AluJb_Transcriptome_Genes_Up_vs_Cantarella_2019_Genes_Up", '.pdf', sep=""), width = 5, height = 5)

    mosaicplot(my.mosaic,
               color = TRUE,
               xlab = "AluJb Upregulated Genes",
               ylab = "AluS Upregulated Genes",
               main = paste('Fishers exact p-value = ', signif(my.Fischer.mosaic$p.value, 3), sep = '')
               )

dev.off()






# COMPARE GENES THAT ARE DOWN IN BOTH MY AND CANTARELLA'S ALU 2019
    
# Set a seed
set.seed(90280)

# Calculate the size of each set *in the plot
size.overlap <- sum(rownames(my.sig.genes.down) %in% Alu.DEGs.sig.down$Ensembl.gene.ID)
size.Alu <- nrow(Alu.DEGs.sig.down) - size.overlap
size.my <- nrow(my.sig.genes.down) - size.overlap
size.universe <- length(union(rownames(my.all.genes.down), Alu.DEGs.sig.down$Ensembl.gene.ID)) # universe = all my downregulated genes + genes in Alu list

# Define the element to plot
my.overlaps <- c('Alu_down_DEGs' = size.Alu, 'my_down_DEGs' = size.my, 'Alu_down_DEGs&my_down_DEGs' = size.overlap)

# Run stats on the overlap

  # Generate a contingency table 
  my.contingency <- matrix(c(size.universe - sum(size.overlap, size.Alu, size.my), size.my, size.Alu, size.overlap), nrow=2)
  
  # Run Fischer's Exact Test for enrichment of overlap
  my.Fischer.res <- fisher.test(my.contingency, alternative = "greater")

  
# Generate the plot
pdf(file = paste(dir.output, "Venn_Diagram_AluJb_Transcriptomic_Genes_Down_vs_Cantarella_2019_Genes_Down", '.pdf', sep=""), width = 4, height = 3)

    # Make the plot
    plot(euler(my.overlaps), quantities = TRUE, main = paste('p = ', signif(my.Fischer.res$p.value, 3), sep = ''))

dev.off()


# MAKE MOSAIC PLOT

# calculate missing values
nonsig_in_list <- sum(rownames(my.nonsig.genes.down) %in% Alu.DEGs.sig.down$Ensembl.gene.ID)
nonsig_not_in_list <- nrow(my.nonsig.genes.down) - nonsig_in_list

# make matrix to hold the data
my.mosaic <- matrix(c(size.overlap, nonsig_in_list, size.my, nonsig_not_in_list), 2, 2)
rownames(my.mosaic) <- c('Sig', 'Nonsignificant')
colnames(my.mosaic) <- c('Yes', 'No')

# Run Fischer's Exact Test comparing frequencies
my.Fischer.mosaic <- fisher.test(my.mosaic)

# Generate the plot
pdf(file = paste(dir.output, "Mosaic_plot_AluJb_Transcriptome_Genes_Down_vs_Cantarella_2019_Genes_Down", '.pdf', sep=""), width = 5, height = 5)

    mosaicplot(my.mosaic,
               color = TRUE,
               xlab = "AluJb Downregulated Genes",
               ylab = "AluS Downregulated Genes",
               main = paste('Fishers exact p-value = ', signif(my.Fischer.mosaic$p.value, 3), sep = '')
               )

dev.off()



# SESSION INFO  ----------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/4_Heatmaps_and_Comparisons/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_Venn_Diagrams_Mosaic_Plots.txt", sep =""))
sessionInfo()
sink()      

    
# Clean the environment
rm(list=ls())

