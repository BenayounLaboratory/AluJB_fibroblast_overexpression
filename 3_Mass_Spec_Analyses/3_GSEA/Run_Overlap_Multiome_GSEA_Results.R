# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Run_GSEA_Functions.R") # ***NOTE: The current script relies on functions prepared in the "2_RNAseq_Analysis" analysis folder

# Load libraries
library(ggplot2) # for bubble plot
library(scales) # for modifying the ggplot colorbar
library(reshape2) # for FC & pval bubble plot
library(DOSE)
library(poolr) # To calculate meta pvalues

# Load GSEA results
load(file = '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T72_Proteomics/All_GSEA_Results_FDR5_Proteomics_AluJb_IMR90.R')
load(file = '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/3_GSEA/Results_Cellular_Proteomics/All_GSEA_Results_FDR5_Cell_Proteome.R')
load(file = '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/3_GSEA/Results_Secreted_Proteomics/All_GSEA_Results_FDR5_Secreted_Proteome.R')


# Multi-ome, no directionality --------------------------------------------------------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/3_GSEA/Results_Overlaps/'

# Define plot labels
my.label.one <- 'RNA'
my.label.two <- 'Cell Protein'
my.label.three <- 'Secreted Protein'


# Find overlapping gene sets   
overlap.Sen_Age <- Overlap_three_GSEA_res_unfiltered(comparison_directionality = c('none', 'none'),
                                                     result.one = Aging_Sen.Proteomics@result,
                                                     result.two = sen_aging.cell.proteome@result,
                                                     result.three = sen_aging.secreted.proteome@result,
                                                     label.one = my.label.one,
                                                     label.two = my.label.two,
                                                     label.three = my.label.three,
                                                     number_to_plot = 10,
                                                     dir.output = my.output.dir,
                                                     gs.label = "Senescence_Aging")



# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
    
# Senescence/Aging Gene Set Collection
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Senescence_Aging", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Sen_Age, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-6,0,6)), guide = "colorbar", limits=c(-6.1, 6.1), aesthetics = 'color') +
             ggtitle("Senescence/Aging Gene Sets GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()



# END CODE --------------------------------------------------------------------------------------------------
  
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/3_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_Overlapping_Multiomics.txt", sep =""))
sessionInfo()
sink()      
    
    
# Clean the environment
rm(list=ls())      