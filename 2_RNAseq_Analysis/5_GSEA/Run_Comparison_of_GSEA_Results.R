# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Run_GSEA_Functions.R") 

# Load libraries
library(ggplot2) # for bubble plot
library(scales) # for modifying the ggplot colorbar
library(reshape2) # for FC & pval bubble plot
library(DOSE)
library(poolr) # To calculate meta pvalues

# Load GSEA results
load(file = '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T24_AluJb/All_GSEA_Results_FDR5_T24_AluJb_IMR90.R')
load(file = '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T72_Proteomics/All_GSEA_Results_FDR5_Proteomics_AluJb_IMR90.R')
load(file = '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T72_Lamivudine/All_GSEA_Results_FDR5_T72_Lamivudine_Experiment.R')
load(file = '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T72_WI38/All_GSEA_Results_FDR5_T72_WI38.R')

    

# IMR90 Proteomic Transcriptome vs IMR90 T24 --------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_Overlap/Proteomics_IMR90_vs_T24_IMR90/'

# Define plot labels
my.label.one <- 'Proteomic_Transcriptome'
my.label.two <- 'T24 IMR-90'


# Find overlapping gene sets   
overlap.Hallmark <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Hallmark.Proteomics@result,
                                         result.two = Hallmark.T24@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Hallmark")

overlap.Reactome <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Reactome.Proteomics@result,
                                         result.two = Reactome.T24@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Reactome")

overlap.GOBP <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                     result.one = GOBP.Proteomics@result,
                                     result.two = GOBP.T24@result,
                                     label.one = my.label.one,
                                     label.two = my.label.two,
                                     number_to_plot = 10,
                                     dir.output = my.output.dir,
                                     gs.label = "GOBP")

overlap.TE_families <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                           result.one = Family.Proteomics@result,
                                           result.two = Family.T24@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "TEs")


    
# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
    
# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Hallmark", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Hallmark, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Hallmark Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

  
# Reactome
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
  

# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


# TE families  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TE_families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()



# IMR90 Proteomic Transcriptome vs IMR90 T72 --------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_Overlap/Proteomics_IMR90_vs_T72_AluEff_Veh/'

# Define plot labels
my.label.one <- 'Proteomic_Transcriptome'
my.label.two <- 'T72 AluEff Veh'


# Find overlapping gene sets   
overlap.Hallmark <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Hallmark.Proteomics@result,
                                         result.two = Hallmark.AluEff.Veh@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Hallmark")

overlap.Reactome <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Reactome.Proteomics@result,
                                         result.two = Reactome.AluEff.Veh@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Reactome")

overlap.GOBP <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                     result.one = GOBP.Proteomics@result,
                                     result.two = GOBP.AluEff.Veh@result,
                                     label.one = my.label.one,
                                     label.two = my.label.two,
                                     number_to_plot = 10,
                                     dir.output = my.output.dir,
                                     gs.label = "GOBP")


overlap.TE_families <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                           result.one = Family.Proteomics@result,
                                           result.two = Family.AluEff.Veh@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "TEs")



# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
    
# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Hallmark", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Hallmark, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Hallmark Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

  
# Reactome
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
  

# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


# TE families  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TE_families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
    


# IMR90 Proteomic Transcriptome vs WI38 T72 --------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_Overlap/Proteomics_IMR90_vs_T72_WI38/'

# Define plot labels
my.label.one <- 'Proteomic_Transcriptome'
my.label.two <- 'T72 WI-38'


# Find overlapping gene sets   
overlap.Hallmark <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Hallmark.Proteomics@result,
                                         result.two = Hallmark.WI38@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Hallmark")

overlap.Reactome <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Reactome.Proteomics@result,
                                         result.two = Reactome.WI38@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Reactome")

overlap.GOBP <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                     result.one = GOBP.Proteomics@result,
                                     result.two = GOBP.WI38@result,
                                     label.one = my.label.one,
                                     label.two = my.label.two,
                                     number_to_plot = 10,
                                     dir.output = my.output.dir,
                                     gs.label = "GOBP")

overlap.TE_families <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                           result.one = Family.Proteomics@result,
                                           result.two = Family.WI38@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "TEs")



# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
    
# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Hallmark", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Hallmark, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Hallmark Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

  
# Reactome
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
  

# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


# TE families  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TE_families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
    


# IMR90 T72 vs WI38 T72 --------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_Overlap/T72_AluEff_Veh_IMR90_vs_T72_WI38/'

# Define plot labels
my.label.one <- 'T72 AluEff Veh'
my.label.two <- 'T72 WI-38'


# Find overlapping gene sets   
overlap.Hallmark <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Hallmark.AluEff.Veh@result,
                                         result.two = Hallmark.WI38@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Hallmark")

overlap.Reactome <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Reactome.AluEff.Veh@result,
                                         result.two = Reactome.WI38@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Reactome")

overlap.GOBP <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                     result.one = GOBP.AluEff.Veh@result,
                                     result.two = GOBP.WI38@result,
                                     label.one = my.label.one,
                                     label.two = my.label.two,
                                     number_to_plot = 10,
                                     dir.output = my.output.dir,
                                     gs.label = "GOBP")


overlap.TE_families <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                           result.one = Family.AluEff.Veh@result,
                                           result.two = Family.WI38@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "TEs")


    
# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
    
# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Hallmark", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Hallmark, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Hallmark Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

  
# Reactome
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
  

# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


# TE families  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TE_families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
    


# All 4 transcriptomes, no directionality --------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_Overlap/All_4_transcriptomes_no_directionality_filters/'

# Define plot labels
my.label.one <- 'Proteomic_Transcriptome'
my.label.two <- 'T24 IMR-90'
my.label.three <- 'T72 IMR-90'
my.label.four <- 'T72 WI-38'


# Find overlapping gene sets   
overlap.Hallmark <- Overlap_four_GSEA_res(comparison_directionality = c('none', 'none', 'none'),
                                         result.one = Hallmark.Proteomics@result,
                                         result.two = Hallmark.T24@result,
                                         result.three = Hallmark.AluEff.Veh@result,
                                         result.four = Hallmark.WI38@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         label.three = my.label.three,
                                         label.four = my.label.four,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Hallmark")

overlap.Reactome <- Overlap_four_GSEA_res(comparison_directionality = c('none', 'none', 'none'),
                                         result.one = Reactome.Proteomics@result,
                                         result.two = Reactome.T24@result,
                                         result.three = Reactome.AluEff.Veh@result,
                                         result.four = Reactome.WI38@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         label.three = my.label.three,
                                         label.four = my.label.four,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Reactome")

overlap.GOBP <- Overlap_four_GSEA_res(comparison_directionality = c('none', 'none', 'none'),
                                     result.one = GOBP.Proteomics@result,
                                     result.two = GOBP.T24@result,
                                     result.three = GOBP.AluEff.Veh@result,
                                     result.four = GOBP.WI38@result,
                                     label.one = my.label.one,
                                     label.two = my.label.two,
                                     label.three = my.label.three,
                                     label.four = my.label.four,
                                     number_to_plot = 10,
                                     dir.output = my.output.dir,
                                     gs.label = "GOBP")

overlap.TE_families <- Overlap_four_GSEA_res(comparison_directionality = c('none', 'none', 'none'),
                                             result.one = Family.Proteomics@result,
                                             result.two = Family.T24@result,
                                             result.three = Family.AluEff.Veh@result,
                                             result.four = Family.WI38@result,
                                             label.one = my.label.one,
                                             label.two = my.label.two,
                                             label.three = my.label.three,
                                             label.four = my.label.four,
                                             number_to_plot = 10,
                                             dir.output = my.output.dir,
                                             gs.label = "TEs")



# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
    
# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Hallmark", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Hallmark, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Hallmark Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

  
# Reactome
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
  

# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


# TE families  
pdf(paste(my.output.dir, "GSEA_BubblePlot_Overlap_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TE_families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
    


# END CODE --------------------------------------------------------------------------------------------------
  
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_Overlapping_Transcriptomes.txt", sep =""))
sessionInfo()
sink()      

    
# Clean the environment
rm(list=ls())      