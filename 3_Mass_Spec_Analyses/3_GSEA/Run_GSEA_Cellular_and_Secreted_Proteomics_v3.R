# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Run_GSEA_Functions.R")

# Load libraries
library(clusterProfiler) # to run GSEA
library(ggplot2) # for split dotplots
library(scales) # for modifying the ggplot colorbar
library(enrichplot) # for gseaplot2

# Load gene sets 
load(file = '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Gene_Set_Collections_for_GSEA.R')

# PLEASE NOTE:
    # For most of the GSEA analyses below, a warning will appear about a small percentage of ranked genes/proteins tying for the same rank because the value for their ranking statistic is equal
    # The order of these genes is expected to slightly affect the resulting statistics but are unlikely to affect the top significant pathways or major themes identified
    # Keeping that in mind, the warning can be ignored


# CELLULAR PROTEOMICS ------------------------------------------------------------------------------------------------------------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/3_GSEA/Results_Cellular_Proteomics/'

# Load differential abundance results
diff.abundance <- read.csv('/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu//Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Cellular_Proteome_Differential_Abundance_Results.txt', header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

        # Split concatenated genes
        # The warning message that " 'as.is' should be specified by the caller" does not affect results and can be ignored
        diff.abundance <- as.data.frame(splitstackshape::cSplit(diff.abundance, splitCols = "Genes", sep = ";", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE)) # warning about specifying 'as.is' by the caller can be ignored
    
        # Remove duplicate gene names
        diff.abundance <- diff.abundance[!duplicated(diff.abundance$Genes), ]
        
        # Define genelists using the product of Qvalue and log2FC
        genelist.Cell_proteome <- diff.abundance$Rank_metric
        
        # Assign gene names to the list 
        names(genelist.Cell_proteome) <- diff.abundance$Genes
        
        # Sort the genelist, largest to smallest
        genelist.Cell_proteome = sort(genelist.Cell_proteome, decreasing = TRUE)


# Run GSEA using GO BP Collection
GOBP.cell.proteome <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.Cell_proteome, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'Cell_Proteome', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.cell.proteome <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.Cell_proteome, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'Cell_Proteome', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.cell.proteome <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.Cell_proteome, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'Cell_Proteome', output.object = 'filtered')

# Run GSEA using senescence/aging collection
sen_aging.cell.proteome <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Aging_Senescence/', my.genelist = genelist.Cell_proteome, my.gs.collection = Aging_Senescence.gs, gs.label = 'Aging_Senescence', condition_label = 'Cell_Proteome', output.object = 'unfiltered', my.max_gs_size = Inf, items_to_plot = 10, plot_non_sig = TRUE)

# Save results
save(GOBP.cell.proteome, Reactome.cell.proteome, Hallmark.cell.proteome, sen_aging.cell.proteome,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Cell_Proteome.R", sep = ''))



# SECRETED PROTEOMICS ------------------------------------------------------------------------------------------------------------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/3_GSEA/Results_Secreted_Proteomics/'

# Load differential abundance results
diff.abundance <- read.csv('/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Secreted_Proteome_Differential_Abundance_Results.txt', header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

        # Split concatenated genes
        # The warning message that " 'as.is' should be specified by the caller" does not affect results and can be ignored
        diff.abundance <- as.data.frame(splitstackshape::cSplit(diff.abundance, splitCols = "Genes", sep = ";", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE)) # warning about specifying 'as.is' by the caller can be ignored
    
        # Remove duplicate gene names
        diff.abundance <- diff.abundance[!duplicated(diff.abundance$Genes), ]
        
        # Define genelists using the product of Qvalue and log2FC
        genelist.secreted_proteome <- diff.abundance$Rank_metric
        
        # Assign gene names to the list 
        names(genelist.secreted_proteome) <- diff.abundance$Genes
        
        # Sort the genelist, largest to smallest
        genelist.secreted_proteome = sort(genelist.secreted_proteome, decreasing = TRUE)
        

# Run GSEA using GO BP Collection
GOBP.cell.proteome <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.secreted_proteome, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'Secreted_Proteome', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.cell.proteome <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.secreted_proteome, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'Secreted_Proteome', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.cell.proteome <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.secreted_proteome, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'Secreted_Proteome', output.object = 'filtered')

# Run GSEA using senescence/aging collection
sen_aging.secreted.proteome <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Aging_Senescence/', my.genelist = genelist.secreted_proteome, my.gs.collection = Aging_Senescence.gs, gs.label = 'Aging_Senescence', condition_label = 'Secreted_Proteome', output.object = 'unfiltered', my.max_gs_size = Inf, items_to_plot = 10, plot_non_sig = TRUE)

# Save results
save(GOBP.cell.proteome, Reactome.cell.proteome, Hallmark.cell.proteome, sen_aging.secreted.proteome,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Secreted_Proteome.R", sep = ''))



# SESSION INFO ------------------------------------------------------------------------------------------------------------------------------------------------------


# Save Session Info 
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/3_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_vs_Cell_and_Secreted_Proteome.txt", sep =""))
sessionInfo()
sink()      
    
    
# Clean the environment
rm(list=ls())
    