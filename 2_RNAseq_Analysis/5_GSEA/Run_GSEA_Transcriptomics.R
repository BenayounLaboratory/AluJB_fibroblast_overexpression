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
    # For most of the GSEA analyses below, a warning will appear about a small percentage of ranked genes tying for the same rank because the value for their ranking statistic is equal
    # The order of these genes is expected to slightly affect the resulting statistics but are unlikely to affect the top significant pathways or major themes identified
    # Keeping that in mind, the warning can be ignored


# Set 1: T24 AluJb ------------------------------------------------------------------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T24_AluJb/'

# Make gene lists
genelist.T24 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T24_AluJb_DEGs_All.txt")



# Run GSEA using GO BP Collection
GOBP.T24 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.T24, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'T24_IMR90', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.T24 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.T24, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'T24_IMR90', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.T24 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.T24, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'T24_IMR90', output.object = 'filtered')

# Run GSEA using Repeat Class Collection
Class.T24 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.T24, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'T24_IMR90', output.object = 'unfiltered')

# Run GSEA using Repeat Family Collection
Family.T24 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.T24, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'T24_IMR90', output.object = 'unfiltered')

    # Generate GSEA Plots for Alu family
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_T24_IMR90", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.T24, geneSetID = c("Alu subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()

# Run GSEA using Aging Fibroblast and senescence genesets
Aging_Sen.T24 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Aging_Senescence/', my.genelist = genelist.T24, my.gs.collection = Aging_Senescence.gs, gs.label = 'Aging_Senescence_Lists', condition_label = 'T24_IMR90', output.object = 'unfiltered', my.max_gs_size = Inf, items_to_plot = 10, plot_non_sig = TRUE)

# Save results
save(GOBP.T24, Reactome.T24, Hallmark.T24, Class.T24, Family.T24, Aging_Sen.T24,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_T24_AluJb_IMR90.R", sep = ''))



# Set 2: T72 Proteomics ------------------------------------------------------------------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T72_Proteomics/'

# Make gene lists
genelist.Proteomics <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_All.txt")



# Run GSEA using GO BP Collection
GOBP.Proteomics <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.Proteomics, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'Proteomics_IMR90', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.Proteomics <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.Proteomics, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'Proteomics_IMR90', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.Proteomics <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.Proteomics, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'Proteomics_IMR90', output.object = 'filtered')

# Run GSEA using Repeat Class Collection
Class.Proteomics <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.Proteomics, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'Proteomics_IMR90', output.object = 'unfiltered')

# Run GSEA using Repeat Family Collection
Family.Proteomics <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.Proteomics, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'Proteomics_IMR90', output.object = 'unfiltered')

    # Generate GSEA Plots for Alu family
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_Proteomics_IMR90", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.Proteomics, geneSetID = c("Alu subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()

# Run GSEA using Aging Fibroblast and senescence genesets
Aging_Sen.Proteomics <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Aging_Senescence/', my.genelist = genelist.Proteomics, my.gs.collection = Aging_Senescence.gs, gs.label = 'Aging_Senescence_Lists', condition_label = 'Proteomics_IMR90', output.object = 'unfiltered', my.max_gs_size = Inf, items_to_plot = 10, plot_non_sig = TRUE)

# Run GSEA using Cantarella 2019
Cantarella.Proteomics <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Cantarella_2019/', my.genelist = genelist.Proteomics, my.gs.collection = Cantarella.gs, gs.label = 'Cantarella_2019', condition_label = 'Proteomics_IMR90', output.object = 'unfiltered')

# Save results
save(GOBP.Proteomics, Reactome.Proteomics, Hallmark.Proteomics, Class.Proteomics, Family.Proteomics, Aging_Sen.Proteomics,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Proteomics_AluJb_IMR90.R", sep = ''))



# Set 3a: T72 Lamivudine ------------------------------------------------------------------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T72_Lamivudine/'

# Make gene lists
genelist.AluEff.Veh <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_LamExperiment_Alu_vs_Empty_Veh_DEGs_All.txt")
genelist.AluEff.3TC <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_LamExperiment_Alu_vs_Empty_3TC_DEGs_All.txt")
genelist.3TCEff.Empty <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_LamExperiment_3TC_vs_Veh_Empty_DEGs_All.txt")
genelist.3TCEff.Alu <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_LamExperiment_3TC_vs_Veh_Alu_DEGs_All.txt")


# Run GSEA using GO BP Collection
GOBP.AluEff.Veh <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.AluEff.Veh, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'Alu_Effect.Veh', output.object = 'filtered')
GOBP.AluEff.3TC <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.AluEff.3TC, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'Alu_Effect.3TC', output.object = 'filtered')
GOBP.3TCEff.Empty <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.3TCEff.Empty, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = '3TC_Effect.Empty', output.object = 'filtered')
GOBP.3TCEff.Alu <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.3TCEff.Alu, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = '3TC_Effect.Alu', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.AluEff.Veh <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.AluEff.Veh, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'Alu_Effect.Veh', output.object = 'filtered')
Reactome.AluEff.3TC <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.AluEff.3TC, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'Alu_Effect.3TC', output.object = 'filtered')
Reactome.3TCEff.Empty <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.3TCEff.Empty, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = '3TC_Effect.Empty', output.object = 'filtered')
Reactome.3TCEff.Alu <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.3TCEff.Alu, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = '3TC_Effect.Alu', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.AluEff.Veh <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.AluEff.Veh, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'Alu_Effect.Veh', output.object = 'filtered')
Hallmark.AluEff.3TC <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.AluEff.3TC, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'Alu_Effect.3TC', output.object = 'filtered')
Hallmark.3TCEff.Empty <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.3TCEff.Empty, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = '3TC_Effect.Empty', output.object = 'filtered')
Hallmark.3TCEff.Alu <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.3TCEff.Alu, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = '3TC_Effect.Alu', output.object = 'filtered')

# Run GSEA using Repeat Class Collection
Class.AluEff.Veh <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.AluEff.Veh, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'Alu_Effect.Veh', output.object = 'unfiltered')
Class.AluEff.3TC <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.AluEff.3TC, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'Alu_Effect.3TC', output.object = 'unfiltered')
Class.3TCEff.Empty <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.3TCEff.Empty, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = '3TC_Effect.Empty', output.object = 'unfiltered')
Class.3TCEff.Alu <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.3TCEff.Alu, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = '3TC_Effect.Alu', output.object = 'unfiltered')

# Run GSEA using Repeat Family Collection
Family.AluEff.Veh <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.AluEff.Veh, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'Alu_Effect.Veh', output.object = 'unfiltered')
Family.AluEff.3TC <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.AluEff.3TC, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'Alu_Effect.3TC', output.object = 'unfiltered')
Family.3TCEff.Empty <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.3TCEff.Empty, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = '3TC_Effect.Empty', output.object = 'unfiltered')
Family.3TCEff.Alu <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.3TCEff.Alu, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = '3TC_Effect.Alu', output.object = 'unfiltered')

    # Generate GSEA Plots for Alu family
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_T72_Alu_Eff_Veh", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.AluEff.Veh, geneSetID = c("Alu subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_T72_Alu_Eff_3TC", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.AluEff.3TC, geneSetID = c("Alu subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_T72_3TC_Eff_Empty", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.3TCEff.Empty, geneSetID = c("Alu subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_T72_3TC_Eff_Alu", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.3TCEff.Alu, geneSetID = c("Alu subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()

# Save results
save(GOBP.AluEff.Veh, Reactome.AluEff.Veh, Hallmark.AluEff.Veh, Class.AluEff.Veh, Family.AluEff.Veh,
     GOBP.AluEff.3TC, Reactome.AluEff.3TC, Hallmark.AluEff.3TC, Class.AluEff.3TC, Family.AluEff.3TC,
     GOBP.3TCEff.Empty, Reactome.3TCEff.Empty, Hallmark.3TCEff.Empty, Class.3TCEff.Empty, Family.3TCEff.Empty,
     GOBP.3TCEff.Alu, Reactome.3TCEff.Alu, Hallmark.3TCEff.Alu, Class.3TCEff.Alu, Family.3TCEff.Alu,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_T72_Lamivudine_Experiment.R", sep = ''))



# Set 3b: T72 Lamivudine Clustering ORA ------------------------------------------------------------------------------------------------------------


# directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T72_Lamivudine_ORA/'

# Load DESeq results 
results.all <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/LRT_LamExperiment_All_Groups_All_Genes.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # keep only genes
    results.genes <- results.all[grepl('ENSG', rownames(results.all)), ]
    
    # Assign EnsembleID to entries with an NA or Empty cell for HGNC Symbol
    results.genes[is.na(results.genes$hgnc_symbol), 'hgnc_symbol'] <- rownames(results.genes[is.na(results.genes$hgnc_symbol), ])
    results.genes[which(results.genes$hgnc_symbol == ''), 'hgnc_symbol'] <- rownames(results.genes[which(results.genes$hgnc_symbol == ''), ])
    
    # Remove duplicate symbols
    results.genes <- results.genes[!duplicated(results.genes$hgnc_symbol), ]
    
    # Assign symbol to rownames
    rownames(results.genes) <- results.genes$hgnc_symbol

# define the universe genes
universe.genes <- rownames(results.genes)

# Load DEGPatterns results
gene.clusters <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/LRT_LamExperiment_All_Groups_degPatterns_Clustering.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # keep only genes
    gene.clusters <- gene.clusters[grepl('ENSG', rownames(gene.clusters)), ]
    
    # Assign symbols to sig genes
    gene.clusters$hgnc_symbol <- results.all[gene.clusters$genes, 'hgnc_symbol']
    
    # Assign EnsembleID to entries with an NA or Empty cell for HGNC Symbol
    gene.clusters[is.na(gene.clusters$hgnc_symbol), 'hgnc_symbol'] <- rownames(gene.clusters[is.na(gene.clusters$hgnc_symbol), ])
    gene.clusters[which(gene.clusters$hgnc_symbol == ''), 'hgnc_symbol'] <- rownames(gene.clusters[which(gene.clusters$hgnc_symbol == ''), ])
    
    # Remove duplicate symbols
    gene.clusters <- gene.clusters[!duplicated(gene.clusters$hgnc_symbol), ]
    
    # Assign symbol to rownames
    rownames(gene.clusters) <- gene.clusters$hgnc_symbol
  
# Extract gene cluster numbers
gene.modules <- unique(gene.clusters$cluster)

# Loop over each cluster and run ORA
for (ith_cluster in gene.modules) {
  
    # Define the ith genelist
    ith_gene_list <- rownames(gene.clusters[which(gene.clusters$cluster == ith_cluster), ])
      
    # using GO BP Collection
    ORA.GO <- Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = universe.genes, my.genelist = ith_gene_list, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = paste('Cluster_', ith_cluster, sep = ''))
    
    # using Hallmark
    ORA.Hallmark <- Run_ORA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', universe.list = universe.genes, my.genelist = ith_gene_list, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = paste('Cluster_', ith_cluster, sep = ''))

    # using Reactome
    ORA.Reactome <- Run_ORA(output.dir = my.output.dir, output_subfolder = 'Reactome/', universe.list = universe.genes, my.genelist = ith_gene_list, my.gs.collection = Reactome_Geneset, gs.label = 'Reactome', condition_label = paste('Cluster_', ith_cluster, sep = ''))

}



# Set 4: T72 WI38 ------------------------------------------------------------------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_T72_WI38/'

# Make gene lists
genelist.WI38 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_WI38_DEGs_All.txt")



# Run GSEA using GO BP Collection
GOBP.WI38 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.WI38, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'T72_WI38', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.WI38 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.WI38, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'T72_WI38', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.WI38 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.WI38, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'T72_WI38', output.object = 'filtered')

# Run GSEA using Repeat Class Collection
Class.WI38 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.WI38, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'T72_WI38', output.object = 'unfiltered')

# Run GSEA using Repeat Family Collection
Family.WI38 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.WI38, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'T72_WI38', output.object = 'unfiltered')

    # Generate GSEA Plots for Alu family
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_T72_WI38", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.WI38, geneSetID = c("Alu subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
# Run GSEA using Aging Fibroblast and senescence genesets
Aging_Sen.WI38 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Aging_Senescence/', my.genelist = genelist.WI38, my.gs.collection = Aging_Senescence.gs, gs.label = 'Aging_Senescence_Lists', condition_label = 'T72_WI38', output.object = 'unfiltered', my.max_gs_size = Inf, items_to_plot = 10, plot_non_sig = TRUE)

# Save results
save(GOBP.WI38, Reactome.WI38, Hallmark.WI38, Class.WI38, Family.WI38, Aging_Sen.WI38,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_T72_WI38.R", sep = ''))



# Set 5: Aging Fibroblasts vs Alu OE genesets ------------------------------------------------------------------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_Aging_Fibroblasts_vs_Alu/'

# Make gene lists
genelist.aging <- make_genelist_Ensemble(DESeq_results_path = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/DESeq_Results/TETranscripts_GSE113957_Linear_Across_Age_All_Genes.txt")

# Run GSEA using Alu OE gene sets
AluOE_gs.aging <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Alu_OE_genesets/', my.genelist = genelist.aging, my.gs.collection = AluOE.gs, gs.label = 'Alu_OE', condition_label = 'Aging_Fibroblasts', output.object = 'filtered', plot_non_sig = TRUE)

# Save results
save(AluOE_gs.aging,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Aging_Fibroblasts.R", sep = ''))





# Set 6: Across Endogenous AluJb Expression in Aging Fibroblasts ------------------------------------------------------------------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Results_Endogenous_AluJb_Expression/'

# Make gene lists
genelist.endoAlu <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/DESeq_Results/TETranscripts_GSE113957_Linear_Across_AluJb_Expression_All_Genes.txt")

# Run GSEA using GO BP Collection
GOBP.endoAlu <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.endoAlu, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'endo_Alu', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.endoAlu <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.endoAlu, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'endo_Alu', output.object = 'filtered')

# Save results
save(GOBP.endoAlu, Reactome.endoAlu,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Endogenous_AluJb_Expression_in_Aging_Fibroblasts.R", sep = ''))


# Save Session Info  ------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_vs_Transcriptomics.txt", sep =""))
sessionInfo()
sink()      
    
    
# Clean the environment
rm(list=ls())