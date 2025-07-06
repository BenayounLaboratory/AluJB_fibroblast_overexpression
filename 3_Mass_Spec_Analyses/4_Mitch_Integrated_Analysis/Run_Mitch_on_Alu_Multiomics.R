# Integrate data from transcriptome, cell proteome, and secreted proteome

# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library("mitch")
library('openxlsx')

# Set the working directory
setwd('/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/4_Mitch_Integrated_Analysis/Results/')


# Load transcriptomic and proteomics results ################################################################################


# Load transcriptome
my.transcriptome <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_All.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
        
        # Assign EnsembleID to entries with an NA or Empty cell for HGNC Symbol
        my.transcriptome[is.na(my.transcriptome$hgnc_symbol), 'hgnc_symbol'] <- rownames(my.transcriptome[is.na(my.transcriptome$hgnc_symbol), ])
        my.transcriptome[which(my.transcriptome$hgnc_symbol == ''), 'hgnc_symbol'] <- rownames(my.transcriptome[which(my.transcriptome$hgnc_symbol == ''), ])
        
        # Remove duplicate gene names
        my.transcriptome <- my.transcriptome[!duplicated(my.transcriptome$hgnc_symbol), ]
        
        # Assign symbol to rownames
        rownames(my.transcriptome) <- my.transcriptome$hgnc_symbol
        
        # only keep necessary columns
        my.transcriptome <- my.transcriptome[, c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")]
        
# Load proteome
my.cell_proteins <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Cellular_Proteome_Differential_Abundance_Results.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

        # Split concatenated genes
        my.cell_proteins <- as.data.frame(splitstackshape::cSplit(my.cell_proteins, splitCols = "Genes", sep = ";", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE)) # warning about specifying 'as.is' by the caller can be ignored
    
        # Remove duplicate gene names
        my.cell_proteins <- my.cell_proteins[!duplicated(my.cell_proteins$Genes), ]
        
        # Assign gene names to the list 
        rownames(my.cell_proteins) <- my.cell_proteins$Genes
        
        # only keep necessary columns
        my.cell_proteins <- my.cell_proteins[, c("AVG.Group.Quantity.Denominator", "AVG.Log2.Ratio","Standard.Error","Rank_metric","Pvalue","Qvalue")]
        
        # reformat colnames to match DESeq2
        colnames(my.cell_proteins) <- c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")
     
# Load secretome           
my.secreted_proteins <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Secreted_Proteome_Differential_Abundance_Results.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

        # Split concatenated genes
        my.secreted_proteins <- as.data.frame(splitstackshape::cSplit(my.secreted_proteins, splitCols = "Genes", sep = ";", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE)) # warning about specifying 'as.is' by the caller can be ignored
    
        # Remove duplicate gene names
        my.secreted_proteins <- my.secreted_proteins[!duplicated(my.secreted_proteins$Genes), ]
        
        # Assign gene names to the list 
        rownames(my.secreted_proteins) <- my.secreted_proteins$Genes
        
        # only keep necessary columns
        my.secreted_proteins <- my.secreted_proteins[, c("AVG.Group.Quantity.Denominator", "AVG.Log2.Ratio","Standard.Error","Rank_metric","Pvalue","Qvalue")]
        
        # reformat colnames to match DESeq2
        colnames(my.secreted_proteins) <- c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")
    

# Load gene sets 
load(file = '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Gene_Set_Collections_for_GSEA.R')

    # Convert gene sets to GMT format
    Reactome_Geneset <- split(x = Reactome_Geneset$gene, f = Reactome_Geneset$gs_name)
    GO_BP_Geneset <- split(x = GO_BP_Geneset$gene, f = GO_BP_Geneset$gs_name)
    Hallmark_Geneset <- split(x = Hallmark_Geneset$gene, f = Hallmark_Geneset$gs_name)
    Aging_Senescence.gs <- split(x = Aging_Senescence.gs$gene, f = Aging_Senescence.gs$gs_name)



# Import data for mitch and run analysis ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Create structure for the mitch data import
cell.types.mitch <- c('IMR90')
mitch.data <- vector(mode = "list", length = length(cell.types.mitch) )
names(mitch.data) <- 'IMR90'

# create structure for the mitch results
mitch.res.gobp  <- vector(mode = "list", length = length(cell.types.mitch) )
mitch.res.reactome <- vector(mode = "list", length = length(cell.types.mitch) )
mitch.res.Hallmark <- vector(mode = "list", length = length(cell.types.mitch) )
mitch.res.Age_Sen <- vector(mode = "list", length = length(cell.types.mitch) )

names(mitch.res.gobp) <- cell.types.mitch
names(mitch.res.reactome) <- cell.types.mitch
names(mitch.res.Hallmark) <- cell.types.mitch
names(mitch.res.Age_Sen) <- cell.types.mitch


# Run Mitch over each cell type
for (i in 1:length(cell.types.mitch) ) {
 
  # create required list
  tmp.mitch.list <- list("RNA"  = my.transcriptome,
                         "Cell_Proteins" = my.cell_proteins,
                         "Secreted_Proteins" = my.secreted_proteins)
  
  
  # generate mitch input
  mitch.data[[i]] <- mitch_import(tmp.mitch.list, DEtype = "DESeq2")

  # prioritisation by significance
  mitch.res.gobp [[i]] <- mitch_calc(mitch.data[[i]], GO_BP_Geneset, priority = "significance", minsetsize = 10, resrows = 25, cores = 7)
  mitch.res.reactome[[i]] <- mitch_calc(mitch.data[[i]], Reactome_Geneset, priority = "significance", minsetsize = 10, resrows = 25, cores = 7)
  mitch.res.Hallmark[[i]] <- mitch_calc(mitch.data[[i]], Hallmark_Geneset, priority = "significance", minsetsize = 10, resrows = 9, cores = 7) 
  mitch.res.Age_Sen[[i]] <- mitch_calc(mitch.data[[i]], Aging_Senescence.gs, priority = "significance", minsetsize = 10, resrows = 25, cores = 7)

  # Generate output files
  mitch_report(mitch.res.gobp[[i]] , paste(Sys.Date(), cell.types.mitch[i],"GOBP_mitch_report.html", sep = "_"))
  mitch_plots(mitch.res.gobp[[i]]  , outfile = paste(Sys.Date(),cell.types.mitch[i],"GOBP_mitch_charts.pdf", sep = "_"))

  mitch_report(mitch.res.reactome[[i]] , paste(Sys.Date(),cell.types.mitch[i],"REACTOME_mitch_report.html", sep = "_"))
  mitch_plots(mitch.res.reactome[[i]]  , outfile = paste(Sys.Date(),cell.types.mitch[i],"REACTOME_mitch_charts.pdf", sep = "_"))
  
  mitch_report(mitch.res.Hallmark[[i]] , paste(Sys.Date(),cell.types.mitch[i],"Hallmark_mitch_report.html", sep = "_"))
  mitch_plots(mitch.res.Hallmark[[i]]  , outfile = paste(Sys.Date(),cell.types.mitch[i],"Hallmark_mitch_charts.pdf", sep = "_"))
  
  mitch_report(mitch.res.Age_Sen[[i]] , paste(Sys.Date(), cell.types.mitch[i], "Aging_Senescence_mitch_report.html", sep = "_"))
  mitch_plots(mitch.res.Age_Sen[[i]]  , outfile = paste(Sys.Date(), cell.types.mitch[i], "Aging_Senescence_mitch_charts.pdf", sep = "_"))

}

save(mitch.data, 
     mitch.res.gobp, mitch.res.reactome, mitch.res.Hallmark, mitch.res.Age_Sen,
     file = paste0(Sys.Date(),"_IMR90_AluOE_RNA_Cellular_and_Secreted_Proteins_mitch_Results.RData"))



# Export to xlsx ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


options(java.parameters = "-Xmx16g" )
require(openxlsx)

# Define function to get enrichments
get_enrich_res <- function (mitch.res) {
  mitch.res$enrichment_result
}

# Output results
GOBP.out <-  paste0(Sys.Date(),"_mitch_Results_GOBP_IMR90_AluOE_RNA_Cellular_and_Secreted_Proteins.xlsx")
write.xlsx(lapply(mitch.res.gobp, get_enrich_res), rowNames = F, file = GOBP.out)


REACT.out <-  paste0(Sys.Date(),"_mitch_Results_REACTOME_IMR90_AluOE_RNA_Cellular_and_Secreted_Proteins.xlsx")
write.xlsx(lapply(mitch.res.reactome, get_enrich_res), rowNames = F, file = REACT.out)


Hallmark.out <-  paste0(Sys.Date(),"_mitch_Results_Hallmark_IMR90_AluOE_RNA_Cellular_and_Secreted_Proteins.xlsx")
write.xlsx(lapply(mitch.res.Hallmark, get_enrich_res), rowNames = F, file = Hallmark.out)

Aging_Sen.out <-  paste0(Sys.Date(),"_mitch_Results_Aging_Senescence_IMR90_AluOE_RNA_Cellular_and_Secreted_Proteins.xlsx")
write.xlsx(lapply(mitch.res.Age_Sen, get_enrich_res), rowNames = F, file = Aging_Sen.out)



# SESSION INFO ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/4_Mitch_Integrated_Analysis/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_Mitch_Transcriptomics_and_Cellular_Secreted_Proteomics.txt", sep =""))
sessionInfo()
sink()      
    
    
# Clean the environment
rm(list=ls())      
