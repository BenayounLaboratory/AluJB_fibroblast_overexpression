# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(stringr) # for capitalization changes
library(splitstackshape)
library(mitch)

# Define output directory
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/'             

  



# LOAD CURATED GENE SETS  

# Read gmt files
Reactome_Geneset <- gmt_import("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Pathway_Gene_Set_Collection_GMTs/ReactomePathways_v92.gmt.txt")
GO_BP_Geneset <- gmt_import("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Pathway_Gene_Set_Collection_GMTs/c5.go.bp.v2024.1.Hs.symbols.gmt.txt")
Hallmark_Geneset <- gmt_import("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Pathway_Gene_Set_Collection_GMTs/h.all.v2024.1.Hs.symbols.gmt.txt")

# Convert to two column dataframe         

Reactome_Geneset <- data.frame(gs_name = rep(names(Reactome_Geneset), 
                               times = lengths(Reactome_Geneset)),
                               gene = unlist(Reactome_Geneset, use.names = FALSE),
                               row.names = NULL)

GO_BP_Geneset <- data.frame(gs_name = rep(names(GO_BP_Geneset), 
                            times = lengths(GO_BP_Geneset)),
                            gene = unlist(GO_BP_Geneset, use.names = FALSE),
                            row.names = NULL)

Hallmark_Geneset <- data.frame(gs_name = rep(names(Hallmark_Geneset), 
                               times = lengths(Hallmark_Geneset)),
                               gene = unlist(Hallmark_Geneset, use.names = FALSE),
                               row.names = NULL)





# CLEANUP GENE SET NAMES 

# remove redundant parts of labels
Hallmark_Geneset$gs_name <- gsub("HALLMARK_", "", Hallmark_Geneset$gs_name)
GO_BP_Geneset$gs_name <- gsub("GOBP_", "", GO_BP_Geneset$gs_name)

# write function to lowercase all letters, then uppercase the first
sentence_casing <- function(input_names){
  
  # lowercase all letters
  input_lowercased <- tolower(input_names)
  
  # uppercase the first letter
  input_sentence_case <- str_to_sentence(input_lowercased)
  
  # output the final string
  return(input_sentence_case)
  
}

# change cases of gene set names, so only the first letter is capitalized
Hallmark_Geneset$gs_name <- sentence_casing(input_names = Hallmark_Geneset$gs_name)
GO_BP_Geneset$gs_name <- sentence_casing(input_names = GO_BP_Geneset$gs_name)

# Update underscores to blank spaces 
Hallmark_Geneset$gs_name <- gsub("_", " ", Hallmark_Geneset$gs_name)
GO_BP_Geneset$gs_name <- gsub("_", " ", GO_BP_Geneset$gs_name)


  


# MAKE AGING AND SENESCENCE GENE SET

# Load gene lists
core.SASP <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Aging_and_Senescence_Gene_Lists/Basisty_2020_Core_SASP_CLEANED_UP_gProfiler.csv", header = TRUE, stringsAsFactors = FALSE, sep = ',')
CellAge.Up <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Aging_and_Senescence_Gene_Lists/CellAge_Build_3_Inducers_gProfiler.csv", header = TRUE, stringsAsFactors = FALSE, sep = ',')
CellAge.Down <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Aging_and_Senescence_Gene_Lists/CellAge_Build_3_Repressors_gProfiler.csv", header = TRUE, stringsAsFactors = FALSE, sep = ',')
SenMayo <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Aging_and_Senescence_Gene_Lists/Senmayo_gProfiler_hsapiens_4-3-2025.csv", header = TRUE, stringsAsFactors = FALSE, sep = ',')
GTEx.Up <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Aging_and_Senescence_Gene_Lists/Jia_2018_GTEx_Aging_Genes_Upregulated_gProfiler.csv", header = TRUE, stringsAsFactors = FALSE, sep = ',')
GTEx.Down <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Aging_and_Senescence_Gene_Lists/Jia_2018_GTEx_Aging_Genes_Downregulated_gProfiler.csv", header = TRUE, stringsAsFactors = FALSE, sep = ',')
aging.sig <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/DESeq_Results/TETranscripts_GSE113957_Linear_Across_Age_FDR5.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')

    # Split all sig genes and repeats 
    aging.sig.repeats <- aging.sig[!grepl('ENSG', rownames(aging.sig)), ]
    aging.sig.genes <- aging.sig[grepl('ENSG', rownames(aging.sig)), ]
    
    # Assign EnsembleID to entries with an NA or Empty cell for HGNC Symbol
    aging.sig.repeats[is.na(aging.sig.repeats$hgnc_symbol), 'hgnc_symbol'] <- rownames(aging.sig.repeats[is.na(aging.sig.repeats$hgnc_symbol), ])
    aging.sig.repeats[which(aging.sig.repeats$hgnc_symbol == ''), 'hgnc_symbol'] <- rownames(aging.sig.repeats[which(aging.sig.repeats$hgnc_symbol == ''), ])
    
    aging.sig.genes[is.na(aging.sig.genes$hgnc_symbol), 'hgnc_symbol'] <- rownames(aging.sig.genes[is.na(aging.sig.genes$hgnc_symbol), ])
    aging.sig.genes[which(aging.sig.genes$hgnc_symbol == ''), 'hgnc_symbol'] <- rownames(aging.sig.genes[which(aging.sig.genes$hgnc_symbol == ''), ])
    
    # Only keep sig genes/repeats split by up/down change
    aging.sig.repeats.up <- aging.sig.repeats[which(aging.sig.repeats$log2FoldChange > 0), ]
    aging.sig.repeats.down <- aging.sig.repeats[which(aging.sig.repeats$log2FoldChange < 0), ]
    aging.sig.genes.up <- aging.sig.genes[which(aging.sig.genes$log2FoldChange > 0), ]
    aging.sig.genes.down <- aging.sig.genes[which(aging.sig.genes$log2FoldChange < 0), ]
    
    # Remove variables that aren't needed
    rm(aging.sig, aging.sig.repeats, aging.sig.genes)

# Generate individual gene sets
core.SASP <- data.frame(gs_name = 'Core_SASP', gene = core.SASP$initial_alias)
CellAge.Up <- data.frame(gs_name = 'CellAge_Inducer', gene = CellAge.Up$initial_alias)
CellAge.Down <- data.frame(gs_name = 'CellAge_Repressor', gene = CellAge.Down$initial_alias)
SenMayo <- data.frame(gs_name = 'SenMayo', gene = SenMayo$initial_alias)
GTEx.Up <- data.frame(gs_name = 'GTEx_Aging_Up', gene = GTEx.Up$initial_alias)
GTEx.Down <- data.frame(gs_name = 'GTEx_Aging_Down', gene = GTEx.Down$initial_alias)
aging.sig.repeats.up <- data.frame(gs_name = 'Aging_Fibroblast_Repeats_Up', gene = aging.sig.repeats.up$hgnc_symbol)
aging.sig.repeats.down <- data.frame(gs_name = 'Aging_Fibroblast_Repeats_Down', gene = aging.sig.repeats.down$hgnc_symbol)
aging.sig.genes.up <- data.frame(gs_name = 'Aging_Fibroblast_Genes_Up', gene = aging.sig.genes.up$hgnc_symbol)
aging.sig.genes.down <- data.frame(gs_name = 'Aging_Fibroblast_Genes_Down', gene = aging.sig.genes.down$hgnc_symbol)

# Generate the gene set collection
Aging_Senescence.gs <- rbind(core.SASP, CellAge.Up, CellAge.Down, SenMayo, GTEx.Up, GTEx.Down, aging.sig.repeats.up, aging.sig.repeats.down, aging.sig.genes.up, aging.sig.genes.down)





# MAKE ALUJb OE (MULTIOMIC SAMPLES) GENE SET (KEEP ENSEMBL NAMES)
Alu.sig <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_FDR5.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')

    # Split all sig genes and repeats 
    Alu.sig.repeats <- Alu.sig[!grepl('ENSG', rownames(Alu.sig)), ]
    Alu.sig.genes <- Alu.sig[grepl('ENSG', rownames(Alu.sig)), ]
    
    # Only keep sig genes/repeats split by up/down change
    Alu.sig.repeats.up <- Alu.sig.repeats[which(Alu.sig.repeats$log2FoldChange > 0), ]
    Alu.sig.repeats.down <- Alu.sig.repeats[which(Alu.sig.repeats$log2FoldChange < 0), ]
    Alu.sig.genes.up <- Alu.sig.genes[which(Alu.sig.genes$log2FoldChange > 0), ]
    Alu.sig.genes.down <- Alu.sig.genes[which(Alu.sig.genes$log2FoldChange < 0), ]
    
    # Remove variables that aren't needed
    rm(Alu.sig, Alu.sig.repeats, Alu.sig.genes)
    
# Generate individual gene sets
Alu.sig.repeats.up <- data.frame(gs_name = 'Alu_OE_Repeats_Up', gene = rownames(Alu.sig.repeats.up))
Alu.sig.repeats.down <- data.frame(gs_name = 'Alu_OE_Repeats_Down', gene = rownames(Alu.sig.repeats.down))
Alu.sig.genes.up <- data.frame(gs_name = 'Alu_OE_Genes_Up', gene = rownames(Alu.sig.genes.up))
Alu.sig.genes.down <- data.frame(gs_name = 'Alu_OE_Genes_Down', gene = rownames(Alu.sig.genes.down))

# Generate the gene set collection
AluOE.gs <- rbind(Alu.sig.repeats.up, Alu.sig.repeats.down, Alu.sig.genes.up, Alu.sig.genes.down)





# MAKE ALUS OE (Cantarella 2019) GENE SET 

# Load AluS OE DEGs
AluSx.sig <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Cantarella_2019_AluS_OE_Gene_Lists/Cantarella_2019_AluSx_DEGs.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
AluSq.sig <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/0_Shared_Metadata_and_External_Resources/Cantarella_2019_AluS_OE_Gene_Lists/Cantarella_2019_AluSq2_DEGs.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Aggregate into 1 list
    Alu.sig <- rbind(AluSx.sig, AluSq.sig)
    
    # Remove transcript or ENSEMBL version info. Keeping this can raise issues with gene name mapping.
    Alu.sig$Ensembl.gene.ID <- sub("\\..*", "", Alu.sig$Ensembl.gene.ID, fixed=FALSE)
  
    # Separate up and down genes
    AluS.DEGs.sig.up <- Alu.sig[which(Alu.sig$log2FoldChange > 0), ]
    AluS.DEGs.sig.down <- Alu.sig[which(Alu.sig$log2FoldChange < 0), ]
    
    # Remove duplicates
    AluS.DEGs.sig.up <- AluS.DEGs.sig.up[which(!duplicated(AluS.DEGs.sig.up$Ensembl.gene.ID)), ]
    AluS.DEGs.sig.down <- AluS.DEGs.sig.down[which(!duplicated(AluS.DEGs.sig.down$Ensembl.gene.ID)), ]

    # Remove variables that aren't needed
    rm(AluSx.sig, AluSq.sig, Alu.sig)
    
# Generate individual gene sets
AluS.DEGs.sig.up <- data.frame(gs_name = 'AluS_OE_Genes_Up', gene = rownames(AluS.DEGs.sig.up))
AluS.DEGs.sig.down <- data.frame(gs_name = 'AluS_OE_Genes_Down', gene = rownames(AluS.DEGs.sig.down))

# Generate the gene set collection
Cantarella.gs <- rbind(AluS.DEGs.sig.up, AluS.DEGs.sig.down)





# MAKE TRANSPOSON GENE SETS (BY CLASS AND FAMILY, USING TETRANSCRIPTS DESIGNATIONS)

# Load TETranscripts raw counts file (which has all the annotated transposons)  
my.counts <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/2_Read_counting/counts_T72_AluJb_IMR90_for_Proteomics/fastp_Alu_2Aligned.sortedByCoord.out.bam.cntTable", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')
  
# Define rows that are TEs
TE.rows <- which(!grepl('ENSG', rownames(my.counts)) & !grepl('U6_AluJb', rownames(my.counts)))

# Make dataframe to hold TE classifications
my.TE.classifications <- data.frame(row.names = rownames(my.counts[TE.rows, , drop = FALSE]),
                                    Full_label = rownames(my.counts[TE.rows, , drop = FALSE]),
                                    To_split = rownames(my.counts[TE.rows, , drop = FALSE]))
  
# Split TE labels into subfamily, family, and class. 
# The warning message that " 'as.is' should be specified by the caller" does not affect results and can be ignored
my.TE.classifications <- as.data.frame(splitstackshape::cSplit(my.TE.classifications, splitCols = "To_split", sep = ":", direction = "wide", fixed = TRUE, drop = TRUE, stripWhite = FALSE, makeEqual = TRUE))

# Update the column names
names(my.TE.classifications)[2:4] <- c('Subfamily', 'Family', 'Class')   

# Append "subfamilies" to the Family/Class names, since these names will be used as the gene set names
my.TE.classifications$Class <- paste(my.TE.classifications$Class, ' subfamilies', sep = '')
my.TE.classifications$Family <- paste(my.TE.classifications$Family, ' subfamilies', sep = '')
  
# Make gene set (ready for GSEA)
Repeat_class.gs <- data.frame(gs_name = my.TE.classifications$Class, gene = my.TE.classifications$Full_label)
Repeat_family.gs <- data.frame(gs_name = my.TE.classifications$Family, gene = my.TE.classifications$Full_label)
           
          



# MAKE Alu-subfamily GENE SETS (BY Alu AGE, USING TETRANSCRIPTS DESIGNATIONS)
  
# Load TETranscripts raw counts file (which has all the annotated transposons)  
my.counts <- read.csv(file = "/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/2_Read_counting/counts_T72_AluJb_IMR90_for_Proteomics/fastp_Alu_2Aligned.sortedByCoord.out.bam.cntTable", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')
  
# Define rows that are L1
L1.rows <- which(grepl(':Alu:', rownames(my.counts)))

# Make dataframe to hold TE classifications
my.Alu.classifications <- data.frame(row.names = rownames(my.counts[L1.rows, , drop = FALSE]),
                                    Full_label = rownames(my.counts[L1.rows, , drop = FALSE]),
                                    To_split = rownames(my.counts[L1.rows, , drop = FALSE]))


# Segregate by age
AluJ_subfamilies <- my.Alu.classifications[which(grepl('AluJ', rownames(my.Alu.classifications))), ]
AluS_subfamilies <- my.Alu.classifications[which(grepl('AluS', rownames(my.Alu.classifications))), ]
AluY_subfamilies <- my.Alu.classifications[which(grepl('AluY', rownames(my.Alu.classifications))), ]

# Add age columns
AluJ_subfamilies$Age <- 'AluJ subfamilies'
AluS_subfamilies$Age <- 'AluS subfamilies'
AluY_subfamilies$Age <- 'AluY subfamilies'

# Rejoin L1 subfamilies
my.Alu.classifications <- rbind(AluJ_subfamilies, AluS_subfamilies, AluY_subfamilies)
  
# Make gene set (ready for GSEA)
Alu_by_age.gs <- data.frame(gs_name = my.Alu.classifications$Age, gene = my.Alu.classifications$Full_label)





# Save gene sets
save(Hallmark_Geneset,
     Reactome_Geneset,
     GO_BP_Geneset,
     Aging_Senescence.gs,
     AluOE.gs,
     Cantarella.gs,
     Repeat_class.gs, 
     Repeat_family.gs,
     Alu_by_age.gs,
     file = paste(dir.output, "Gene_Set_Collections_for_GSEA.R", sep = ''))
           
            



# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/5_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_Prepare_Gene_Sets.txt", sep =""))
sessionInfo()
sink()   

      
      
# Clean the environment
rm(list=ls())      

