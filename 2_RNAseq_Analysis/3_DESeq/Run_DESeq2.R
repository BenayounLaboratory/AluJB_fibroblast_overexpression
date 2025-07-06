# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Process_RNASeq_functions.R")

# Load libraries
library(DESeq2) # For differential expression
library(BiocParallel) # Used by DESeq2 for parallelization
  register(MulticoreParam(6)) # Use six cores
library(sva) # Correcting for batch effects
library(limma) # Correcting for batch effects
library(biomaRt) # To map Ensembl IDs to gene symbols
library(DEGreport)
library(beeswarm)

# Define the output directories
counts.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/'
DESeq.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/'
MDS.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/MDS_PCA/'
    
    

# Set 1: T24 IMR-90 AluJb OE ------------------------------------------------------------------------------------------------------------------------


# CALCULATE FRACTION OF READS FROM REPEATS 

# Load counts
gene_TE_counts <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/Filtered_counts_T24_AluJb_OE.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Define repeat-only counts
repeat_counts <- gene_TE_counts[!grepl('ENSG', rownames(gene_TE_counts)), ]

# Sum repeat-only and total counts
repeat_sum <- colSums(repeat_counts)
total_sum <- colSums(gene_TE_counts)

# Calculate fraction of reads from repeats
fraction_repeats <- repeat_sum/total_sum

    # Multiply by 100 to get percents
    fraction_repeats <- fraction_repeats * 100

# Define group names
treatment.group <- c(rep('Empty', 4), rep('AluJb', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(gene_TE_counts),
                         Treatment = treatment.group,
                         Fraction = fraction_repeats)

# define factors for boxplots
SampleInfo$Treatment <- factor(SampleInfo$Treatment, levels = c("Empty", "AluJb"))

# get p-values
Alu.stats <- wilcox.test(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'Fraction'],
                         SampleInfo[which(SampleInfo$Treatment == 'AluJb'), 'Fraction'])

# Make a list of the data for this metric
fraction.plot <- list("Empty" = SampleInfo$Fraction[SampleInfo$Treatment == "Empty"],
                      "AluJb" = SampleInfo$Fraction[SampleInfo$Treatment == "AluJb"]
                      )

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'Fraction'])

# Generate the plot
pdf(paste(counts.dir, "Boxplot_T24_AluJb_OE_Repeat_Read_Percents.pdf", sep = ''), height = 5, width = 5)

    boxplot(fraction.plot, 
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(0, 5),
            ylab    = "% of repeat-mapping reads",
            outline = F ,
            main    = "Repeat Abundance")
    
    beeswarm(fraction.plot, add = T, pch = 16)
    
    # Add pvalues
    text(1.5,  5, paste('p = ', signif(Alu.stats$p.value,3), sep = ''))

    # Add line at y = 1
    abline(h = Empty.median, lty = 2)

dev.off()




    
# DIFFERENTIAL EXPRESSION ANALYSIS (DESEQ2)
    
# Load counts
gene_TE_counts <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/Filtered_counts_T24_AluJb_OE.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# General parameters
padj_limit <- 0.05 # Define alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

# Define covariates
treatment.group <- c(rep('Empty', 4), rep('AluJb', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(gene_TE_counts),
                         Treatment = treatment.group)

# Create DESeq2 object
# NOTE: the warning message about converting variables from characters to factors does not affect results and can be ignored
dds <- DESeqDataSetFromMatrix(countData = gene_TE_counts,
                              colData = SampleInfo,
                              design = ~ Treatment) 

# run DESeq2
dds <- DESeq(dds, parallel = TRUE)

# Get VST data (includes size factor normalization)
normalized_counts <- assay(varianceStabilizingTransformation(dds, blind = FALSE))

     # Save VST data
     write.table(normalized_counts, file = paste(counts.dir, 'VST_Expression_T24_AluJb', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

#Extract DESeq results. 
AluOE_vs_Control <- results(dds, contrast = c("Treatment", "AluJb", "Empty"), alpha = padj_limit, independentFiltering = TRUE)

# DESeq Stats
summary(AluOE_vs_Control, alpha = padj_limit)

# Extract significant results and save.
AluOE_vs_Control.sig <- Extract_DESeq_stats(DESeq_results = AluOE_vs_Control, padj_limit = padj_limit, organism = 'hs', 
                                            output.dir = DESeq.dir,
                                            output_file_prefix_all = 'T24_AluJb_DEGs_All',
                                            output_file_prefix_sig = 'T24_AluJb_DEGs_FDR5')
  
  
  
  
  
# MDS/PCA ANALYSIS

# Define point colors and shapes
my.colors <- c(rep('dodgerblue', 4), rep('dodgerblue4', 4))
my.shapes <- c(16, 16, 16, 16, 16, 16, 16, 16)

# Run MDS 
Run_MDS(VST_expression = normalized_counts, 
        treatment_label = 'T24_AluJb', 
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = MDS.dir)

# Run PCA 
Run_PCA(VST_expression = normalized_counts, 
        treatment_label = 'T24_AluJb', 
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = MDS.dir)
    
    

    


# Set 2: T72 IMR-90 Proteomics ------------------------------------------------------------------------------------------------------------------------
    
    
# CALCULATE FRACTION OF READS FROM REPEATS 

# Load counts
gene_TE_counts <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/Filtered_counts_T72_Proteomics.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Define repeat-only counts
repeat_counts <- gene_TE_counts[!grepl('ENSG', rownames(gene_TE_counts)), ]

# Sum repeat-only and total counts
repeat_sum <- colSums(repeat_counts)
total_sum <- colSums(gene_TE_counts)

# Calculate fraction of reads from repeats
fraction_repeats <- repeat_sum/total_sum

    # Multiply by 100 to get percents
    fraction_repeats <- fraction_repeats * 100

# Define group names
treatment.group <- c(rep('Empty', 4), rep('AluJb', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(gene_TE_counts),
                         Treatment = treatment.group,
                         Fraction = fraction_repeats)

# define factors for boxplots
SampleInfo$Treatment <- factor(SampleInfo$Treatment, levels = c("Empty", "AluJb"))

# get p-values
Alu.stats <- wilcox.test(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'Fraction'],
                         SampleInfo[which(SampleInfo$Treatment == 'AluJb'), 'Fraction'])

# Make a list of the data for this metric
fraction.plot <- list("Empty" = SampleInfo$Fraction[SampleInfo$Treatment == "Empty"],
                      "AluJb" = SampleInfo$Fraction[SampleInfo$Treatment == "AluJb"]
                      )

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'Fraction'])

# Generate the plot
pdf(paste(counts.dir, "Boxplot_T72_Proteomics_Samples_Repeat_Read_Percents.pdf", sep = ''), height = 5, width = 5)

    boxplot(fraction.plot, 
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(0, 5),
            ylab    = "% of repeat-mapping reads",
            outline = F ,
            main    = "Repeat Abundance")
    
    beeswarm(fraction.plot, add = T, pch = 16)
    
    # Add pvalues
    text(1.5,  5, paste('p = ', signif(Alu.stats$p.value,3), sep = ''))

    # Add line at y = 1
    abline(h = Empty.median, lty = 2)

dev.off()





# DIFFERENTIAL EXPRESSION ANALYSIS (DESEQ2)
    
# Load counts
gene_TE_counts <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/Filtered_counts_T72_Proteomics.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# General parameters
padj_limit <- 0.05 # Define alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

# Define covariates
treatment.group <- c(rep('Empty', 4), rep('AluJb', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(gene_TE_counts),
                         Treatment = treatment.group)

# Create DESeq2 object
# NOTE: the warning message about converting variables from characters to factors does not affect results and can be ignored
dds <- DESeqDataSetFromMatrix(countData = gene_TE_counts,
                              colData = SampleInfo,
                              design = ~ Treatment) 

# run DESeq2
dds <- DESeq(dds, parallel = TRUE)

# Get VST data (includes size factor normalization)
normalized_counts <- assay(varianceStabilizingTransformation(dds, blind = FALSE))

     # Save VST data
     write.table(normalized_counts, file = paste(counts.dir, 'VST_Expression_T72_Proteomics', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

#Extract DESeq results. 
AluOE_vs_Control <- results(dds, contrast = c("Treatment", "AluJb", "Empty"), alpha = padj_limit, independentFiltering = TRUE)

# DESeq Stats
summary(AluOE_vs_Control, alpha = padj_limit)

# Extract significant results and save.
AluOE_vs_Control.sig <- Extract_DESeq_stats(DESeq_results = AluOE_vs_Control, padj_limit = padj_limit, organism = 'hs', 
                                            output.dir = DESeq.dir,
                                            output_file_prefix_all = 'T72_Proteomics_DEGs_All',
                                            output_file_prefix_sig = 'T72_Proteomics_DEGs_FDR5')
  
  
  
  
  
# MDS/PCA ANALYSIS

# Define point colors and shapes
my.colors <- c(rep('dodgerblue', 4), rep('dodgerblue4', 4))
my.shapes <- c(16, 16, 16, 16, 16, 16, 16, 16)

# Run MDS 
Run_MDS(VST_expression = normalized_counts, 
        treatment_label = 'T72_Proteomics', 
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = MDS.dir)

# Run PCA 
Run_PCA(VST_expression = normalized_counts, 
        treatment_label = 'T72_Proteomics', 
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = MDS.dir)
    
    

    

# Set 3: T72 IMR-90 Lamivudine ------------------------------------------------------------------------------------------------------------------------
    
    
# SURROGATE VARIABLE ANALYSIS (SVA)

# Load counts
filtered_counts <- as.matrix(read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/Filtered_counts_T72_Lamivudine.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1))

# Define treatments and covariates (!!!!!! MAKE SURE VALUES MATCH IF COLUMN ORDER IS CHANGED !!!!!!)
treatment.group <- c(rep('Empty_Veh', 4), 
                     rep('Empty_3TC', 3),
                     rep('AluJb_Veh', 4), 
                     rep('AluJb_3TC', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(filtered_counts),
                         Treatment = treatment.group
                         )

# Set the full models (adjustment + variables of interest).
mod_full = model.matrix(~ Treatment, data = SampleInfo)

# Set null models (adjustment variables)
mod_null = model.matrix(~ 1, data = SampleInfo)

# estimate the # of SVs
# n.sv.be = num.sv(dat = filtered_counts, mod = mod_full, method ="be", B = 50, seed = 12345) # 2 detected

# apply SVAseq algorithm
my_svaseq = svaseq(dat = filtered_counts, mod = mod_full, mod0 = mod_null, n.sv = 2, constant = 1, B = 50)

# Remove batch effects with Limma
my.corrected.data <- removeBatchEffect(x = log2(filtered_counts + 1.0),
                                    batch = NULL,
                                    covariates = cbind(my_svaseq$sv), # !!!!!!! Change with SVs
                                    design = mod_full)

# Delog and round data.
gene_TE_counts <- ceiling(2^my.corrected.data-1.0)

    # Save batch corrected counts (to visualize downstream results)
    write.table(gene_TE_counts, file = paste(counts.dir, 'Filtered_counts_T72_Lamivudine_SVA', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)  
  





    
# DIFFERENTIAL EXPRESSION ANALYSIS (DESEQ2)

# General parameters
padj_limit <- 0.05 # Define alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

# Define covariates
treatment.group <- c(rep('U6_Veh', 4), rep('U6_3TC', 3), rep('Alu_Veh', 4), rep('Alu_3TC', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(gene_TE_counts),
                         Treatment = treatment.group)

# Create DESeq2 object
# NOTE: the warning message about converting variables from characters to factors does not affect results and can be ignored
dds <- DESeqDataSetFromMatrix(countData = gene_TE_counts,
                              colData = SampleInfo,
                              design = ~ Treatment) 

# run DESeq2 (for WALD TEST)
dds <- DESeq(dds, parallel = TRUE)

# Get VST data (includes size factor normalization)
normalized_counts <- assay(varianceStabilizingTransformation(dds, blind = FALSE))

     # Save VST data
     write.table(normalized_counts, file = paste(counts.dir, 'VST_Expression_T72_Lamivudine', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

#Extract DESeq results. 
Alu_vs_Empty.Veh <- results(dds, contrast = c("Treatment", "Alu_Veh", "U6_Veh"), alpha = padj_limit, independentFiltering = TRUE)     
Alu_vs_Empty.Lam <- results(dds, contrast = c("Treatment", "Alu_3TC", "U6_3TC"), alpha = padj_limit, independentFiltering = TRUE)     
Lam_vs_Veh.Empty <- results(dds, contrast = c("Treatment", "U6_3TC", "U6_Veh"), alpha = padj_limit, independentFiltering = TRUE)     
Lam_vs_Veh.Alu <- results(dds, contrast = c("Treatment", "Alu_3TC", "Alu_Veh"), alpha = padj_limit, independentFiltering = TRUE)     
     
# DESeq Stats
summary(Alu_vs_Empty.Veh, alpha = padj_limit)
summary(Alu_vs_Empty.Lam, alpha = padj_limit)
summary(Lam_vs_Veh.Empty, alpha = padj_limit)
summary(Lam_vs_Veh.Alu, alpha = padj_limit)

# Extract significant results and save.
Alu_vs_Empty.Veh <- Extract_DESeq_stats(DESeq_results = Alu_vs_Empty.Veh, padj_limit = padj_limit, organism = 'hs', 
                                        output.dir = DESeq.dir,
                                        output_file_prefix_all = 'T72_LamExperiment_Alu_vs_Empty_Veh_DEGs_All',
                                        output_file_prefix_sig = 'T72_LamExperiment_Alu_vs_Empty_Veh_DEGs_FDR5')
  
Alu_vs_Empty.Lam <- Extract_DESeq_stats(DESeq_results = Alu_vs_Empty.Lam, padj_limit = padj_limit, organism = 'hs', 
                                        output.dir = DESeq.dir,
                                        output_file_prefix_all = 'T72_LamExperiment_Alu_vs_Empty_3TC_DEGs_All',
                                        output_file_prefix_sig = 'T72_LamExperiment_Alu_vs_Empty_3TC_DEGs_FDR5')

Lam_vs_Veh.Empty <- Extract_DESeq_stats(DESeq_results = Lam_vs_Veh.Empty, padj_limit = padj_limit, organism = 'hs', 
                                        output.dir = DESeq.dir,
                                        output_file_prefix_all = 'T72_LamExperiment_3TC_vs_Veh_Empty_DEGs_All',
                                        output_file_prefix_sig = 'T72_LamExperiment_3TC_vs_Veh_Empty_DEGs_FDR5')

Lam_vs_Veh.Alu <- Extract_DESeq_stats(DESeq_results = Lam_vs_Veh.Alu, padj_limit = padj_limit, organism = 'hs', 
                                        output.dir = DESeq.dir,
                                        output_file_prefix_all = 'T72_LamExperiment_3TC_vs_Veh_Alu_DEGs_All',
                                        output_file_prefix_sig = 'T72_LamExperiment_3TC_vs_Veh_Alu_DEGs_FDR5')
  



# run DESeq2 (Likelihood Ratio Test)
dds_lrt <- DESeq(dds, test = 'LRT', reduced = ~1, parallel = TRUE)

# Extract results
res_LRT <- results(dds_lrt)

# DESeq Stats
summary(res_LRT, alpha = padj_limit)

# Extract significant results and save.
LRT.All_Groups.sig <- Extract_DESeq_stats(DESeq_results = res_LRT, padj_limit = padj_limit, organism = 'hs', 
                                          output.dir = DESeq.dir,
                                          output_file_prefix_all = 'LRT_LamExperiment_All_Groups_All_Genes',
                                          output_file_prefix_sig = 'LRT_LamExperiment_All_Groups_FDR5')

# Define significant LRT gene names
LRT.sig.names <- rownames(LRT.All_Groups.sig)

# Extract VST data for sig genes
vst.LRT.sig <- normalized_counts[LRT.sig.names, ]   
      
# Define sample group names (note: these will be appear in alphabetical/numerical order on the plot)
group.labels <- c(rep('1_U6_Veh', 4), rep('3_U6_3TC', 3), rep('2_Alu_Veh', 4), rep('4_Alu_3TC', 4))

# Format sample group labels as a df
metadata.Alu <- data.frame(row.names = colnames(vst.LRT.sig),
                           Treatment = group.labels)

    # Update the label columns as a factor
    metadata.Alu$Treatment <- as.factor(metadata.Alu$Treatment)
    
# Run degPatterns
pdf(paste(DESeq.dir, "Plot_DEGPatterns_Lamivudine_LRT_FDR5", ".pdf", sep=""), height = 10, width = 10)

Alu.clusters <- degPatterns(vst.LRT.sig, 
                            metadata = metadata.Alu, 
                            time = "Treatment", 
                            col = NULL,
                            minc = 15,
                            reduce = FALSE,
                            plot = TRUE,
                            consensusCluster = FALSE)

dev.off()

    # Save gene-cluster table
    write.table(Alu.clusters[["df"]], file = paste(DESeq.dir, 'LRT_LamExperiment_All_Groups_degPatterns_Clustering', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)


  
      
      
# MDS/PCA ANALYSIS

# Define point colors and shapes
my.colors <- c(rep('dodgerblue', 7), rep('dodgerblue4', 8))
my.shapes <- c(16, 16, 16, 16, 17, 17, 17, 16, 16, 16, 16, 17, 17, 17, 17)

# Run MDS 
Run_MDS(VST_expression = normalized_counts, 
        treatment_label = 'T72_Lamivudine', 
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = MDS.dir)

# Run PCA 
Run_PCA(VST_expression = normalized_counts, 
        treatment_label = 'T72_Lamivudine', 
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = MDS.dir)
    
    

    

# Set 4: T72 WI-38 AluJb OE ------------------------------------------------------------------------------------------------------------------------

    
# SURROGATE VARIABLE ANALYSIS (SVA)

# Load counts
filtered_counts <- as.matrix(read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/Filtered_counts_T72_WI38.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1))

# Define treatments (!!!!!! MAKE SURE VALUES MATCH IF COLUMN ORDER IS CHANGED !!!!!!)
treatment.group <- c(rep('Empty', 4), 
                   rep('AluJb', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(filtered_counts),
                         Treatment = treatment.group
                         )

# Set the full models (adjustment + variables of interest).
mod_full = model.matrix(~ Treatment, data = SampleInfo)

# Set null models (adjustment variables)
mod_null = model.matrix(~ 1, data = SampleInfo)

# estimate the # of SVs
n.sv.be = num.sv(dat = filtered_counts, mod = mod_full, method ="be", B = 50, seed = 12345) # 2 detected

# apply SVAseq algorithm
my_svaseq = svaseq(dat = filtered_counts, mod = mod_full, mod0 = mod_null, n.sv = n.sv.be, constant = 1, B = 50)

# Remove batch effects with Limma
my.corrected.data <- removeBatchEffect(x = log2(filtered_counts + 1.0),
                                    batch = NULL,
                                    covariates = cbind(my_svaseq$sv), # !!!!!!! Change with SVs
                                    design = mod_full)

# Delog and round data.
gene_TE_counts <- ceiling(2^my.corrected.data-1.0)

    # Save batch corrected counts (to visualize downstream results)
    write.table(gene_TE_counts, file = paste(counts.dir, 'Filtered_counts_T72_WI38_SVA', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)  





    
# DIFFERENTIAL EXPRESSION ANALYSIS (DESEQ2)

# General parameters
padj_limit <- 0.05 # Define alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

# Define covariates
treatment.group <- c(rep('Empty', 4), rep('AluJb', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(gene_TE_counts),
                         Treatment = treatment.group)

# Create DESeq2 object
# NOTE: the warning message about converting variables from characters to factors does not affect results and can be ignored
dds <- DESeqDataSetFromMatrix(countData = gene_TE_counts,
                              colData = SampleInfo,
                              design = ~ Treatment) 

# run DESeq2
dds <- DESeq(dds, parallel = TRUE)

# Get VST data (includes size factor normalization)
normalized_counts <- assay(varianceStabilizingTransformation(dds, blind = FALSE))

     # Save VST data
     write.table(normalized_counts, file = paste(counts.dir, 'VST_Expression_T72_WI38', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

#Extract DESeq results. 
AluOE_vs_Control <- results(dds, contrast = c("Treatment", "AluJb", "Empty"), alpha = padj_limit, independentFiltering = TRUE)

# DESeq Stats
summary(AluOE_vs_Control, alpha = padj_limit)

# Extract significant results and save.
AluOE_vs_Control.sig <- Extract_DESeq_stats(DESeq_results = AluOE_vs_Control, padj_limit = padj_limit, organism = 'hs', 
                                            output.dir = DESeq.dir,
                                            output_file_prefix_all = 'T72_WI38_DEGs_All',
                                            output_file_prefix_sig = 'T72_WI38_DEGs_FDR5')
  
  
  
  
  
# MDS/PCA ANALYSIS

# Define point colors and shapes
my.colors <- c(rep('dodgerblue', 4), rep('dodgerblue4', 4))
my.shapes <- c(16, 16, 16, 16, 16, 16, 16, 16)

# Run MDS 
Run_MDS(VST_expression = normalized_counts, 
        treatment_label = 'T72_WI38', 
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = MDS.dir)

# Run PCA 
Run_PCA(VST_expression = normalized_counts, 
        treatment_label = 'T72_WI38', 
        plot_colors = my.colors,
        plot_shapes = my.shapes,
        MDS_dir = MDS.dir)
    
    

    

# Session Info ------------------------------------------------------------------------------------------------------------------------



# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_DESeq2_Transcriptomics.txt", sep =""))
sessionInfo()
sink()      
    


    
# Clean the environment
rm(list=ls())


