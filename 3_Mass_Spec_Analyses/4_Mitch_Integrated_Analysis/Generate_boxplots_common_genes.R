# Generate expression boxplots for genes/proteins differentially regulated across -omes

# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(beeswarm)

# Set the working directory
output.dir <- c('/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/4_Mitch_Integrated_Analysis/Common_Gene_Boxplots/')


# BOXPLOTS FOR TRANSCRIPTOME ----------------------------------------------------------------------------------------------------------------------------------------------------------------


# Load expressions values and differential analysis results
expr.VST <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/Processed_counts/VST_Expression_T72_Proteomics.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
diff.expr <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_FDR5.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Define common genes/proteins
common_genes <- c('DCN', 'NES', 'POSTN', 'DBI')

# Subset differential results for common genes
diff.expr.common <- diff.expr[which(diff.expr$hgnc_symbol %in% common_genes), ]

# Subset expression for common genes
expr.VST.common <- expr.VST[rownames(diff.expr.common), ]
rownames(expr.VST.common) <- diff.expr.common$hgnc_symbol

# Define group names
treatment.group <- c(rep('Empty', 4), rep('AluJb', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(expr.VST.common),
                         Treatment = treatment.group,
                         DCN = t(expr.VST.common['DCN', ]),
                         NES = t(expr.VST.common['NES', ]),
                         POSTN = t(expr.VST.common['POSTN', ]),
                         DBI = t(expr.VST.common['DBI', ]))

# define factors for boxplots
SampleInfo$Treatment <- factor(SampleInfo$Treatment, levels = c("Empty", "AluJb"))





# DCN 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'DCN'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Transcriptome_DCN.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(DCN ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log2(counts)",
            outline = F ,
            main    = "DCN")
    
    beeswarm(DCN ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()





# NES 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'NES'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Transcriptome_NES.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(NES ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log2(counts)",
            outline = F ,
            main    = "NES")
    
    beeswarm(NES ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()





# POSTN 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'POSTN'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Transcriptome_POSTN.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(POSTN ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log2(counts)",
            outline = F ,
            main    = "POSTN")
    
    beeswarm(POSTN ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()





# DBI 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'DBI'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Transcriptome_DBI.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(DBI ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log2(counts)",
            outline = F ,
            main    = "DBI")
    
    beeswarm(DBI ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()



# BOXPLOTS FOR CELL PROTEOME ----------------------------------------------------------------------------------------------------------------------------------------------------------------


# Load protein abundance values, remove duplicate genes, and assign gene names to rownames
expr.log <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Cellular_Proteome_Abundances_log10.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Define common genes/proteins
common_genes <- c('DCN', 'NES', 'POSTN', 'DBI')

# Subset expression for common genes
expr.cellular.common <- expr.log[which(expr.log$Genes %in% common_genes), ]
rownames(expr.cellular.common) <- expr.cellular.common$Genes

    # Remove unneeded columns
    expr.cellular.common <- expr.cellular.common[, -c(1:2)]

# Define group names
treatment.group <- c(rep('Empty', 5), rep('AluJb', 5))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(expr.cellular.common),
                         Treatment = treatment.group,
                         DCN = t(expr.cellular.common['DCN', ]),
                         NES = t(expr.cellular.common['NES', ]),
                         POSTN = t(expr.cellular.common['POSTN', ]),
                         DBI = t(expr.cellular.common['DBI', ]))

# define factors for boxplots
SampleInfo$Treatment <- factor(SampleInfo$Treatment, levels = c("Empty", "AluJb"))





# DCN 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'DCN'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Proteome_DCN.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(DCN ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log10(counts)",
            outline = F ,
            main    = "DCN")
    
    beeswarm(DCN ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()





# NES 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'NES'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Proteome_NES.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(NES ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(6.6, 7.7),
            ylab    = "Normalized log10(counts)",
            outline = F ,
            main    = "NES")
    
    beeswarm(NES ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()





# POSTN 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'POSTN'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Proteome_POSTN.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(POSTN ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(4.90, 5.25),
            ylab    = "Normalized log10(counts)",
            outline = F ,
            main    = "POSTN")
    
    beeswarm(POSTN ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()





# DBI 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'DBI'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Proteome_DBI.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(DBI ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log10(counts)",
            outline = F ,
            main    = "DBI")
    
    beeswarm(DBI ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()



# BOXPLOTS FOR SECRETOME ----------------------------------------------------------------------------------------------------------------------------------------------------------------


# Load protein abundance values, remove duplicate genes, and assign gene names to rownames
expr.log <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/Secreted_Proteome_Abundances_log10.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Define common genes/proteins
common_genes <- c('DCN', 'NES', 'POSTN', 'DBI')

# Subset expression for common genes
expr.secreted.common <- expr.log[which(expr.log$PG.Genes %in% common_genes), ]
rownames(expr.secreted.common) <- expr.secreted.common$PG.Genes

    # Remove unneeded columns
    expr.secreted.common <- expr.secreted.common[, -c(1:7)]

# Define group names
treatment.group <- c(rep('Empty', 4), rep('AluJb', 4))

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(expr.secreted.common),
                         Treatment = treatment.group,
                         DCN = t(expr.secreted.common['DCN', ]),
                         NES = t(expr.secreted.common['NES', ]),
                         POSTN = t(expr.secreted.common['POSTN', ]),
                         DBI = t(expr.secreted.common['DBI', ]))

# define factors for boxplots
SampleInfo$Treatment <- factor(SampleInfo$Treatment, levels = c("Empty", "AluJb"))





# DCN 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'DCN'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Secretome_DCN.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(DCN ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log10(counts)",
            outline = F ,
            main    = "DCN")
    
    beeswarm(DCN ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()





# NES 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'NES'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Secretome_NES.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(NES ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log10(counts)",
            outline = F ,
            main    = "NES")
    
    beeswarm(NES ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()





# POSTN 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'POSTN'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Secretome_POSTN.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(POSTN ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log10(counts)",
            outline = F ,
            main    = "POSTN")
    
    beeswarm(POSTN ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()





# DBI 

# Get median for control group
Empty.median <- median(SampleInfo[which(SampleInfo$Treatment == 'Empty'), 'DBI'])

# Generate the plot
pdf(paste(output.dir, "Boxplot_Secretome_DBI.pdf", sep = ''), height = 5, width = 2.5)

    boxplot(DBI ~ Treatment, 
            data = SampleInfo,
            col     = c("dodgerblue", "dodgerblue4"), 
            #ylim    = c(0, 5),
            ylab    = "Normalized log10(counts)",
            outline = F ,
            main    = "DBI")
    
    beeswarm(DBI ~ Treatment, data = SampleInfo, add = T, pch = 16)

    # Add line at control median
    abline(h = Empty.median, lty = 2)

dev.off()


