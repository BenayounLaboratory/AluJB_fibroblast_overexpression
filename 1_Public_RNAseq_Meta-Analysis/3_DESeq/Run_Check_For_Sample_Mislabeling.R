# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/Process_RNASeq_functions.R")

# Load libraries
library(DESeq2)


    

  
# Define the output directory
dir.output <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/Processed_counts/'

# Define the directory with the counts tables and collect individual file paths in a list.
count_Table.dir <- c("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/2_Read_counting/TETranscripts_counts_GSE113957_Human_Dermal_Fibroblasts/")
cntTable_list <- list.files(count_Table.dir, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables into one df.
gene_TE_counts_raw <- aggregate_counts(cntTable_complete.list = cntTable_list)

# Remove useless info. from sample names
colnames(gene_TE_counts_raw) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub(".sra*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub("*fastp_", "", colnames(gene_TE_counts_raw))

# Load sample metadata
my.meta_data <- read.csv("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/0_Metadata_and_External_Resources/GSE113957_Metadata_Aging_Subset2_White-Caucasian_ArmSkin.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
    
    # Update instrument names to remove non-letter characters
    my.meta_data[which(my.meta_data$Instrument == c("NextSeq 500")), 'Instrument'] <- 'NextSeq'
    my.meta_data[which(my.meta_data$Instrument == c("Illumina HiSeq 2500")), 'Instrument'] <- 'HiSeq'
    
    # Add column for age group (by 10s)
    my.meta_data$Age_Group1 <- 'Fill'
    
        # Assign samples to groups
        my.meta_data[which(my.meta_data$Age >= 20 & my.meta_data$Age <= 29), 'Age_Group1'] <- 'A_20s'
        my.meta_data[which(my.meta_data$Age >= 30 & my.meta_data$Age <= 39), 'Age_Group1'] <- 'B_30s'
        my.meta_data[which(my.meta_data$Age >= 40 & my.meta_data$Age <= 49), 'Age_Group1'] <- 'C_40s'
        my.meta_data[which(my.meta_data$Age >= 50 & my.meta_data$Age <= 59), 'Age_Group1'] <- 'D_50s'
        my.meta_data[which(my.meta_data$Age >= 60 & my.meta_data$Age <= 69), 'Age_Group1'] <- 'E_60s'
        my.meta_data[which(my.meta_data$Age >= 70 & my.meta_data$Age <= 79), 'Age_Group1'] <- 'F_70s'
        my.meta_data[which(my.meta_data$Age >= 80 & my.meta_data$Age <= 89), 'Age_Group1'] <- 'G_80s'
        my.meta_data[which(my.meta_data$Age >= 90 & my.meta_data$Age <= 99), 'Age_Group1'] <- 'H_90s'

        # Update the column as a factor
        my.meta_data$Age_Group1 <- as.factor(my.meta_data$Age_Group1)
        
    # Add column for age group (by broader brackets)
    my.meta_data$Age_Group2 <- 'Fill'
    
        # Assign samples to groups
        my.meta_data[which(my.meta_data$Age >= 20 & my.meta_data$Age <= 39), 'Age_Group2'] <- 'A_Young'
        my.meta_data[which(my.meta_data$Age >= 40 & my.meta_data$Age <= 59), 'Age_Group2'] <- 'B_Middle_Aged'
        my.meta_data[which(my.meta_data$Age >= 60 & my.meta_data$Age <= 74), 'Age_Group2'] <- 'C_Old'
        my.meta_data[which(my.meta_data$Age >= 75 & my.meta_data$Age <= 99), 'Age_Group2'] <- 'D_Geriatric'

        # Update the column as a factor
        my.meta_data$Age_Group2 <- as.factor(my.meta_data$Age_Group2)
        
    # Keep metadata for covariates only
    my.meta_data <- my.meta_data[, c('Age', 'sex', 'Instrument', 'source_name', 'ETHNICITY', 'RACE', 'Age_Group1', 'Age_Group2')]
    
    # Keep sample metadata only for samples that were aligned
    my.meta_data <- my.meta_data[which(rownames(my.meta_data) %in% colnames(gene_TE_counts_raw)), ]
    
# Remove transcript or ENSEMBL gene version info.
gene_TE_counts_raw[[1]] <- sub("\\..*", "", gene_TE_counts_raw[[1]], fixed=FALSE)    

# Assign gene names to rownames
rownames(gene_TE_counts_raw) <- gene_TE_counts_raw$gene.TE

# Reorder gene expression columns (keep only the final sample set, organized by age)
gene_TE_counts_raw <- gene_TE_counts_raw[, c(rownames(my.meta_data))]

# Remove rows with all 0
gene_TE_counts_raw <- gene_TE_counts_raw[which(rowSums(gene_TE_counts_raw) > 0), ]





# QUALITY CONTROL: Check male/female specific gene expression

# Define sample names
my.samples <- rownames(my.meta_data)

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = rownames(my.meta_data),
                         Age = my.meta_data$Age_Group2,
                         Sex = my.meta_data$sex,
                         Instrument = my.meta_data$Instrument
                         )

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = gene_TE_counts_raw,
                              colData = SampleInfo,
                              design = ~ Age + Sex + Instrument) 

# run DESeq2 normalization
dds <- DESeq(dds, parallel = TRUE)

# Get VST data
normalized_counts <- getVarianceStabilizedData(dds)

# Extract Xist/Ddx3y expression
my.xist <- normalized_counts['ENSG00000229807', my.samples]
my.ddx3y <- normalized_counts['ENSG00000067048', my.samples]

# Combine into a dataframe
sex.marker.expression <- data.frame(Female_marker = my.xist,
                                    Male_marker = my.ddx3y,
                                    Sex = my.meta_data[my.samples, 'sex'])

# Update the colnames
colnames(sex.marker.expression) <- c('Xist', 'Ddx3y', 'Sex')

# Organize Sex column alphabetically
sex.marker.expression <- sex.marker.expression[order(sex.marker.expression$Sex), ]

# Lock in factor level order
sex.marker.expression$Sex <- factor(sex.marker.expression$Sex, levels = unique(sex.marker.expression$Sex))

# Xist Plot (PRIOR TO QC FILTER)
pdf(paste(dir.output, "GSE113957_Boxplot_Expression_Xist_PREFILTER", ".pdf", sep=""), width = 4.4, height = 4.5)

        # make boxplots 
        boxplot(Xist ~ Sex, 
                data = sex.marker.expression,
                ylim = c(0, 15),
                outline = F, 
                main = paste('Xist Expression', sep = ''),
                ylab = paste('Normalized log2(counts)', sep = ''),
                xlab = c(''),
                xaxt = 'n')
        
        # overlay the points
        stripchart(Xist ~ Sex, data = sex.marker.expression, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.5, col = 'red') 
        
        # Specify axis labels
        axis(1, at = c(1, 2), labels = c('Female', 'Male'))

# End pdf
dev.off()

# Ddx3y Plot (PRIOR TO QC FILTER)
pdf(paste(dir.output, "GSE113957_Boxplot_Expression_Ddx3y_PREFILTER", ".pdf", sep=""), width = 4.4, height = 4.5)

        # make boxplots 
        boxplot(Ddx3y ~ Sex, 
                data = sex.marker.expression,
                ylim = c(0, 15),
                outline = F, 
                main = paste('Ddx3y Expression', sep = ''),
                ylab = paste('Normalized log2(counts)', sep = ''),
                xlab = c(''),
                xaxt = 'n')
        
        # overlay the points
        stripchart(Ddx3y ~ Sex, data = sex.marker.expression, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.5, col = 'red') 
        
        # Specify axis labels
        axis(1, at = c(1, 2), labels = c('Female', 'Male'))

# End pdf
dev.off()



# Calculate outlier bounds (1.5*IQR for Xist and Ddx3y)
xist_female_bounds <- calculate_IQR_bounds(sex.marker.expression[which(sex.marker.expression$Sex == 'female'), 'Xist'])
xist_male_bounds <- calculate_IQR_bounds(sex.marker.expression[which(sex.marker.expression$Sex == 'male'), 'Xist'])
ddx3y_female_bounds <- calculate_IQR_bounds(sex.marker.expression[which(sex.marker.expression$Sex == 'female'), 'Ddx3y'])
ddx3y_male_bounds <- calculate_IQR_bounds(sex.marker.expression[which(sex.marker.expression$Sex == 'male'), 'Ddx3y'])

# Define samples to remove
removable.samples <- c('SRR7093907', # labeled female, but low Xist and high Ddx3y. Also not Caucasian.
                       'SRR7093870', # labeled male, but high Xist expression between male and female
                       'SRR7093857', # labeled male, but high Xist expression between male and female
                       'SRR7093830', # labeled female, but high Ddx3y between male and female
                       'SRR7093859', # labeled male, but low Ddx3y between male and female
                       'SRR7093866', # borderline outside 1.5*IQR
                       'SRR7093818') # borderline outside 1.5*IQR


# Remove bad samples from the metadata
my.meta_data_filtered <- my.meta_data[!(rownames(my.meta_data) %in% removable.samples), ]
    
    # Save modified metadata
    write.table(my.meta_data_filtered, file = paste("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/0_Metadata_and_External_Resources/GSE113957_Metadata_Aging_Subset2_White-Caucasian_ArmSkin_CLEAN.txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)

# Update the QC dataframe to recheck plots
sex.marker.expression <- sex.marker.expression[!(rownames(sex.marker.expression) %in% removable.samples), ]

# Lock in factor level order
sex.marker.expression$Sex <- factor(sex.marker.expression$Sex, levels = unique(sex.marker.expression$Sex))

# Xist Plot (POST QC FILTER)
pdf(paste(dir.output, "GSE113957_Boxplot_Expression_Xist_POSTFILTER", ".pdf", sep=""), width = 4.4, height = 4.5)

        # make boxplots 
        boxplot(Xist ~ Sex, 
                data = sex.marker.expression,
                ylim = c(0, 15),
                outline = F, 
                main = paste('Xist Expression', sep = ''),
                ylab = paste('Normalized log2(counts)', sep = ''),
                xlab = c(''),
                xaxt = 'n')
        
        # overlay the points
        stripchart(Xist ~ Sex, data = sex.marker.expression, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.5, col = 'red') 
        
        # Specify axis labels
        axis(1, at = c(1, 2), labels = c('Female', 'Male'))


# End pdf
dev.off()

# Ddx3y Plot (POST QC FILTER)
pdf(paste(dir.output, "GSE113957_Boxplot_Expression_Ddx3y_POSTFILTER", ".pdf", sep=""), width = 4.4, height = 4.5)

        # make boxplots 
        boxplot(Ddx3y ~ Sex, 
                data = sex.marker.expression,
                ylim = c(0, 15),
                outline = F, 
                main = paste('Ddx3y Expression', sep = ''),
                ylab = paste('Normalized log2(counts)', sep = ''),
                xlab = c(''),
                xaxt = 'n')
        
        # overlay the points
        stripchart(Ddx3y ~ Sex, data = sex.marker.expression, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.5, col = 'red') 
        
        # Specify axis labels
        axis(1, at = c(1, 2), labels = c('Female', 'Male'))


# End pdf
dev.off()
    
    

    



# Session Info ------------------------------------------------------------------------------------------------------------------------



# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/1_Public_RNAseq_Meta-Analysis/3_DESeq/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_Check_Mislabeling.txt", sep =""))
sessionInfo()
sink()      
    


    
# Clean the environment
rm(list=ls())


