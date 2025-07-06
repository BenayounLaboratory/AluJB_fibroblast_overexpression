# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(readxl) # To read excel
library(splitstackshape)
library(limma) # For variance stabilizing normalization (VSN)
    
# Define the output directories
abundance.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Processed_Abundances/'
MDS.dir <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/MDS_PCA/'


# CELLULAR PROTEINS ----------------------------------------------------------------------------------------------------------------------------------------------------------------


# ADD GSEA METRIC TO DIFFERENTIAL ABUNDANCE RESULTS

# Load the differential protein results excel file (skip the first 4 lines)
diff.abundances <- as.data.frame(read_excel("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Collaborator_Files_Cellular_Proteins/Candidates_JB1_Pellets_2022_0407_directDIA_v6_Processed_2022_0407_v1.xlsx", sheet = "All_Quantifiable", skip = 4))
    
# Add column with GSEA ranking metric (log2FC * -log(Q))
diff.abundances$Rank_metric <- diff.abundances$`AVG Log2 Ratio` * -log10(diff.abundances$Qvalue)

# Assign uniprotID to rownames
rownames(diff.abundances) <- diff.abundances$UniProtIds
    
# Save data
write.table(diff.abundances, file = paste(abundance.dir, "Cellular_Proteome_Differential_Abundance_Results", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)   



      

# TRANSFORM RAW ABUNDANCES

# Load protein abundance excel file (skip the first 6 lines)
my.abundances <- as.data.frame(read_excel("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Collaborator_Files_Cellular_Proteins/ProteinQuant_JB1_Pellets_2022_0407_directDIA_v6_Processed.xlsx", sheet = "Protein_Abundances", skip = 6))
    
    # Assign uniprotID to rownames
    rownames(my.abundances) <- my.abundances$UniProtIds
    
    # Save raw abundances in txt format
    write.table(my.abundances, file = paste(abundance.dir, "Cellular_Proteome_Abundances_RAW", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)

# Carry out VSN transform
my.abundances.vsn <- as.data.frame(normalizeVSN(as.matrix(my.abundances[, -c(1, 2)])))

    # Add gene IDs back
    my.abundances.vsn <- cbind(my.abundances[, c("UniProtIds", "Genes")], 
                               my.abundances.vsn)
    
    # Save data
    write.table(my.abundances.vsn, file = paste(abundance.dir, "Cellular_Proteome_Abundances_VSN", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    
# Carry out log10 transform
my.abundances.log10 <- as.data.frame(log10(as.matrix(my.abundances[, -c(1, 2)])))

    # Add gene IDs back
    my.abundances.log10 <- cbind(my.abundances[, c("UniProtIds", "Genes")], 
                                 my.abundances.log10)
    
    # Save data
    write.table(my.abundances.log10, file = paste(abundance.dir, "Cellular_Proteome_Abundances_log10", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    



    
# MULTIDIMENTIONAL SCALING (MDS) AND PRINCIPAL COMPONENT ANALYSIS (PCA)

# Define matrix to analyze
abundance.matrix <- my.abundances.log10[, -c(1, 2)]

# Define MDS parameters
point.cex <- 1
label.cex <- 1
axis.cex <- 1
main.cex <- 1
my.colors <- c(rep('dodgerblue', 5), rep('dodgerblue4', 5))
my.shapes <- c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16)
    
# Run MDS
mds.proteins <- cmdscale(1-cor(abundance.matrix, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

# Generate and save MDS plot
pdf(paste(MDS.dir, "Plot_MDS_Cellular_Proteins_AluJb_OE_IMR90", ".pdf", sep=""), width = 8, height = 8)
par(mfrow=c(3,3))

      # Protein data
      plot(x = mds.proteins[, 1],
           y = mds.proteins[, 2],
           pch = my.shapes,
           col = my.colors,
           xlab = "MDS dimension 1",
           ylab = "MDS dimension 2",
           main = paste("AluJb OE, Cellular Protein Abundance MDS", sep = ''),
           #xlim = c(-0.015, 0.015),
           #ylim = c(-0.010, 0.010),
           cex.main = main.cex,
           cex.axis = axis.cex,
           cex.lab = label.cex,
           cex = point.cex)

dev.off()
    

# Check that var for each protein is > 0
var.proteins <- apply(abundance.matrix, 1, var) > 0

# Perform PCA using proteins with variation
PCA.proteins <- prcomp(t(abundance.matrix[var.proteins, ]), scale = TRUE)

# Extract pca stats
summary.PCA.proteins <- summary(PCA.proteins)

# Generate and save plot
pdf(paste(MDS.dir, "Plot_PCA_Cellular_Proteins_AluJb_OE_IMR90", ".pdf", sep=""), width = 8, height = 8)
par(mfrow=c(3,3)) 

    # Protein data
    plot(PCA.proteins$x[, 1],
         PCA.proteins$x[, 2],
         col = my.colors,
         pch = my.shapes,
         main = paste("AluJb OE, Cellular Protein Abundance PCA", sep = ''),
         xlab = paste('PC1 (', round(100*summary.PCA.proteins$importance[,1][2],1),"%)", sep=""),
         ylab = paste('PC2 (', round(100*summary.PCA.proteins$importance[,2][2],1),"%)", sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.main = main.cex,
         cex.lab = label.cex,
         cex.axis = axis.cex,
         cex = point.cex)

dev.off()



# SECRETED PROTEINS ----------------------------------------------------------------------------------------------------------------------------------------------------------------


# ADD GSEA METRIC TO DIFFERENTIAL ABUNDANCE RESULTS

# Load the differential protein results excel file (skip the first 4 lines)
diff.abundances <- as.data.frame(read_excel("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Collaborator_Files_Secreted_Proteins/2022_02_11_Secretome_Candidates_21_1115_JB1_Allsamples_nonorm_directDIA_v5.xlsx", sheet = "JB1_2pept", skip = 4))
    
# Add column with GSEA ranking metric (log2FC * -log(Q))
diff.abundances$Rank_metric <- diff.abundances$`AVG Log2 Ratio` * -log10(diff.abundances$Qvalue)

# Assign uniprotID to rownames
rownames(diff.abundances) <- diff.abundances$UniProtIds
    
# Save data
write.table(diff.abundances, file = paste(abundance.dir, "Secreted_Proteome_Differential_Abundance_Results", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)   



      

# TRANSFORM RAW ABUNDANCES

# Load protein abundance excel file (skip the first 6 lines)
my.abundances <- as.data.frame(read_excel("/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Collaborator_Files_Secreted_Proteins/2022_02_15_Report_Birgit_Protein Quant_Pivot (Pivot)_21_1115_JB1_Allsamples_nonorm_directDIA_v5.xlsx", sheet = "Report_Birgit_Protein Quant_Piv", skip = 4))
    
    # Assign uniprotID to rownames
    rownames(my.abundances) <- my.abundances$PG.UniProtIds
    
    # Save raw abundances in txt format
    write.table(my.abundances, file = paste(abundance.dir, "Secreted_Proteome_Abundances_RAW", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)

# Carry out VSN transform
my.abundances.vsn <- as.data.frame(normalizeVSN(as.matrix(my.abundances[, -c(1:7)])))

    # Add gene IDs back
    my.abundances.vsn <- cbind(my.abundances[, c(1:7)], 
                               my.abundances.vsn)
    
    # Save data
    write.table(my.abundances.vsn, file = paste(abundance.dir, "Secreted_Proteome_Abundances_VSN", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    
# Carry out log10 transform
my.abundances.log10 <- as.data.frame(log10(as.matrix(my.abundances[, -c(1:7)])))

    # Add gene IDs back
    my.abundances.log10 <- cbind(my.abundances[, c(1:7)], 
                                 my.abundances.log10)
    
    # Save data
    write.table(my.abundances.log10, file = paste(abundance.dir, "Secreted_Proteome_Abundances_log10", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    
    

    
    
# MULTIDIMENTIONAL SCALING (MDS) AND PRINCIPAL COMPONENT ANALYSIS (PCA)

# Define matrix to analyze
abundance.matrix <- my.abundances.log10[, -c(1:7)]

# Define MDS parameters
point.cex <- 1
label.cex <- 1
axis.cex <- 1
main.cex <- 1
my.colors <- c(rep('dodgerblue', 4), rep('dodgerblue4', 4))
my.shapes <- c(16, 16, 16, 16, 16, 16, 16, 16)
    
# Run MDS
mds.proteins <- cmdscale(1-cor(abundance.matrix, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

# Generate and save plot
pdf(paste(MDS.dir, "Plot_MDS_Secreted_Proteins_AluJb_OE_IMR90", ".pdf", sep=""), width = 8, height = 8)
par(mfrow=c(3,3))

      # Protein data
      plot(x = mds.proteins[, 1],
           y = mds.proteins[, 2],
           pch = my.shapes,
           col = my.colors,
           xlab = "MDS dimension 1",
           ylab = "MDS dimension 2",
           main = paste("AluJb OE, Secreted Protein Abundance MDS", sep = ''),
           #xlim = c(-0.015, 0.015),
           #ylim = c(-0.010, 0.010),
           cex.main = main.cex,
           cex.axis = axis.cex,
           cex.lab = label.cex,
           cex = point.cex)

dev.off()
    

# Check that var for each protein is > 0
var.proteins <- apply(abundance.matrix, 1, var) > 0

# Perform PCA using proteins with variation
PCA.proteins <- prcomp(t(abundance.matrix[var.proteins, ]), scale = TRUE)

# Extract pca stats
summary.PCA.proteins <- summary(PCA.proteins)

# Generate and save plot
pdf(paste(MDS.dir, "Plot_PCA_Secreted_Proteins_AluJb_OE_IMR90", ".pdf", sep=""), width = 8, height = 8)
par(mfrow=c(3,3)) 

    # Protein data
    plot(PCA.proteins$x[, 1],
         PCA.proteins$x[, 2],
         col = my.colors,
         pch = my.shapes,
         main = paste("AluJb OE, Secreted Protein Abundance PCA", sep = ''),
         xlab = paste('PC1 (', round(100*summary.PCA.proteins$importance[,1][2],1),"%)", sep=""),
         ylab = paste('PC2 (', round(100*summary.PCA.proteins$importance[,2][2],1),"%)", sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.main = main.cex,
         cex.lab = label.cex,
         cex.axis = axis.cex,
         cex = point.cex)

dev.off()



# SESSION INFO ----------------------------------------------------------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/3_Mass_Spec_Analyses/1_Differential_Abundance_Analysis/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_Process_Mass_Spec_and_MDS_PCA.txt", sep =""))
sessionInfo()
sink()      
    
    
# Clean the environment
rm(list=ls())