aggregate_counts <- function(cntTable_complete.list) {
  
  # cntTable_complete.list should be a list of counts files
  
  
  
  # Read 1 count table and keep only the gene column. In the FOR loop below, extend by adding counts from each sample. 
  count_data <- read.csv(cntTable_complete.list[1], header=TRUE, sep="")
  count_data <- count_data[1]
  
  
  for (count_table_file in cntTable_complete.list) {
    
    # Load individual table temporarily. 
    count_data_temp <- read.csv(count_table_file, header=TRUE, sep="")
    
    # Check whether the GENE/TE names column matches that of the main count_data table.  
    column_1_match <- all.equal(count_data_temp[1], count_data[1])
    
    if (column_1_match == TRUE) {
      
      # Add the temporary counts data to the final counts object.
      count_data <- data.frame(count_data, count_data_temp[2])
      
      # Remove count_data_temp since it is not needed
      rm(count_data_temp)
      
    } else {
      
      count_data <- c('Count table gene names do not match in some samples')
      break
      
    } # End if else.
    
  } #End FOR loop
  
  
  
  # Output the final count data table.
  return(count_data)
  
  
  
  
} # END FUNCTION

calculate_IQR_bounds <- function(values) {
  
  # Calculate the first and third quartiles
  Q1 <- quantile(values, 0.25)
  Q3 <- quantile(values, 0.75)
  
  # Calculate the Interquartile Range (IQR)
  IQR <- Q3 - Q1
  
  # Calculate the 1.5 * IQR bounds
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # Return the bounds as a list
  return(list(lower_bound = lower_bound, upper_bound = upper_bound))
}

prefilter_TE_aggregation <- function(gene_TE_counts, aggregate_all, intergenic_nearby_repeats, intergenic_distal_repeats, intronic_repeats, exonic_repeats) {
  
  # FUNCTION INFO: THIS FUNCTION AGGREGATES TE COUNTS BY the specified TE group (subfamily, family, class).
  # Aggregation with this function is prior to the filtering step, since filtering removes too many loci


  
  
  # Assign name to column with gene/TE names.
  names(gene_TE_counts)[1] <- 'TEs'  
  
  # Subset genes
  gene_counts <- gene_TE_counts[grepl('ENSG', gene_TE_counts[,1]), ]
  
  # Subset TEs (shouldn't have ensembl prefix)
  TE_counts <- gene_TE_counts[!grepl('ENSG', gene_TE_counts[,1]), ]
  
      # Print out how many TEs remain
      message(paste(nrow(TE_counts), 'TE loci with at least 1 count in 1 sample'))
  
  # Split TE column. Split columns are placed at end of the matrix.
  TE_counts <- as.data.frame(splitstackshape::cSplit(TE_counts, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = TRUE, stripWhite = FALSE, makeEqual = TRUE))
  
      # Re-add the TE columns to the start of the df and delete the TE columns at the end of the df
      TE_counts <- cbind(TE_counts[, c('TEs_1', 'TEs_2', 'TEs_3', 'TEs_4')], TE_counts[, 1:(ncol(TE_counts)-4)])
    
      # Re-label TE columns
      names(TE_counts)[1:4] <- c('Locus', 'Subfamily', 'Family', 'Class')
    
      # Define label indices (columns that have label)
      label_cols <- c(1:4)
      
      # Assign locus name to the rownames
      rownames(TE_counts) <- TE_counts$Locus
      
      
  # Subset loci of interest and aggregate by subfamily
  if (aggregate_all == TRUE) { # aggregate all loci without filtering
    
      # Generate new names in format Subfamily:Family:Class
      final_names <- paste(TE_counts$Subfamily, TE_counts$Family, TE_counts$Class, sep = ':')
      
      # Aggregate counts by subfamily
      TE_counts.family.sum <- aggregate(x = TE_counts[, -label_cols], by = list(final_names), FUN = 'sum')
      
      # Name TE column; rbind doesn't work without this.
      names(TE_counts.family.sum)[1] <- 'TEs' 
      
      # Recombine gene and TE counts
      gene_TE_counts.sum <- rbind(gene_counts, TE_counts.family.sum)
    
  } else {
    
    # Group expressed repeats by genic location
    TE_counts.intergenic.nearby <- TE_counts[which(rownames(TE_counts) %in% rownames(intergenic_nearby_repeats)), ]
    TE_counts.intergenic.distal <- TE_counts[which(rownames(TE_counts) %in% rownames(intergenic_distal_repeats)), ]
    TE_counts.intronic <- TE_counts[which(rownames(TE_counts) %in% rownames(intronic_repeats)), ]
    TE_counts.exonic <- TE_counts[which(rownames(TE_counts) %in% rownames(exonic_repeats)), ]
    
    # Generate new names in format group:Subfamily:Family:Class
    intergenic.nearby.names <- paste('intergenic_near', TE_counts.intergenic.nearby$Subfamily, TE_counts.intergenic.nearby$Family, TE_counts.intergenic.nearby$Class, sep = ':')
    intergenic.distal.names <- paste('intergenic_distal', TE_counts.intergenic.distal$Subfamily, TE_counts.intergenic.distal$Family, TE_counts.intergenic.distal$Class, sep = ':')
    intronic.names <- paste('intronic', TE_counts.intronic$Subfamily, TE_counts.intronic$Family, TE_counts.intronic$Class, sep = ':')
    exonic.names <- paste('exonic', TE_counts.exonic$Subfamily, TE_counts.exonic$Family, TE_counts.exonic$Class, sep = ':')
    
    # Aggregate counts by subfamily
    TE_counts.intergenic.nearby <- aggregate(x = TE_counts.intergenic.nearby[, -label_cols], by = list(intergenic.nearby.names), FUN = 'sum')
    TE_counts.intergenic.distal <- aggregate(x = TE_counts.intergenic.distal[, -label_cols], by = list(intergenic.distal.names), FUN = 'sum')
    TE_counts.intronic <- aggregate(x = TE_counts.intronic[, -label_cols], by = list(intronic.names), FUN = 'sum')
    TE_counts.exonic <- aggregate(x = TE_counts.exonic[, -label_cols], by = list(exonic.names), FUN = 'sum')
    
    # Name the TE column; rbind doesn't work without this.
    names(TE_counts.intergenic.nearby)[1] <- 'TEs' 
    names(TE_counts.intergenic.distal)[1] <- 'TEs' 
    names(TE_counts.intronic)[1] <- 'TEs' 
    names(TE_counts.exonic)[1] <- 'TEs' 
    
    # Recombine gene and TE counts
    gene_TE_counts.sum <- rbind(gene_counts, TE_counts.intergenic.nearby, TE_counts.intergenic.distal, TE_counts.intronic, TE_counts.exonic)
  
  }


      
  
  # Return the updated counts
  return(gene_TE_counts.sum)

      
  
} # END FUNCTION

cleanup_and_filter_counts <- function(count_data, filter, min_counts_per_sample, fraction_of_samples) {
  
  # FUNCTION INFO:
  # CLEANUP: REMOVE GENE/TRANSCRIPT VERSION AND COMBINE COUNTS FROM SIMILAR GENES/TRANSCRIPTS
  # FILTER: LOWLY EXPRESSED GENES
  
  

  # Assign name to columns with gene/TE names. Name needs to specifically be 'gene.TE' since that is referenced.
  names(count_data)[1] <- 'gene.TE'    

  
  # Remove transcript or ENSEMBL version info. Keeping this can raise issues with gene name mapping.
  count_data[[1]] <- sub("\\..*", "", count_data[[1]], fixed=FALSE)
  
  
  # Collect same transcripts into one row. NOTE, pseudoautosomal PAR_Y genes will be combined with genes in homologous regions.
  gene_counts_summed <- aggregate(count_data[,-c(1)],
                                  by = list(count_data$gene.TE),
                                  FUN = 'sum')
  
  # Update the counts matrix by moving the gene/TE column to rownames and deleting the original column
  rownames(gene_counts_summed) <- gene_counts_summed[, 1]
  gene_counts_summed <- gene_counts_summed[, c(-1)]
  
  
  # Remove lowly expressed genes, if desired. Filter by CPM (about 10 reads across all samples) to control for differences in library depth.
  if (filter == TRUE) {
    
    # Calculate the number of samples in the counts table
    num_of_samples <- ncol(count_data) - 1
    
    # Define the minimum number of samples that have to meet the CPM threshold
    min.samples <- round(fraction_of_samples * num_of_samples)
    
    # Calculate library sizes for each sample 
    library_sizes <- colSums(gene_counts_summed)
    
    # Calculate the median library size
    MedianLibSize <- median(library_sizes)
    
    # calculate the CPM cutoff for the input libraries (changes with the median library size)
    CPM.Cutoff <- min_counts_per_sample / MedianLibSize*1e6
    
    # Convert counts to CPMs
    CPM_table <- edgeR::cpm(gene_counts_summed, lib.size = library_sizes)
    
    # Make a filter function, where "min.samples" have to exceed "CPM.Cutoff"
    cpm.filter.func <- genefilter::kOverA(min.samples, CPM.Cutoff)
    
    # Bind the filtering function to filterfun
    flist <- genefilter::filterfun(cpm.filter.func)
    
    # Identify genes that pass the filtering threshold
    kept_genes <- genefilter::genefilter(CPM_table, flist)
    
    # Get counts for the genes that pass filtering
    gene_counts_summed_filtered <- gene_counts_summed[kept_genes, ]
    
    # Return the cleaned and filtered counts
    return(gene_counts_summed_filtered)
    
  } else {
    
    # Return the cleaned counts
    return(gene_counts_summed)
    
  }
  
  

  
} # END FUNCTION

Ensembl_to_symbol <- function(organism, Ensembl_genes){
  
  # THIS FUNCTION MAPS AN ENSEMBL GENE NAME TO ITS SYMBOL
  
  
  
  # Define organism specific parameters
  if (organism == 'hs') {
    organism_dataset <- "hsapiens_gene_ensembl"
    Ensembl_version <- 110
    
  } else if (organism == 'mm') {
    organism_dataset <- "mmusculus_gene_ensembl"
    Ensembl_version <- 99 # DOUBLE CHECK
  }
  
  
  
  # Retrieve ensembl database
  ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset = organism_dataset, host = 'https://www.ensembl.org', version = Ensembl_version, verbose = TRUE)
  
  # Map Ensembl names to gene symbols
  Mapped_gene_info <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "entrezgene_id"), 
                                     filters = c("ensembl_gene_id"), 
                                     values = Ensembl_genes,
                                     mart = ensembl,
                                     verbose = FALSE,
                                     uniqueRows = TRUE)  
  
  # If genes have multiple mappings, keep only the first. 
  Mapped_gene_info <- Mapped_gene_info[!duplicated(Mapped_gene_info$ensembl_gene_id), ]
  
  # Assign Ensembl names to the rownames
  rownames(Mapped_gene_info) <- Mapped_gene_info$ensembl_gene_id
  
  # Remove ensembl gene id column since they're in rownames
  Mapped_gene_info <- Mapped_gene_info[, -c(1)] 
  
  # Return table with the mapping info
  return(Mapped_gene_info)
  
  
} # END FUNCTION

Run_DESeq <- function(filtered_counts, padj_limit, control_label, control_reps, treatment_label, treatment_reps, VST_output, DESeq_output){
  
# padj_limit defines alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

  
  
  


# Define the treatments
treatment.group <- c(rep(control_label, control_reps), 
                     rep(treatment_label, treatment_reps)
                     )

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(filtered_counts),
                         Treatment = treatment.group)


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                                  colData = SampleInfo,
                                  design = ~ Treatment) 

# run DESeq2
dds <- DESeq(dds, parallel = TRUE)

# Get VST data (includes size factor normalization)
my.VST <- assay(varianceStabilizingTransformation(dds, blind = FALSE))

     # Save VST data
     write.table(my.VST, file = paste(VST_output, 'All_counts_filtered_', treatment_label, '_VST', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
     



# Extract DESeq results. 
DESeq.Effects <- results(dds, contrast = c("Treatment", treatment_label, control_label), alpha = padj_limit, independentFiltering = TRUE)


# DESeq Stats
summary(DESeq.Effects, alpha = padj_limit)



# Generate MA plots
pdf(paste(DESeq_output,"Plot_MA_DESeq_Results_", treatment_label, ".pdf", sep=""), height = 10, width = 10)

    DESeq2::plotMA(object = DESeq.Effects, alpha = padj_limit, main = 'MA Plot', colSig = 'orange', ylim = c(-4,4))

dev.off()





# Run function to map ensembl labels, extract significant results, and save them.
my.sig <- Extract_DESeq_stats(DESeq_results = DESeq.Effects, padj_limit = padj_limit, organism = 'hs', 
                              output.dir = DESeq_output,
                              output_file_prefix_all = paste('All_Genes_', treatment_label, sep = ''),
                              output_file_prefix_sig = paste('FDR_', padj_limit, '_Genes_', treatment_label, sep = ''))

  
  

# Return the VST data and significant results as a list

return(list(my.VST, my.sig))
  
  
  
} # END FUNCTION

Extract_DESeq_stats <- function(DESeq_results, padj_limit, organism, output.dir, output_file_prefix_all, output_file_prefix_sig){
  
  # FUNCTION NOTE: Extract statistics from DESeq Results object, remove genes with NA for padj, save padj filtered and unfiltered list. The unfiltered list can be used as gene background.
  
  
  
  # Check if output directory exists. If not, make it. 
  if ( dir.exists(output.dir) == FALSE) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  
  
  # Extract statistics table from DESeq results
  DESeq_stats <- as.data.frame(DESeq_results@listData)
  
  # Assign rownames (genes) to statistics table
  rownames(DESeq_stats) <- DESeq_results@rownames
  
  # Remove rows with NAs
  DESeq_stats <- na.omit(DESeq_stats)
  
  # Getting alternative gene names for remaining genes
  alternative_names <- Ensembl_to_symbol(organism = organism, Ensembl_genes = rownames(DESeq_stats))
  
  # Make columns to hold alternative gene names
  DESeq_stats$external_gene_name <- NA 
  DESeq_stats$hgnc_symbol <- NA
  DESeq_stats$entrezgene_id <- NA
  
  # Assign alternative names
  DESeq_stats[rownames(alternative_names),c('external_gene_name', 'hgnc_symbol', 'entrezgene_id')] <- alternative_names[,]

  # Filter out DEGs
  DESeq_stats.sig <- DESeq_stats[DESeq_stats$padj < padj_limit, ]
  
  # Save stats for the full and significant gene names
  write.table(DESeq_stats.sig, file = paste(output.dir, output_file_prefix_sig, '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
  write.table(DESeq_stats, file = paste(output.dir, output_file_prefix_all, '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

  
  
  # Output the significant gene stats
  return(DESeq_stats.sig)
  
} # END FUNCTION

Run_MDS <-function(VST_expression, plot_colors, plot_shapes, treatment_label, MDS_dir){
  
  
  
    
# Generate Gene-Only, TE-Only, LINE-only, and L1-only expression to assess their ability to stratify samples
Gene_expr <- VST_expression[grepl('ENSG', rownames(VST_expression)), ]
TE_expr <- VST_expression[!grepl('ENSG', rownames(VST_expression)), ]

# Run MDS
mds.genes <- cmdscale(1-cor(Gene_expr, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
mds.TE <- cmdscale(1-cor(TE_expr, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

# Define MDS parameters
mds.point.cex <- 1
mds.label.cex <- 1
mds.axis.cex <- 1
mds.main.cex <- 1

# Generate and save plot
pdf(paste(MDS_dir, "Plot_MDS_", treatment_label, ".pdf", sep=""), width = 8, height = 8)
par(mfrow=c(3,3))

      # Gene data
      plot(x = mds.genes[, 1],
           y = mds.genes[, 2],
           pch = plot_shapes,
           col = plot_colors,
           xlab = "MDS dimension 1",
           ylab = "MDS dimension 2",
           main = paste(treatment_label, " Gene Expression MDS", sep = ''),
           #xlim = c(-0.015, 0.015),
           #ylim = c(-0.010, 0.010),
           cex.main = mds.main.cex,
           cex.axis = mds.axis.cex,
           cex.lab = mds.label.cex,
           cex = mds.point.cex)

      # TE data
      plot(x = mds.TE[, 1],
           y = mds.TE[, 2],
           pch = plot_shapes,
           col = plot_colors,
           xlab = "MDS dimension 1",
           ylab = "MDS dimension 2",
           main = paste(treatment_label, " Repeat Expression MDS", sep = ''),
           #xlim = c(-0.02, 0.02),
           #ylim = c(-0.02, 0.02),
           cex.main = mds.main.cex,
           cex.axis = mds.axis.cex,
           cex.lab = mds.label.cex,
           cex = mds.point.cex)

dev.off()
  
  

# Return nothing
return(NULL)
  
  
  
} # END FUNCTION

Run_PCA <-function(VST_expression, plot_colors, plot_shapes, treatment_label, MDS_dir){
  
  
  
    
# Generate Gene-Only, TE-Only, LINE-only, L1-only expression to assess their ability to stratify samples
Gene_expr <- VST_expression[grepl('ENSG', rownames(VST_expression)), ]
TE_expr <- VST_expression[!grepl('ENSG', rownames(VST_expression)), ]

# Check that var for each gene is > 0
var.genes <- apply(Gene_expr, 1, var) > 0
var.TEs <- apply(TE_expr, 1, var) > 0

# Perform PCA using genes with variation
PCA.genes <- prcomp(t(Gene_expr[var.genes, ]), scale = TRUE)
PCA.TEs <- prcomp(t(TE_expr[var.TEs, ]), scale = TRUE)

# Extract pca stats
summary.PCA.genes <- summary(PCA.genes)
summary.PCA.TEs <- summary(PCA.TEs)

# Define PCA parameters
mds.point.cex <- 1
mds.label.cex <- 1
mds.axis.cex <- 1
mds.main.cex <- 1


# Generate and save plot
pdf(paste(MDS_dir, "Plot_PCA_", treatment_label, ".pdf", sep=""), width = 8, height = 8)
par(mfrow=c(3,3)) # Added to produce images of similar shape across analyses

    # PCA for Genes ONLY
    plot(PCA.genes$x[, 1],
         PCA.genes$x[, 2],
         col = plot_colors,
         pch = plot_shapes,
         main = paste(treatment_label, " Gene Expression PCA", sep = ''),
         xlab = paste('PC1 (', round(100*summary.PCA.genes$importance[,1][2],1),"%)", sep=""),
         ylab = paste('PC2 (', round(100*summary.PCA.genes$importance[,2][2],1),"%)", sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.main = mds.main.cex,
         cex.lab = mds.label.cex,
         cex.axis = mds.axis.cex,
         cex = mds.point.cex)

    # PCA for TEs ONLY
    plot(PCA.TEs$x[, 1],
         PCA.TEs$x[, 2],
         col = plot_colors,
         pch = plot_shapes,
         main = paste(treatment_label, " Repeat Expression PCA", sep = ''),
         xlab = paste('PC1 (', round(100*summary.PCA.TEs$importance[,1][2],1),"%)", sep=""),
         ylab = paste('PC2 (', round(100*summary.PCA.TEs$importance[,2][2],1),"%)", sep=""),
         #xlim = c(-75, 75),
         #ylim = c(-25, 25),
         cex.main = mds.main.cex,
         cex.lab = mds.label.cex,
         cex.axis = mds.axis.cex,
         cex = mds.point.cex)

dev.off()





  

# Return nothing
return(NULL)
  
  
  
} # END FUNCTION
