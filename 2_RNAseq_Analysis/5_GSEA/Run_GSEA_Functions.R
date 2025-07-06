make_genelist <- function(DESeq_results_path) {
  
  
  # Load DESeq results
  my.DESeq.res <- read.csv(DESeq_results_path, header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

  # Assign EnsembleID to entries with an NA or Empty cell for HGNC Symbol
  my.DESeq.res[is.na(my.DESeq.res$hgnc_symbol), 'hgnc_symbol'] <- rownames(my.DESeq.res[is.na(my.DESeq.res$hgnc_symbol), ])
  my.DESeq.res[which(my.DESeq.res$hgnc_symbol == ''), 'hgnc_symbol'] <- rownames(my.DESeq.res[which(my.DESeq.res$hgnc_symbol == ''), ])
  
  # Remove duplicate gene names
  my.DESeq.res <- my.DESeq.res[!duplicated(my.DESeq.res$hgnc_symbol), ]
  
  # Define genelists using the Wald stat
  my.genelist <- my.DESeq.res$stat
  
  # Assign gene name or repeat subfamily to the list
  names(my.genelist) <- my.DESeq.res$hgnc_symbol
   
  # Sort the genelist, largest to smallest
  my.genelist = sort(my.genelist, decreasing = TRUE)
  
  
  # Output ranked gene list
  return(my.genelist)
  
  
  
} # END FUNCTION

make_genelist_Ensemble <- function(DESeq_results_path) {
  
  
  # Load DESeq results
  my.DESeq.res <- read.csv(DESeq_results_path, header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

  # Define genelists using the Wald stat
  my.genelist <- my.DESeq.res$stat
  
  # Assign gene names to the list
  names(my.genelist) <- rownames(my.DESeq.res)
  
  # Sort the genelist, largest to smallest
  my.genelist = sort(my.genelist, decreasing = TRUE)
  
  # Output ranked gene list
  return(my.genelist)
  
  
  
} # END FUNCTION

Run_GSEA <- function(output.dir, output_subfolder, my.genelist, my.gs.collection, gs.label, condition_label, output.object, my.max_gs_size = 500, items_to_plot = 5, plot_non_sig = FALSE) {
  
  
  # Define the output directory
  dir.output <- paste(output.dir, output_subfolder, sep = '')
  
  # Set seed for reproducibility
  set.seed(90007)
  
  # Run GSEA (ignore warnings about 1: ties in the ranking stat, assuming they are only a very small fraction)
  my.GSEA <- GSEA(geneList = my.genelist, 
                   TERM2GENE = my.gs.collection,  
                   minGSSize = 15, 
                   maxGSSize = my.max_gs_size, 
                   eps = 1e-50,
                   nPermSimple = 1e5,
                   pvalueCutoff = 1,
                   pAdjustMethod = 'BH',
                   by = 'fgsea',
                   seed = TRUE)
      
  # Save *ALL* results, regardless of significance
  write.table(my.GSEA@result, file = paste(dir.output, 'ClusterProfiler_GSEA_Table_', gs.label, '_', condition_label, '_ALL', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

  # Make an object to hold results filtered on significance
  my.GSEA.sig <- my.GSEA
  
  # Filter results on FDR < 0.05
  my.GSEA.sig@result <- my.GSEA.sig@result[which(my.GSEA.sig@result$p.adjust < 0.05), ]
  
  # Define the number of significant results
  my.sig.gs.num <- nrow(my.GSEA.sig@result)

  # if there are significant results, save files
  if (my.sig.gs.num > 0) {

    # write results to file
    write.table(my.GSEA.sig@result, file = paste(dir.output, 'ClusterProfiler_GSEA_Table_', gs.label, '_', condition_label, '_FDR5_', my.sig.gs.num, '_Significant' ,'.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

    # Make and save dotplot with top genesets
    pdf(paste(dir.output, 'ClusterProfiler_GSEA_Dotplot_FDR5_', gs.label, '_', condition_label, ".pdf", sep=""), width = 9, height = 8)

      print(dotplot(my.GSEA.sig, x = TRUE, showCategory = items_to_plot, size = '-log10(p.adjust)', split = '.sign', color = 'NES', title = paste("GSEA ", gs.label, sep = ''), font.size = 12) +
            scale_size(range = c(4, 12)) + 
            aes(shape = I(16)) + 
            aes(color = `NES`) + 
            scale_fill_gradientn(colors = c("blue","white","red"), values = rescale(c(-6,0,6)), guide = "colorbar", limits=c(-6.1, 6.1), aesthetics = 'color')
           ) 

    dev.off()

  }
  


  # Make and save a second dotplot with all top genesets, including non-significant
  if (plot_non_sig == TRUE) {
    
    pdf(paste(dir.output, 'ClusterProfiler_GSEA_Dotplot_ALL_', gs.label, '_', condition_label, ".pdf", sep=""), width = 9, height = 8)

      print(dotplot(my.GSEA, x = TRUE, showCategory = items_to_plot, size = '-log10(p.adjust)', split = '.sign', color = 'NES', title = paste("GSEA ", gs.label, sep = ''), font.size = 12) +
            scale_size(range = c(4, 12)) + 
            aes(shape = I(16)) + 
            aes(color = `NES`) + 
            scale_fill_gradientn(colors = c("blue","white","red"), values = rescale(c(-6,0,6)), guide = "colorbar", limits=c(-6.1, 6.1), aesthetics = 'color')
           ) 
      
    dev.off()
    
  }
    
  
  
  
  
  # Output either the filtered or unfiltered GSEA object
  if (output.object == 'unfiltered') {
    return(my.GSEA)
  } else if (output.object == 'filtered') {
    return(my.GSEA.sig)
  }
  
  
  
  
} # END FUNCTION

Run_ORA <- function(output.dir, output_subfolder, universe.list, my.genelist, my.gs.collection, gs.label, condition_label, output.object) {
  
  
  # Define the output directory
  dir.output <- paste(output.dir, output_subfolder, sep = '')
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT
  # enricher uses a weird universe, which is defined as the intersection of the supplied universe list with gene set collection genes. This is problematic when the gene sets don't cover the entire genome
  # Here, add the supplied universe list to the gene set collection. The intersection of the supplied universe list with the gene sets should now be the supplied universe list (which is the background we want)
  # Adjust the max GSSize so it does not include the universe geneset
  
      # Create a df to hold the universe geneset
      universe.geneset <- data.frame(gs_name = 'Universe', gene = universe.list)
      
      # Update the colnames to match the geneset collection
      colnames(universe.geneset) <- colnames(my.gs.collection)
      
      # Add the universe geneset to the geneset collection of interest
      my.gs.collection <- rbind(my.gs.collection, universe.geneset)
  
  # Run ORA
  my.ORA <- enricher(gene = my.genelist,
                     universe      = universe.list,
                     TERM2GENE     = my.gs.collection,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1,
                     minGSSize     = 15,
                     maxGSSize     = 10000)
  
  # if no genes mapped to the gene list, end the function
  if (is.null(my.ORA)) {
    return()
  }
  
  # Save *ALL* results, regardless of significance
  write.table(my.ORA@result, file = paste(dir.output, 'ClusterProfiler_ORA_Table_', gs.label, '_', condition_label, '_ALL', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
  
  # Make an object to hold results filtered on significance
  my.ORA.sig <- my.ORA
  
  # Filter results on FDR < 0.05
  my.ORA.sig@result <- my.ORA.sig@result[which(my.ORA.sig@result$p.adjust < 0.05 & my.ORA.sig@result$qvalue < 0.20), ]
  
  # Define the number of significant results
  my.sig.gs.num <- nrow(my.ORA.sig@result)
  
  # if there are significant results, save files
  if (my.sig.gs.num > 0) {

      # Define the number of terms to plot. Do a maximum of 10 *significant* terms
      if (my.sig.gs.num > 10) {
        terms_to_plot <- 10
      } else {
        terms_to_plot <- my.sig.gs.num
      }
      
      # Define the maximum gene ratio
      generatio_limit <- ceiling(max(as.numeric(my.ORA.sig@result[1:terms_to_plot, 'GeneRatio'])))
      
      # write results to file
      write.table(my.ORA.sig@result, file = paste(dir.output, 'ClusterProfiler_ORA_Table_', gs.label, '_', condition_label, '_FDR5_', my.sig.gs.num, '_Significant' ,'.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
  
      # Make and save dotplot with top genesets
      pdf(paste(dir.output, 'ClusterProfiler_ORA_Dotplot_', gs.label, '_', condition_label, '_FDR5', ".pdf", sep=""), width = 9, height = 8)
  
        print(dotplot(my.ORA.sig, x = TRUE, showCategory = terms_to_plot, size = '-log10(p.adjust)', color = NULL, title = paste(condition_label, " ORA ", gs.label, sep = ''), font.size = 12) + scale_size(range = c(4, 12)) + aes(shape = I(16)) + aes(color = `GeneRatio`) + set_enrichplot_color(type = "color", name = "GeneRatio", colors = c("white","red"), limits=c(-0.01, generatio_limit))) 
  
      dev.off()

  }
  
  # No output needed
  return()
  
  
} # END FUNCTION

find_genesets_with_directionality <- function(input.one, input.two, input_directionality){
  
  # This function is called by the overlap functions below
  # This function compares two columns and outputs a list of the genesets that meet the specified directionality (ie they change in the same or opposite direction)
  
  
  
  # Calculate the product of the NES scores (which will inform whether results are in the same direction or opposite)
  NES.product <- input.one$NES * input.two$NES
  
  # Extract gene sets changing in the same direction, sets changing in the opposite direction, or all gene sets
  if (input_directionality == c('same')) {
    
      # Find gene set indices changing in the same directions
      consistent.genesets <- which(NES.product > 0)
      
      # Extract the id names
      consistent.genesets <- input.one[consistent.genesets, 'ID']
      
      # Output the geneset IDs
      return(consistent.genesets)
    
  } else if (input_directionality == c('opposite')) {
    
      # Find gene set indices changing in the opposite directions
      opposite.genesets <- which(NES.product < 0)
      
      # Extract the id names
      opposite.genesets <- input.one[opposite.genesets, 'ID']
      
      # Output the geneset IDs
      return(opposite.genesets)
      
    
  } else {
    
      # Apply no filter, and return the original gene set ID list
      return(input.one$ID)
    
  }
  
  
  
} # END FUNCTION

Overlap_two_GSEA_res <- function(comparison_directionality,
                                 result.one, result.two,
                                 label.one, label.two,
                                 number_to_plot,
                                 dir.output, gs.label){
  
  # NOTE: This function will keep gene sets that are significant in two sets of GSEA results.
  # comparison_directionality can be 'same' or 'opposite'. Otherwise, gene sets will not be filtered based on directionality
  # Something to potentially include: If no overlap or if no consistent overlap, end function immediately and maybe specify message
  
  

  
  # -log10 transform the adjusted pvalues
  result.one$p.adjust <- -log10(result.one$p.adjust)
  result.two$p.adjust <- -log10(result.two$p.adjust)

  # Add sample group label column
  result.one$Group <- label.one
  result.two$Group <- label.two
  
  # Define sig gene sets in each analysis
  a <- result.one$ID
  b <- result.two$ID
  
  # Find gene set overlaps
  common_genesets <- Reduce(intersect, list(a, b))

  # Extract overlapping genesets and columns necessary for plotting
  result.one <- result.one[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.two <- result.two[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  
  # Run function to find gene set IDs with the desired directionality
  genesets_of_interest <- find_genesets_with_directionality(input.one = result.one, input.two = result.two, input_directionality = comparison_directionality)
  
  # Keep gene sets changing in the desired way
  result.one <- result.one[genesets_of_interest, ]
  result.two <- result.two[genesets_of_interest, ]
  
  # Make a subfunction to calculate Fischer meta-analysis pvalues; this will be used on conjunction with 'apply' to calculate a meta pvalue for each row
  calculate.fisher.p <- function(my.pvalues){
    
      # Extract pvalues from the ith row
      my.pvalue.A = my.pvalues[1]
      my.pvalue.B = my.pvalues[2]

      # Calculate the meta pvalue using Fishers method
      my.meta.pvalue <- poolr::fisher(p = c(my.pvalue.A, my.pvalue.B))
      
      # Output the meta pvalue
      return(my.meta.pvalue$p)
    
  }
  
  # Calculate a combined Fischer pvalue for each gene set
  all.my.meta.pvalues <- apply(cbind(result.one$pvalue, result.two$pvalue), 
                               MARGIN = 1, 
                               calculate.fisher.p)
  
  # Combine the results into a table that will be saved as an output (and that will be different than the table used to plot results)
  my.output.table <- cbind(result.one, 
                           result.two[, c('Group', 'NES', 'pvalue', 'p.adjust')])
  
  # Add the meta pvalue to the output table
  my.output.table$meta.pvalue <- all.my.meta.pvalues
  
  # Reorder gene sets from biggest to smallest meta pvalue (***this is needed if we want the final plots to show the more significant values at the top)
  my.output.table <- my.output.table[order(my.output.table$meta.pvalue, decreasing = TRUE), ]
  
      # write results to file
      write.table(my.output.table, file = paste(dir.output, 'ClusterProfiler_GSEA_Table_Overlapping_Gene_Sets_', gs.label, '_FDR5_', nrow(my.output.table), '_Significant' ,'.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

  
  # Extract the data for the first comparison, which will be used to define reference positive and negative terms
  result.one <- result.one[rownames(my.output.table),]
  
  
  
  
  
  # DEFINE THE NUMBER OF POS AND NEG NES TERMS TO PLOT
  
  # Extract the positive NES gene sets
  positive.genesets <- result.one[which(result.one$NES > 0, ), ]
  
  # Extract the negative NES gene sets
  negative.genesets <- result.one[which(result.one$NES < 0, ), ]
  
  # Define the maximum number of positive or negative terms to plot
  max_to_plot <- number_to_plot / 2
  
  # Define the actual number of positive or negative terms to plot
  if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- max_to_plot
      
      # Define the number of negative terms to plot
      negative_to_plot <- max_to_plot
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the maximum number of negative terms to plot
      max_negative <- number_to_plot - positive_to_plot
      
      # Define the number negative terms to plot, either the maximum possible or the number available (whichever is less)
      negative_to_plot <- min(max_negative, nrow(negative.genesets))
      
  } else if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
      
      # Define the maximum number of positive terms to plot
      max_positive <- number_to_plot - negative_to_plot
      
      # Define the number positive terms to plot, either the maximum possible or the number available (whichever is less)
      positive_to_plot <- min(max_positive, nrow(positive.genesets))
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
  }
  
  # Define empty vectors of gene sets to plot, to update if there are terms to plot
  pos.genesets.to.plot <- c()
  neg.genesets.to.plot <- c()
  
  
  
  
  
  # IDENTIFY GENE SET NAMES TO PLOT
  
  if (positive_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      pos_starting_index <- nrow(positive.genesets) - positive_to_plot + 1
      
      # Define the index for the last geneset to plot
      pos_ending_index <- nrow(positive.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      pos.genesets.to.plot <- positive.genesets[pos_starting_index:pos_ending_index, 'ID']
      
  } 
  
  if (negative_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      neg_starting_index <- nrow(negative.genesets) - negative_to_plot + 1
      
      # Define the index for the last geneset to plot
      neg_ending_index <- nrow(negative.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      neg.genesets.to.plot <- negative.genesets[neg_starting_index:neg_ending_index, 'ID']
      
  } 
  
  
  
  
  
  # Collect the positive and negative gene set names to plot
  genesets.to.plot <- c(neg.genesets.to.plot, pos.genesets.to.plot)
  
  # Combine results into 1 dataframe in the *melted* format ggplot needs
  melted_format_data <- rbind(result.one[genesets.to.plot, ], 
                              result.two[genesets.to.plot, ])
  
  # Change non-significant FDRs to NAs so they don't have a dot in the dotplot (this is only important for TEs, which are currently not filtered by significance)
  melted_format_data[which(melted_format_data$p.adjust <= -log10(0.05)), 'p.adjust'] <- NA

  # lock in factor level order (*******IMPORTANT, otherwise GGplot will arrange the y-axis and x-axis by alphabetical order************)
  melted_format_data$ID <- factor(melted_format_data$ID, levels = unique(melted_format_data$ID))
  melted_format_data$Group <- factor(melted_format_data$Group, levels = unique(melted_format_data$Group))
  
  # Remove rows with NAs
  melted_format_data <- na.omit(melted_format_data)


  
  
  
  
  # Return
  return(melted_format_data) 
  
    
    
} # END FUNCTION

Overlap_three_GSEA_res <- function(result.one, result.two, result.three,
                                   label.one, label.two, label.three,
                                   comparison_directionality,
                                   number_to_plot,
                                   dir.output, gs.label){
  
  # NOTE: This function will keep gene sets that are significant in three sets of GSEA results.
  # Something to potentially include: If no overlap or if no consistent overlap, end function immediately and maybe specify message
  # Comparison directionality here specifies the directionalities 1) comparing result 1 vs 2 and 2) comparing result 1 vs 3
  
  
  
  

  
  # -log10 transform the adjusted pvalues
  result.one$p.adjust <- -log10(result.one$p.adjust)
  result.two$p.adjust <- -log10(result.two$p.adjust)
  result.three$p.adjust <- -log10(result.three$p.adjust)

  # Add sample group label column
  result.one$Group <- label.one
  result.two$Group <- label.two
  result.three$Group <- label.three
  
  # Define sig gene sets in each analysis
  a <- result.one$ID
  b <- result.two$ID
  c <- result.three$ID

  # Find gene set overlaps
  common_genesets <- Reduce(intersect, list(a, b, c)
                            )

  # Extract overlapping genesets and columns necessary for plotting
  result.one <- result.one[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.two <- result.two[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.three <- result.three[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  
  # Run function to find gene set IDs with the desired directionality (COMPARING RESULT 1 VS RESULT 2, AND RESULT 1 VS RESULT 3)
  genesets_of_interest_1 <- find_genesets_with_directionality(input.one = result.one, input.two = result.two, input_directionality = comparison_directionality[1])
  genesets_of_interest_2 <- find_genesets_with_directionality(input.one = result.one, input.two = result.three, input_directionality = comparison_directionality[2])
  
  # Define the gene sets that pass all directionality filters
  genesets_of_interest <- Reduce(intersect, list(genesets_of_interest_1, genesets_of_interest_2))
  
  # Keep gene sets changing in the desired way
  result.one <- result.one[genesets_of_interest, ]
  result.two <- result.two[genesets_of_interest, ]
  result.three <- result.three[genesets_of_interest, ]
  
  # Make a subfunction to calculate Fischer meta-analysis pvalues; this will be used on conjunction with 'apply' to calculate a meta pvalue for each row
  calculate.fisher.p <- function(my.pvalues){
    
      # Extract pvalues from the ith row
      my.pvalue.A = my.pvalues[1]
      my.pvalue.B = my.pvalues[2]
      my.pvalue.C = my.pvalues[3]
      
      # Calculate the meta pvalue using Fishers method
      my.meta.pvalue <- poolr::fisher(p = c(my.pvalue.A, my.pvalue.B, my.pvalue.C))
      
      # Output the meta pvalue
      return(my.meta.pvalue$p)
    
  }
  
  # Calculate a combined Fischer pvalue for each gene set
  all.my.meta.pvalues <- apply(cbind(result.one$pvalue, result.two$pvalue, result.three$pvalue), 
                               MARGIN = 1, 
                               calculate.fisher.p)
  
  # Combine the results into a table that will be saved as an output (and that will be different than the table used to plot results)
  my.output.table <- cbind(result.one, 
                           result.two[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.three[, c('Group', 'NES', 'pvalue', 'p.adjust')])
  
  # Add the meta pvalue to the output table
  my.output.table$meta.pvalue <- all.my.meta.pvalues
  
  # Reorder gene sets from biggest to smallest meta pvalue (***this is needed if we want the final plots to show the more significant values at the top)
  my.output.table <- my.output.table[order(my.output.table$meta.pvalue, decreasing = TRUE), ]
  
      # write results to file
      write.table(my.output.table, file = paste(dir.output, 'ClusterProfiler_GSEA_Table_Overlapping_Gene_Sets_', gs.label, '_FDR5_', nrow(my.output.table), '_Significant' ,'.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
  
  # Extract the data for the first comparison, which will be used to define reference positive and negative terms
  result.one <- result.one[rownames(my.output.table),]
  
  
  
  # DEFINE THE NUMBER OF POS AND NEG NES TERMS TO PLOT
  
  # Extract the positive NES gene sets
  positive.genesets <- result.one[which(result.one$NES > 0, ), ]
  
  # Extract the negative NES gene sets
  negative.genesets <- result.one[which(result.one$NES < 0, ), ]
  
  # Define the maximum number of positive or negative terms to plot
  max_to_plot <- number_to_plot / 2
  
  # Define the actual number of positive or negative terms to plot
  if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- max_to_plot
      
      # Define the number of negative terms to plot
      negative_to_plot <- max_to_plot
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the maximum number of negative terms to plot
      max_negative <- number_to_plot - positive_to_plot
      
      # Define the number negative terms to plot, either the maximum possible or the number available (whichever is less)
      negative_to_plot <- min(max_negative, nrow(negative.genesets))
      
  } else if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
      
      # Define the maximum number of positive terms to plot
      max_positive <- number_to_plot - negative_to_plot
      
      # Define the number positive terms to plot, either the maximum possible or the number available (whichever is less)
      positive_to_plot <- min(max_positive, nrow(positive.genesets))
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
  }
  
  # Define empty vectors of gene sets to plot, to update if there are terms to plot
  pos.genesets.to.plot <- c()
  neg.genesets.to.plot <- c()
  
  
  
  
  
  # IDENTIFY GENE SET NAMES TO PLOT
  
  if (positive_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      pos_starting_index <- nrow(positive.genesets) - positive_to_plot + 1
      
      # Define the index for the last geneset to plot
      pos_ending_index <- nrow(positive.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      pos.genesets.to.plot <- positive.genesets[pos_starting_index:pos_ending_index, 'ID']
      
  } 
  
  if (negative_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      neg_starting_index <- nrow(negative.genesets) - negative_to_plot + 1
      
      # Define the index for the last geneset to plot
      neg_ending_index <- nrow(negative.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      neg.genesets.to.plot <- negative.genesets[neg_starting_index:neg_ending_index, 'ID']
      
  } 
  
  
  
  
  
  # Collect the positive and negative gene set names to plot
  genesets.to.plot <- c(neg.genesets.to.plot, pos.genesets.to.plot)
      
  # Combine results into 1 dataframe in the *melted* format ggplot needs
  melted_format_data <- rbind(result.one[genesets.to.plot, ], 
                              result.two[genesets.to.plot, ], 
                              result.three[genesets.to.plot, ])
  
  # Change non-significant FDRs to NAs so they don't have a dot in the dotplot (this is only important for TEs, which are currently not filtered by significance)
  melted_format_data[which(melted_format_data$p.adjust <= -log10(0.05)), 'p.adjust'] <- NA

  # lock in factor level order (*******IMPORTANT, otherwise GGplot will arrange the y-axis and x-axis by alphabetical order************)
  melted_format_data$ID <- factor(melted_format_data$ID, levels = unique(melted_format_data$ID))
  melted_format_data$Group <- factor(melted_format_data$Group, levels = unique(melted_format_data$Group))
  
  # Remove rows with NAs
  melted_format_data <- na.omit(melted_format_data)


  
  
  
  #return(melted_format_data) 
  return(melted_format_data) 
  
    
    
} # END FUNCTION

Overlap_three_GSEA_res_unfiltered <- function(result.one, result.two, result.three,
                                             label.one, label.two, label.three,
                                             comparison_directionality,
                                             number_to_plot,
                                             dir.output, gs.label){
            
            # NOTE: This function will keep gene sets that are present in three sets of GSEA results.
            # Something to potentially include: If no overlap or if no consistent overlap, end function immediately and maybe specify message
            # Comparison directionality here specifies the directionalities 1) comparing result 1 vs 2 and 2) comparing result 1 vs 3
            
            
            
            
          
            
            # -log10 transform the adjusted pvalues
            result.one$p.adjust <- -log10(result.one$p.adjust)
            result.two$p.adjust <- -log10(result.two$p.adjust)
            result.three$p.adjust <- -log10(result.three$p.adjust)
          
            # Add sample group label column
            result.one$Group <- label.one
            result.two$Group <- label.two
            result.three$Group <- label.three
            
            # Define sig gene sets in each analysis
            a <- result.one$ID
            b <- result.two$ID
            c <- result.three$ID
          
            # Find gene set overlaps
            common_genesets <- Reduce(intersect, list(a, b, c)
                                      )
          
            # Extract overlapping genesets and columns necessary for plotting
            result.one <- result.one[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
            result.two <- result.two[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
            result.three <- result.three[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
            
            # Run function to find gene set IDs with the desired directionality (COMPARING RESULT 1 VS RESULT 2, AND RESULT 1 VS RESULT 3)
            genesets_of_interest_1 <- find_genesets_with_directionality(input.one = result.one, input.two = result.two, input_directionality = comparison_directionality[1])
            genesets_of_interest_2 <- find_genesets_with_directionality(input.one = result.one, input.two = result.three, input_directionality = comparison_directionality[2])
            
            # Define the gene sets that pass all directionality filters
            genesets_of_interest <- Reduce(intersect, list(genesets_of_interest_1, genesets_of_interest_2))
            
            # Keep gene sets changing in the desired way
            result.one <- result.one[genesets_of_interest, ]
            result.two <- result.two[genesets_of_interest, ]
            result.three <- result.three[genesets_of_interest, ]
            
            # Make a subfunction to calculate Fischer meta-analysis pvalues; this will be used on conjunction with 'apply' to calculate a meta pvalue for each row
            calculate.fisher.p <- function(my.pvalues){
              
                # Extract pvalues from the ith row
                my.pvalue.A = my.pvalues[1]
                my.pvalue.B = my.pvalues[2]
                my.pvalue.C = my.pvalues[3]
                
                # Calculate the meta pvalue using Fishers method
                my.meta.pvalue <- poolr::fisher(p = c(my.pvalue.A, my.pvalue.B, my.pvalue.C))
                
                # Output the meta pvalue
                return(my.meta.pvalue$p)
              
            }
            
            # Calculate a combined Fischer pvalue for each gene set
            all.my.meta.pvalues <- apply(cbind(result.one$pvalue, result.two$pvalue, result.three$pvalue), 
                                         MARGIN = 1, 
                                         calculate.fisher.p)
            
            # Combine the results into a table that will be saved as an output (and that will be different than the table used to plot results)
            my.output.table <- cbind(result.one, 
                                     result.two[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                                     result.three[, c('Group', 'NES', 'pvalue', 'p.adjust')])
            
            # Add the meta pvalue to the output table
            my.output.table$meta.pvalue <- all.my.meta.pvalues
            
            # Reorder gene sets from biggest to smallest meta pvalue (***this is needed if we want the final plots to show the more significant values at the top)
            my.output.table <- my.output.table[order(my.output.table$meta.pvalue, decreasing = TRUE), ]
            
                # write results to file
                write.table(my.output.table, file = paste(dir.output, 'ClusterProfiler_GSEA_Table_Overlapping_Gene_Sets_', gs.label, '_FDR100_', nrow(my.output.table), '_Analyzed' ,'.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
            
            # Extract the data for the first comparison, which will be used to define reference positive and negative terms
            result.one <- result.one[rownames(my.output.table),]
            
            
            
            # DEFINE THE NUMBER OF POS AND NEG NES TERMS TO PLOT
            
            # Extract the positive NES gene sets
            positive.genesets <- result.one[which(result.one$NES > 0, ), ]
            
            # Extract the negative NES gene sets
            negative.genesets <- result.one[which(result.one$NES < 0, ), ]
            
            # Define the maximum number of positive or negative terms to plot
            max_to_plot <- number_to_plot / 2
            
            # Define the actual number of positive or negative terms to plot
            if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) >= max_to_plot) {
              
                # Define the number of positive terms to plot
                positive_to_plot <- max_to_plot
                
                # Define the number of negative terms to plot
                negative_to_plot <- max_to_plot
                
            } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) >= max_to_plot) {
              
                # Define the number of positive terms to plot
                positive_to_plot <- nrow(positive.genesets)
                
                # Define the maximum number of negative terms to plot
                max_negative <- number_to_plot - positive_to_plot
                
                # Define the number negative terms to plot, either the maximum possible or the number available (whichever is less)
                negative_to_plot <- min(max_negative, nrow(negative.genesets))
                
            } else if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) < max_to_plot) {
              
                # Define the number of negative terms to plot
                negative_to_plot <- nrow(negative.genesets)
                
                # Define the maximum number of positive terms to plot
                max_positive <- number_to_plot - negative_to_plot
                
                # Define the number positive terms to plot, either the maximum possible or the number available (whichever is less)
                positive_to_plot <- min(max_positive, nrow(positive.genesets))
                
            } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) < max_to_plot) {
              
                # Define the number of positive terms to plot
                positive_to_plot <- nrow(positive.genesets)
                
                # Define the number of negative terms to plot
                negative_to_plot <- nrow(negative.genesets)
            }
            
            # Define empty vectors of gene sets to plot, to update if there are terms to plot
            pos.genesets.to.plot <- c()
            neg.genesets.to.plot <- c()
            
            
            
            
            
            # IDENTIFY GENE SET NAMES TO PLOT
            
            if (positive_to_plot > 0) {
              
                # Define the index for the first geneset to plot
                pos_starting_index <- nrow(positive.genesets) - positive_to_plot + 1
                
                # Define the index for the last geneset to plot
                pos_ending_index <- nrow(positive.genesets) 
                
                # Extract the top N gene sets (remember, most significant are at the bottom)
                pos.genesets.to.plot <- positive.genesets[pos_starting_index:pos_ending_index, 'ID']
                
            } 
            
            if (negative_to_plot > 0) {
              
                # Define the index for the first geneset to plot
                neg_starting_index <- nrow(negative.genesets) - negative_to_plot + 1
                
                # Define the index for the last geneset to plot
                neg_ending_index <- nrow(negative.genesets) 
                
                # Extract the top N gene sets (remember, most significant are at the bottom)
                neg.genesets.to.plot <- negative.genesets[neg_starting_index:neg_ending_index, 'ID']
                
            } 
            
            
            
            
            
            # Collect the positive and negative gene set names to plot
            genesets.to.plot <- c(neg.genesets.to.plot, pos.genesets.to.plot)
                
            # Combine results into 1 dataframe in the *melted* format ggplot needs
            melted_format_data <- rbind(result.one[genesets.to.plot, ], 
                                        result.two[genesets.to.plot, ], 
                                        result.three[genesets.to.plot, ])
            
            # Change non-significant FDRs to NAs so they don't have a dot in the dotplot (this is only important for TEs, which are currently not filtered by significance)
            #melted_format_data[which(melted_format_data$p.adjust <= -log10(0.05)), 'p.adjust'] <- NA
          
            # lock in factor level order (*******IMPORTANT, otherwise GGplot will arrange the y-axis and x-axis by alphabetical order************)
            melted_format_data$ID <- factor(melted_format_data$ID, levels = unique(melted_format_data$ID))
            melted_format_data$Group <- factor(melted_format_data$Group, levels = unique(melted_format_data$Group))
            
            # Remove rows with NAs
            melted_format_data <- na.omit(melted_format_data)
          
          
            
            
            
            #return(melted_format_data) 
            return(melted_format_data) 
            
              
              
          } # END FUNCTION

Overlap_four_GSEA_res <- function(result.one, result.two, result.three, result.four,
                                  label.one, label.two, label.three, label.four,
                                  comparison_directionality,
                                  number_to_plot,
                                  dir.output, gs.label){
  
  # NOTE: This function will keep gene sets that are significant in four sets of GSEA results.
  # Something to potentially include: If no overlap or if no consistent overlap, end function immediately and maybe specify message
  # Comparison directionality here specifies the directionalities 1) comparing result 1 vs 2, 2) comparing result 1 vs 3, 3) comparing result 1 vs 4, etc......
  
  
  
  

  
  # -log10 transform the adjusted pvalues
  result.one$p.adjust <- -log10(result.one$p.adjust)
  result.two$p.adjust <- -log10(result.two$p.adjust)
  result.three$p.adjust <- -log10(result.three$p.adjust)
  result.four$p.adjust <- -log10(result.four$p.adjust)

  # Add sample group label column
  result.one$Group <- label.one
  result.two$Group <- label.two
  result.three$Group <- label.three
  result.four$Group <- label.four
  
  # Define sig gene sets in each analysis
  a <- result.one$ID
  b <- result.two$ID
  c <- result.three$ID
  d <- result.four$ID
  
  # Find gene set overlaps
  common_genesets <- Reduce(intersect, list(a, b, c, d))

  # Extract overlapping genesets and columns necessary for plotting
  result.one <- result.one[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.two <- result.two[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.three <- result.three[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.four <- result.four[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  
  # Run function to find gene set IDs with the desired directionality (COMPARING RESULT 1 VS RESULT 2, AND RESULT 1 VS RESULT 3)
  genesets_of_interest_1 <- find_genesets_with_directionality(input.one = result.one, input.two = result.two, input_directionality = comparison_directionality[1])
  genesets_of_interest_2 <- find_genesets_with_directionality(input.one = result.one, input.two = result.three, input_directionality = comparison_directionality[2])
  genesets_of_interest_3 <- find_genesets_with_directionality(input.one = result.one, input.two = result.four, input_directionality = comparison_directionality[3])
  
  # Define the gene sets that pass all directionality filters
  genesets_of_interest <- Reduce(intersect, list(genesets_of_interest_1, genesets_of_interest_2, genesets_of_interest_3))
  
  # Keep gene sets changing in the desired way
  result.one <- result.one[genesets_of_interest, ]
  result.two <- result.two[genesets_of_interest, ]
  result.three <- result.three[genesets_of_interest, ]
  result.four <- result.four[genesets_of_interest, ]
  
  # Make a subfunction to calculate Fischer meta-analysis pvalues; this will be used on conjunction with 'apply' to calculate a meta pvalue for each row
  calculate.fisher.p <- function(my.pvalues){
    
      # Extract pvalues from the ith row
      my.pvalue.A = my.pvalues[1]
      my.pvalue.B = my.pvalues[2]
      my.pvalue.C = my.pvalues[3]
      my.pvalue.D = my.pvalues[4]
      
      # Calculate the meta pvalue using Fishers method
      my.meta.pvalue <- poolr::fisher(p = c(my.pvalue.A, my.pvalue.B, my.pvalue.C, my.pvalue.D))
      
      # Output the meta pvalue
      return(my.meta.pvalue$p)
    
  }
  
  # Calculate a combined Fischer pvalue for each gene set
  all.my.meta.pvalues <- apply(cbind(result.one$pvalue, result.two$pvalue, result.three$pvalue, result.four$pvalue), 
                               MARGIN = 1, 
                               calculate.fisher.p)
  
  # Combine the results into a table that will be saved as an output (and that will be different than the table used to plot results)
  my.output.table <- cbind(result.one, 
                           result.two[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.three[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.four[, c('Group', 'NES', 'pvalue', 'p.adjust')])
  
  # Add the meta pvalue to the output table
  my.output.table$meta.pvalue <- all.my.meta.pvalues
  
  # Reorder gene sets from biggest to smallest meta pvalue (***this is needed if we want the final plots to show the more significant values at the top)
  my.output.table <- my.output.table[order(my.output.table$meta.pvalue, decreasing = TRUE), ]
    
      # write results to file
      write.table(my.output.table, file = paste(dir.output, 'ClusterProfiler_GSEA_Table_Overlapping_Gene_Sets_', gs.label, '_FDR5_', nrow(my.output.table), '_Significant' ,'.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
  
  # Extract the data for the first comparison, which will be used to define reference positive and negative terms
  result.one <- result.one[rownames(my.output.table),]
  
  
  
  # DEFINE THE NUMBER OF POS AND NEG NES TERMS TO PLOT
  
  # Extract the positive NES gene sets
  positive.genesets <- result.one[which(result.one$NES > 0, ), ]
  
  # Extract the negative NES gene sets
  negative.genesets <- result.one[which(result.one$NES < 0, ), ]
  
  # Define the maximum number of positive or negative terms to plot
  max_to_plot <- number_to_plot / 2
  
  # Define the actual number of positive or negative terms to plot
  if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- max_to_plot
      
      # Define the number of negative terms to plot
      negative_to_plot <- max_to_plot
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the maximum number of negative terms to plot
      max_negative <- number_to_plot - positive_to_plot
      
      # Define the number negative terms to plot, either the maximum possible or the number available (whichever is less)
      negative_to_plot <- min(max_negative, nrow(negative.genesets))
      
  } else if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
      
      # Define the maximum number of positive terms to plot
      max_positive <- number_to_plot - negative_to_plot
      
      # Define the number positive terms to plot, either the maximum possible or the number available (whichever is less)
      positive_to_plot <- min(max_positive, nrow(positive.genesets))
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
  }
  
  # Define empty vectors of gene sets to plot, to update if there are terms to plot
  pos.genesets.to.plot <- c()
  neg.genesets.to.plot <- c()
  
  
  
  
  
  # IDENTIFY GENE SET NAMES TO PLOT
  
  if (positive_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      pos_starting_index <- nrow(positive.genesets) - positive_to_plot + 1
      
      # Define the index for the last geneset to plot
      pos_ending_index <- nrow(positive.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      pos.genesets.to.plot <- positive.genesets[pos_starting_index:pos_ending_index, 'ID']
      
  } 
  
  if (negative_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      neg_starting_index <- nrow(negative.genesets) - negative_to_plot + 1
      
      # Define the index for the last geneset to plot
      neg_ending_index <- nrow(negative.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      neg.genesets.to.plot <- negative.genesets[neg_starting_index:neg_ending_index, 'ID']
      
  } 
  
  
  
  
  
  # Collect the positive and negative gene set names to plot
  genesets.to.plot <- c(neg.genesets.to.plot, pos.genesets.to.plot)
      
  # Combine results into 1 dataframe in the *melted* format ggplot needs
  melted_format_data <- rbind(result.one[genesets.to.plot, ], 
                              result.two[genesets.to.plot, ], 
                              result.three[genesets.to.plot, ],
                              result.four[genesets.to.plot, ])
  
  # Change non-significant FDRs to NAs so they don't have a dot in the dotplot (this is only important for TEs, which are currently not filtered by significance)
  melted_format_data[which(melted_format_data$p.adjust <= -log10(0.05)), 'p.adjust'] <- NA

  # lock in factor level order (*******IMPORTANT, otherwise GGplot will arrange the y-axis and x-axis by alphabetical order************)
  melted_format_data$ID <- factor(melted_format_data$ID, levels = unique(melted_format_data$ID))
  melted_format_data$Group <- factor(melted_format_data$Group, levels = unique(melted_format_data$Group))
  
  # Remove rows with NAs
  melted_format_data <- na.omit(melted_format_data)


  
  
  # Return
  return(melted_format_data) 
  
    
    
} # END FUNCTION

Overlap_five_GSEA_res <- function(result.one, result.two, result.three, result.four, result.five, 
                                  label.one, label.two, label.three, label.four, label.five,
                                  comparison_directionality,
                                  number_to_plot,
                                  dir.output, gs.label){
  
  # NOTE: This function will keep gene sets that are significant in five sets of GSEA results.
  # Something to potentially include: If no overlap or if no consistent overlap, end function immediately and maybe specify message
  # Comparison directionality here specifies the directionalities 1) comparing result 1 vs 2, 2) comparing result 1 vs 3, 3) comparing result 1 vs 4, etc......
  
  
  
  

  
  # -log10 transform the adjusted pvalues
  result.one$p.adjust <- -log10(result.one$p.adjust)
  result.two$p.adjust <- -log10(result.two$p.adjust)
  result.three$p.adjust <- -log10(result.three$p.adjust)
  result.four$p.adjust <- -log10(result.four$p.adjust)
  result.five$p.adjust <- -log10(result.five$p.adjust)

  # Add sample group label column
  result.one$Group <- label.one
  result.two$Group <- label.two
  result.three$Group <- label.three
  result.four$Group <- label.four
  result.five$Group <- label.five

  # Define sig gene sets in each analysis
  a <- result.one$ID
  b <- result.two$ID
  c <- result.three$ID
  d <- result.four$ID
  e <- result.five$ID
  
  # Find gene set overlaps
  common_genesets <- Reduce(intersect, list(a, b, c, d, e))

  # Extract overlapping genesets and columns necessary for plotting
  result.one <- result.one[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.two <- result.two[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.three <- result.three[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.four <- result.four[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.five <- result.five[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
 
  # Run function to find gene set IDs with the desired directionality (COMPARING RESULT 1 VS RESULT 2, AND RESULT 1 VS RESULT 3)
  genesets_of_interest_1 <- find_genesets_with_directionality(input.one = result.one, input.two = result.two, input_directionality = comparison_directionality[1])
  genesets_of_interest_2 <- find_genesets_with_directionality(input.one = result.one, input.two = result.three, input_directionality = comparison_directionality[2])
  genesets_of_interest_3 <- find_genesets_with_directionality(input.one = result.one, input.two = result.four, input_directionality = comparison_directionality[3])
  genesets_of_interest_4 <- find_genesets_with_directionality(input.one = result.one, input.two = result.five, input_directionality = comparison_directionality[4])
  
  # Define the gene sets that pass all directionality filters
  genesets_of_interest <- Reduce(intersect, list(genesets_of_interest_1, genesets_of_interest_2, genesets_of_interest_3, genesets_of_interest_4))
  
  # Keep gene sets changing in the desired way
  result.one <- result.one[genesets_of_interest, ]
  result.two <- result.two[genesets_of_interest, ]
  result.three <- result.three[genesets_of_interest, ]
  result.four <- result.four[genesets_of_interest, ]
  result.five <- result.five[genesets_of_interest, ]
  
  # Make a subfunction to calculate Fischer meta-analysis pvalues; this will be used on conjunction with 'apply' to calculate a meta pvalue for each row
  calculate.fisher.p <- function(my.pvalues){
    
      # Extract pvalues from the ith row
      my.pvalue.A = my.pvalues[1]
      my.pvalue.B = my.pvalues[2]
      my.pvalue.C = my.pvalues[3]
      my.pvalue.D = my.pvalues[4]
      my.pvalue.E = my.pvalues[5]
      
      # Calculate the meta pvalue using Fishers method
      my.meta.pvalue <- poolr::fisher(p = c(my.pvalue.A, my.pvalue.B, my.pvalue.C, my.pvalue.D, my.pvalue.E))
      
      # Output the meta pvalue
      return(my.meta.pvalue$p)
    
  }
  
  # Calculate a combined Fischer pvalue for each gene set
  all.my.meta.pvalues <- apply(cbind(result.one$pvalue, result.two$pvalue, result.three$pvalue, result.four$pvalue, result.five$pvalue), 
                               MARGIN = 1, 
                               calculate.fisher.p)
  
  # Combine the results into a table that will be saved as an output (and that will be different than the table used to plot results)
  my.output.table <- cbind(result.one, 
                           result.two[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.three[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.four[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.five[, c('Group', 'NES', 'pvalue', 'p.adjust')])
  
  # Add the meta pvalue to the output table
  my.output.table$meta.pvalue <- all.my.meta.pvalues
  
  # Reorder gene sets from biggest to smallest meta pvalue (***this is needed if we want the final plots to show the more significant values at the top)
  my.output.table <- my.output.table[order(my.output.table$meta.pvalue, decreasing = TRUE), ]
  
      # write results to file
      write.table(my.output.table, file = paste(dir.output, 'ClusterProfiler_GSEA_Table_Overlapping_Gene_Sets_', gs.label, '_FDR5_', nrow(my.output.table), '_Significant' ,'.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

  # Extract the data for the first comparison, which will be used to define reference positive and negative terms
  result.one <- result.one[rownames(my.output.table),]
  
  
  
  # DEFINE THE NUMBER OF POS AND NEG NES TERMS TO PLOT
  
  # Extract the positive NES gene sets
  positive.genesets <- result.one[which(result.one$NES > 0, ), ]
  
  # Extract the negative NES gene sets
  negative.genesets <- result.one[which(result.one$NES < 0, ), ]
  
  # Define the maximum number of positive or negative terms to plot
  max_to_plot <- number_to_plot / 2
  
  # Define the actual number of positive or negative terms to plot
  if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- max_to_plot
      
      # Define the number of negative terms to plot
      negative_to_plot <- max_to_plot
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the maximum number of negative terms to plot
      max_negative <- number_to_plot - positive_to_plot
      
      # Define the number negative terms to plot, either the maximum possible or the number available (whichever is less)
      negative_to_plot <- min(max_negative, nrow(negative.genesets))
      
  } else if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
      
      # Define the maximum number of positive terms to plot
      max_positive <- number_to_plot - negative_to_plot
      
      # Define the number positive terms to plot, either the maximum possible or the number available (whichever is less)
      positive_to_plot <- min(max_positive, nrow(positive.genesets))
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
  }
  
  # Define empty vectors of gene sets to plot, to update if there are terms to plot
  pos.genesets.to.plot <- c()
  neg.genesets.to.plot <- c()
  
  
  
  
  
  # IDENTIFY GENE SET NAMES TO PLOT
  
  if (positive_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      pos_starting_index <- nrow(positive.genesets) - positive_to_plot + 1
      
      # Define the index for the last geneset to plot
      pos_ending_index <- nrow(positive.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      pos.genesets.to.plot <- positive.genesets[pos_starting_index:pos_ending_index, 'ID']
      
  } 
  
  if (negative_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      neg_starting_index <- nrow(negative.genesets) - negative_to_plot + 1
      
      # Define the index for the last geneset to plot
      neg_ending_index <- nrow(negative.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      neg.genesets.to.plot <- negative.genesets[neg_starting_index:neg_ending_index, 'ID']
      
  } 
  
  
  
  
  
  # Collect the positive and negative gene set names to plot
  genesets.to.plot <- c(neg.genesets.to.plot, pos.genesets.to.plot)
      
  # Combine results into 1 dataframe in the *melted* format ggplot needs
  melted_format_data <- rbind(result.one[genesets.to.plot, ], 
                              result.two[genesets.to.plot, ], 
                              result.three[genesets.to.plot, ],
                              result.four[genesets.to.plot, ],
                              result.five[genesets.to.plot, ])
  
  # Change non-significant FDRs to NAs so they don't have a dot in the dotplot (this is only important for TEs, which are currently not filtered by significance)
  melted_format_data[which(melted_format_data$p.adjust <= -log10(0.05)), 'p.adjust'] <- NA

  # lock in factor level order (*******IMPORTANT, otherwise GGplot will arrange the y-axis and x-axis by alphabetical order************)
  melted_format_data$ID <- factor(melted_format_data$ID, levels = unique(melted_format_data$ID))
  melted_format_data$Group <- factor(melted_format_data$Group, levels = unique(melted_format_data$Group))

  # Remove rows with NAs
  melted_format_data <- na.omit(melted_format_data)

  
  
  
  # Return
  return(melted_format_data) 
  
    
    
} # END FUNCTION

Overlap_seven_GSEA_res <- function(result.one, result.two, result.three, result.four, result.five, result.six, result.seven,
                                   label.one, label.two, label.three, label.four, label.five, label.six, label.seven,
                                   comparison_directionality,
                                   number_to_plot,
                                   dir.output, gs.label){
  
  # NOTE: This function will keep gene sets that are significant in five sets of GSEA results.
  # Something to potentially include: If no overlap or if no consistent overlap, end function immediately and maybe specify message
  # Comparison directionality here specifies the directionalities 1) comparing result 1 vs 2, 2) comparing result 1 vs 3, 3) comparing result 1 vs 4, etc......
  
  
  
  

  
  # -log10 transform the adjusted pvalues
  result.one$p.adjust <- -log10(result.one$p.adjust)
  result.two$p.adjust <- -log10(result.two$p.adjust)
  result.three$p.adjust <- -log10(result.three$p.adjust)
  result.four$p.adjust <- -log10(result.four$p.adjust)
  result.five$p.adjust <- -log10(result.five$p.adjust)
  result.six$p.adjust <- -log10(result.six$p.adjust)
  result.seven$p.adjust <- -log10(result.seven$p.adjust)

  # Add sample group label column
  result.one$Group <- label.one
  result.two$Group <- label.two
  result.three$Group <- label.three
  result.four$Group <- label.four
  result.five$Group <- label.five
  result.six$Group <- label.six
  result.seven$Group <- label.seven
  
  # Define sig gene sets in each analysis
  a <- result.one$ID
  b <- result.two$ID
  c <- result.three$ID
  d <- result.four$ID
  e <- result.five$ID
  f <- result.six$ID
  g <- result.seven$ID
  
  # Find gene set overlaps
  common_genesets <- Reduce(intersect, list(a, b, c, d, e, f, g))

  # Extract overlapping genesets and columns necessary for plotting
  result.one <- result.one[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.two <- result.two[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.three <- result.three[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.four <- result.four[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.five <- result.five[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.six <- result.six[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  result.seven <- result.seven[common_genesets, c('ID', 'Group', 'NES', 'pvalue', 'p.adjust')]
  
  # Run function to find gene set IDs with the desired directionality (COMPARING RESULT 1 VS RESULT 2, AND RESULT 1 VS RESULT 3)
  genesets_of_interest_1 <- find_genesets_with_directionality(input.one = result.one, input.two = result.two, input_directionality = comparison_directionality[1])
  genesets_of_interest_2 <- find_genesets_with_directionality(input.one = result.one, input.two = result.three, input_directionality = comparison_directionality[2])
  genesets_of_interest_3 <- find_genesets_with_directionality(input.one = result.one, input.two = result.four, input_directionality = comparison_directionality[3])
  genesets_of_interest_4 <- find_genesets_with_directionality(input.one = result.one, input.two = result.five, input_directionality = comparison_directionality[4])
  genesets_of_interest_5 <- find_genesets_with_directionality(input.one = result.one, input.two = result.six, input_directionality = comparison_directionality[5])
  genesets_of_interest_6 <- find_genesets_with_directionality(input.one = result.one, input.two = result.seven, input_directionality = comparison_directionality[6])
  
  # Define the gene sets that pass all directionality filters
  genesets_of_interest <- Reduce(intersect, list(genesets_of_interest_1, genesets_of_interest_2, genesets_of_interest_3, genesets_of_interest_4, genesets_of_interest_5, genesets_of_interest_6))
  
  # Keep gene sets changing in the desired way
  result.one <- result.one[genesets_of_interest, ]
  result.two <- result.two[genesets_of_interest, ]
  result.three <- result.three[genesets_of_interest, ]
  result.four <- result.four[genesets_of_interest, ]
  result.five <- result.five[genesets_of_interest, ]
  result.six <- result.six[genesets_of_interest, ]
  result.seven <- result.seven[genesets_of_interest, ]
  
  # Make a subfunction to calculate Fischer meta-analysis pvalues; this will be used on conjunction with 'apply' to calculate a meta pvalue for each row
  calculate.fisher.p <- function(my.pvalues){
    
      # Extract pvalues from the ith row
      my.pvalue.A = my.pvalues[1]
      my.pvalue.B = my.pvalues[2]
      my.pvalue.C = my.pvalues[3]
      my.pvalue.D = my.pvalues[4]
      my.pvalue.E = my.pvalues[5]
      my.pvalue.F = my.pvalues[6]
      my.pvalue.G = my.pvalues[7]
      
      # Calculate the meta pvalue using Fishers method
      my.meta.pvalue <- poolr::fisher(p = c(my.pvalue.A, my.pvalue.B, my.pvalue.C, my.pvalue.D, my.pvalue.E, my.pvalue.F, my.pvalue.G))
      
      # Output the meta pvalue
      return(my.meta.pvalue$p)
    
  }
  
  # Calculate a combined Fischer pvalue for each gene set
  all.my.meta.pvalues <- apply(cbind(result.one$pvalue, result.two$pvalue, result.three$pvalue, result.four$pvalue, result.five$pvalue, result.six$pvalue, result.seven$pvalue), 
                               MARGIN = 1, 
                               calculate.fisher.p)
  
  # Combine the results into a table that will be saved as an output (and that will be different than the table used to plot results)
  my.output.table <- cbind(result.one, 
                           result.two[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.three[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.four[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.five[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.six[, c('Group', 'NES', 'pvalue', 'p.adjust')], 
                           result.seven[, c('Group', 'NES', 'pvalue', 'p.adjust')])
  
  # Add the meta pvalue to the output table
  my.output.table$meta.pvalue <- all.my.meta.pvalues
  
  # Reorder gene sets from biggest to smallest meta pvalue (***this is needed if we want the final plots to show the more significant values at the top)
  my.output.table <- my.output.table[order(my.output.table$meta.pvalue, decreasing = TRUE), ]
  
      # write results to file
      write.table(my.output.table, file = paste(dir.output, 'ClusterProfiler_GSEA_Table_Overlapping_Gene_Sets_', gs.label, '_FDR5_', nrow(my.output.table), '_Significant' ,'.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
  
  # Extract the data for the first comparison, which will be used to define reference positive and negative terms
  result.one <- result.one[rownames(my.output.table),]
  
  
  
  # DEFINE THE NUMBER OF POS AND NEG NES TERMS TO PLOT
  
  # Extract the positive NES gene sets
  positive.genesets <- result.one[which(result.one$NES > 0, ), ]
  
  # Extract the negative NES gene sets
  negative.genesets <- result.one[which(result.one$NES < 0, ), ]
  
  # Define the maximum number of positive or negative terms to plot
  max_to_plot <- number_to_plot / 2
  
  # Define the actual number of positive or negative terms to plot
  if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- max_to_plot
      
      # Define the number of negative terms to plot
      negative_to_plot <- max_to_plot
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) >= max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the maximum number of negative terms to plot
      max_negative <- number_to_plot - positive_to_plot
      
      # Define the number negative terms to plot, either the maximum possible or the number available (whichever is less)
      negative_to_plot <- min(max_negative, nrow(negative.genesets))
      
  } else if (nrow(positive.genesets) >= max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
      
      # Define the maximum number of positive terms to plot
      max_positive <- number_to_plot - negative_to_plot
      
      # Define the number positive terms to plot, either the maximum possible or the number available (whichever is less)
      positive_to_plot <- min(max_positive, nrow(positive.genesets))
      
  } else if (nrow(positive.genesets) < max_to_plot & nrow(negative.genesets) < max_to_plot) {
    
      # Define the number of positive terms to plot
      positive_to_plot <- nrow(positive.genesets)
      
      # Define the number of negative terms to plot
      negative_to_plot <- nrow(negative.genesets)
  }
  
  # Define empty vectors of gene sets to plot, to update if there are terms to plot
  pos.genesets.to.plot <- c()
  neg.genesets.to.plot <- c()
  
  
  
  
  
  # IDENTIFY GENE SET NAMES TO PLOT
  
  if (positive_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      pos_starting_index <- nrow(positive.genesets) - positive_to_plot + 1
      
      # Define the index for the last geneset to plot
      pos_ending_index <- nrow(positive.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      pos.genesets.to.plot <- positive.genesets[pos_starting_index:pos_ending_index, 'ID']
      
  } 
  
  if (negative_to_plot > 0) {
    
      # Define the index for the first geneset to plot
      neg_starting_index <- nrow(negative.genesets) - negative_to_plot + 1
      
      # Define the index for the last geneset to plot
      neg_ending_index <- nrow(negative.genesets) 
      
      # Extract the top N gene sets (remember, most significant are at the bottom)
      neg.genesets.to.plot <- negative.genesets[neg_starting_index:neg_ending_index, 'ID']
      
  } 
  
  
  
  
  
  # Collect the positive and negative gene set names to plot
  genesets.to.plot <- c(neg.genesets.to.plot, pos.genesets.to.plot)
      
  # Combine results into 1 dataframe in the *melted* format ggplot needs
  melted_format_data <- rbind(result.one[genesets.to.plot, ], 
                              result.two[genesets.to.plot, ], 
                              result.three[genesets.to.plot, ],
                              result.four[genesets.to.plot, ],
                              result.five[genesets.to.plot, ],
                              result.six[genesets.to.plot, ],
                              result.seven[genesets.to.plot, ])
  
  # Change non-significant FDRs to NAs so they don't have a dot in the dotplot (this is only important for TEs, which are currently not filtered by significance)
  melted_format_data[which(melted_format_data$p.adjust <= -log10(0.05)), 'p.adjust'] <- NA

  # lock in factor level order (*******IMPORTANT, otherwise GGplot will arrange the y-axis and x-axis by alphabetical order************)
  melted_format_data$ID <- factor(melted_format_data$ID, levels = unique(melted_format_data$ID))
  melted_format_data$Group <- factor(melted_format_data$Group, levels = unique(melted_format_data$Group))
  
  # Remove rows with NAs
  melted_format_data <- na.omit(melted_format_data)

  
  
  
  # Return
  return(melted_format_data) 
  
    
    
} # END FUNCTION
