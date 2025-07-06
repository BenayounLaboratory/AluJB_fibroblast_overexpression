README -- Analysis of my AluJb overexpression mass spectrometry data
############################

  1. Take protein abundances and transform to log10 or VSN scale and carry out MDS analysis for both cell proteomics and secretomics results (Process_protein_abundances.R)

  2. Generate differential protein expression heatmaps (Generate_Heatmaps.R) and Venn diagrams + mosaic plots comparing the proteomics/secretomics results with aging/senescence gene lists (Generate_Venn_Diagrams_and_Mosaic_Plots.R)

  3. Run GSEA enrichment analysis for cell proteomics/secretomics:
      - Run GSEA and plot the top 5 gene sets in each direction (Run_GSEA_Cellular_and_Secreted_Proteomics_v3.R)
          - Note: The gene sets prepared for the GSEA analysis in '2_RNASeq_Analysis' are used here
          - Note: The functions ('Run_GSEA_Functions.R') prepared for the GSEA analysis in '2_RNASeq_Analysis' are used here
      - Plot the top overlapping gene sets across the transcriptome/proteome/secretome (Run_Overlap_Multiome_GSEA_Results.R)

  4. Run mitch to identify common dysregulated gene sets across the transcriptome/proteome/secretome (Run_Mitch_on_Alu_Multiomics.R)
      - Generate boxplots for differentially altered genes/proteins across all three -omic analyses (Generate_boxplots_common_genes.R)
