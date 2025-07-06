README -- Analysis of my AluJb overexpression transcriptomic data
############################

  1. Trim reads with fastp, check quality with FastQC, and map to hg38 with STAR (Trim_and_Map_AluJb_OE.sh)

  2. Count reads mapping to genes or repeats using TETranscripts or TElocal (TETranscripts_AluJb_OE.sh)

  3. Run DESeq2 analysis for DE genes/TEs as a function of Alu overexpression levels
      - Filter lowly expressed genes (Run_Filter_Lowly_Expressed_Genes_v2.R, Process_RNASeq_functions.R)
      - Calculate fraction of reads from repeats, remove batch effects, run DESeq2, VST transform counts, and generate MDS/PCA plots (Run_DESeq2.R, Process_RNASeq_functions.R)

  4. Generate differential gene expression heatmaps (Generate_Heatmaps.R) and Venn diagrams + mosaic plots comparing the transcriptomic data with aging/senescence gene lists (Generate_Venn_Diagrams_and_Mosaic_Plots.R)

  5. Run GSEA enrichment analysis:
      - Prepare gene set collections for GSEA (Prepare_Gene_Sets.R)
      - Run GSEA and plot the top 5 gene sets in each direction (Run_GSEA_Transcriptomics.R, Run_GSEA_Functions.R)
      - Find overlapping gene sets across transcriptomic analyses (Run_Comparison_of_GSEA_Results.R, Run_GSEA_Functions.R)

  6. Run decoupleR transcription factor analysis (Run_DecoupleR_AluOE.R)
