README -- Meta-analysis of publicly available aging RNASeq
############################

  1. Trim reads with fastp, check quality with FastQC, and map to hg38 with STAR (Trim_and_Map_GSE113957_Aging_Dermal_Fibroblasts_ShortSpots.sh, Trim_and_Map_GSE113957_Aging_Dermal_Fibroblasts_LongSpots.sh)

  2. Count reads mapping to genes or repeats using TETranscripts or TElocal (TETranscripts_GSE113957_Aging_Dermal_Fibroblasts.sh, TElocal_GSE113957_Aging_Dermal_Fibroblasts.sh)
      -  Prepare GTF file to analyze repeat subfamilies stratified by genomic context (TElocal_to_intergenic_intronic_v1.sh, TElocal_locations_to_GTF_v1.R)

  3. Run DESeq2 analysis for DE genes/TEs as a function of age
      - Check whether the biological sex listed for each sample is consistent with expression of male/female-specific genes (Run_Check_For_Sample_Mislabeling.R)
      - Filter lowly expressed genes, remove batch effects, run DESeq2, VST transform counts, and generate transposon heatmaps and scatterplots (Run_Low_Expression_Filter_and_DESeq2_TETranscripts_v4_linear.R, Run_Low_Expression_Filter_and_DESeq2_TElocal_v2_linear.R, Process_RNASeq_functions.R)
