setwd('/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/6_decoupleR/Results_serum_deprived/')
options(stringsAsFactors = F)

library("DESeq2")        #
library("decoupleR")     #
library("OmnipathR")     #

library(dplyr)
library(tibble)
library(ggplot2)
library(beeswarm)
library(ggrepel)

library("biomaRt")


library(bitops)

theme_set(theme_bw())   

# 2025-05-15
# Use decoupleR bioconductor package (and human TF target DB)
# to infer differential TF activity in Alu OE data
# https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/tf_bk.html


#####################################################################################################################
#### 0. Load up annotated DEseq2/VST objects

# read DEseq2 result table
my.deseq <- read.csv('/Users/juanb/Desktop/2024_NSC_Fibroblasts_Alu/Code/2_RNAseq_Analysis/3_DESeq/DESeq_Results/T72_Proteomics_DEGs_All.txt', sep = "\t") # 15183

head(my.deseq)
# X   baseMean log2FoldChange     lfcSE     stat       pvalue       padj external_gene_name hgnc_symbol entrezgene_id
# 1 ALR/Alpha:centr:Satellite   145.5397      0.5951499 0.2925336 2.034466 0.0419045857 0.27434504               <NA>        <NA>            NA
# 2              Alu:Alu:SINE   232.8252      0.7278283 0.2751594 2.645115 0.0081663148 0.13475525               <NA>        <NA>            NA
# 3            AluJb:Alu:SINE 20769.5319      0.5918941 0.1805678 3.277960 0.0010456028 0.05023857               <NA>        <NA>            NA
# 4            AluJo:Alu:SINE 12688.9835      0.6389903 0.1872740 3.412060 0.0006447387 0.03990498               <NA>        <NA>            NA
# 5            AluJr:Alu:SINE 14391.8985      0.5588897 0.1814658 3.079862 0.0020709632 0.07065940               <NA>        <NA>            NA
# 6           AluJr4:Alu:SINE  3039.0023      0.6609228 0.1903975 3.471280 0.0005179843 0.03591121               <NA>        <NA>            NA

deseq.for.dR <- my.deseq[!is.na(my.deseq$hgnc_symbol)      ,] # 14488
deseq.for.dR <- deseq.for.dR[deseq.for.dR$hgnc_symbol != "",] # 12966
# length(unique(deseq.for.dR$hgnc_symbol))
# [1] 12947

# collapse duplicates - when 2 entries, keep highest expressed
my.genes                  <- unique(deseq.for.dR$hgnc_symbol)
deseq.for.dR.cl           <- data.frame(matrix(0,length(my.genes),6))
colnames(deseq.for.dR.cl) <- colnames(deseq.for.dR[,c(2:7)])
rownames(deseq.for.dR.cl) <- my.genes

for (i in 1:length(my.genes)) {
  # get current gene
  my.gene <- my.genes[i]
  
  # get corresponding lines 
  my.subset <- deseq.for.dR[deseq.for.dR$hgnc_symbol %in% my.gene,]
  
  # keep highest expressed line when 2 entries have same HGNC symbol
  deseq.for.dR.cl[i,] <- my.subset[which.max(my.subset$baseMean),2:7]
  
}


#####################################################################################################################


#####################################################################################################################
#### 1. Process TF network data from human model
hs.net <- get_collectri(organism='human', split_complexes=FALSE)
hs.net
# # A tibble: 42,595 × 3
#    source target   mor
#    <chr>  <chr>  <dbl>
#  1 MYC    TERT       1
#  2 SPI1   BGLAP      1
#  3 SMAD3  JUN        1
#  4 SMAD4  JUN        1
#  5 STAT5A IL2        1
#  6 STAT5B IL2        1
#  7 RELA   FAS        1
#  8 WT1    NR0B1      1
#  9 NR0B2  CASP1      1
# 10 SP1    ALDOA      1
#####################################################################################################################

#####################################################################################################################
#### 2. Calculate TF network data

n_tfs <- 12

# infer pathway activities from the t-values of the DEGs with aging
# for each cell type

Top_TFs_res <- data.frame(matrix(0,0,6))
colnames(Top_TFs_res) <- c("statistic", "source","condition","score","p_value","Reg_Rank")

All_TFs_res <- data.frame(matrix(0,0,5))
colnames(All_TFs_res) <- c("statistic", "source","condition","score","p_value")

# Run fgsea scoring
contrast_acts <- run_fgsea(mat     = deseq.for.dR.cl[, 'stat', drop=FALSE],
                           net     = hs.net                                  , 
                           minsize = 5                                             )

All_TFs_res  <- contrast_acts[contrast_acts$statistic == 'norm_fgsea',]

# We select the norm_fgsea activities and then we show changes in activity with Alu OE
# Filter norm_fgsea
tf.scores         <- contrast_acts[contrast_acts$statistic == 'norm_fgsea',]
tf.sig            <- tf.scores[tf.scores$p_value < 0.05,]
tf.sig$Reg_Rank   <- NA

# Only keep TF regulons who TF is expressed
tf.sig <- tf.sig[tf.sig$source %in% row.names(deseq.for.dR.cl),]

# Add rank information
msk <- tf.sig$score > 0
tf.sig$Reg_Rank[msk ]  <- round(rank(-tf.sig[msk, 'score']))
tf.sig$Reg_Rank[!msk]  <- round(rank(-abs(tf.sig[!msk, 'score'])))

Top_TFs_res  <- tf.sig

# Filter top significant TFs in both signs
tfs <- tf.sig %>%
  arrange(Reg_Rank) %>%
  head(n_tfs) %>%
  pull(source)

f_contrast_acts <- tf.sig %>% filter(source %in% tfs)

# Plot top each direction
tf.reg.plot <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + geom_bar(aes(fill = score), stat = "identity") 
tf.reg.plot <- tf.reg.plot + scale_fill_gradient2(low = "#333399", mid = "whitesmoke", high = "#CC3333", midpoint = 0, limits = c(-2.5,2.5)) 
tf.reg.plot <- tf.reg.plot  + theme(axis.title = element_text(face = "bold", size = 12),
                                    axis.text.x = element_text(angle = 45, hjust = 1, size =10),
                                    axis.text.y = element_text(size =10, face= "bold"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank() ) 
tf.reg.plot <- tf.reg.plot + xlab("Enriched TF regulons") + ylab("decoupleR score") + theme(plot.title = element_text(hjust = 0.5))
tf.reg.plot <- tf.reg.plot + ylim(c(-3,3)) + ggtitle("IMR90 AluJB OE / DecoupleR") + coord_flip()

pdf(paste0(Sys.Date(),"_Barplot_decoupleR_fgsea_Top_12_TF_regulons_AluJB_OE_FDR5.pdf"), width = 4, height = 5)
print(tf.reg.plot)
dev.off()

write.table(Top_TFs_res, file = paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_AluJB_OE_FDR5.txt"), sep = "\t", quote = F, row.names = F)
write.table(All_TFs_res, file = paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_AluJB_OE_ALL_results.txt"), sep = "\t", quote = F, row.names = F)

########################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_R_session_Info_decoupleR_ALuJBOE.txt", sep =""))
sessionInfo()
sink()

