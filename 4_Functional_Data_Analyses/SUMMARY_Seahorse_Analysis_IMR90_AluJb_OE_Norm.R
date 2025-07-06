# Set strings as factors
options(stringsAsFactors = F)

# Load needed libraries
library(beeswarm)
library(bitops)

# Set a working directory
setwd('/Users/juanb/Library/CloudStorage/Dropbox/Research/2018-2024_Benayoun_Lab/Research/Data/Seahorse/2024_12_17_AluJb_OE_IMR90_WI38_Run_1_and_2_combined/Plot_IMR90_only/Seahorse_Results_AluJb_OE/')



# Load the parsed data
c1.data   <- read.table('2025-04-08_IMR90_AluJb_OE_Run_1_Parsed_Seahorse_Data.txt', sep = "\t", row.names = 1, header = TRUE)
c2.data   <- read.table('2025-04-08_IMR90_AluJb_OE_Run_2_Parsed_Seahorse_Data.txt', sep = "\t", row.names = 1, header = TRUE)

# Normalize each cohort to the IMR90 control group
c1.data$Norm_Basal_Respiration  <- c1.data$Basal_Respiration/median(c1.data$Basal_Respiration[c1.data$Group == "IMR90_U6"])
c1.data$Norm_ATP_Linked_Respiration  <- c1.data$ATP_Linked_Respiration/median(c1.data$ATP_Linked_Respiration[c1.data$Group == "IMR90_U6"])
c1.data$Norm_Proton_Leak  <- c1.data$Proton_Leak/median(c1.data$Proton_Leak[c1.data$Group == "IMR90_U6"])
c1.data$Norm_Maximal_Respiration  <- c1.data$Maximal_Respiration/median(c1.data$Maximal_Respiration[c1.data$Group == "IMR90_U6"])
c1.data$Norm_Spare_Respiratory_Capacity  <- c1.data$Spare_Respiratory_Capacity/median(c1.data$Spare_Respiratory_Capacity[c1.data$Group == "IMR90_U6"])
c1.data$Norm_Non_mitochondrial_OCR  <- c1.data$Non_mitochondrial_OCR/median(c1.data$Non_mitochondrial_OCR[c1.data$Group == "IMR90_U6"])
c1.data$Norm_Basal_ECAR  <- c1.data$Basal_ECAR/median(c1.data$Basal_ECAR[c1.data$Group == "IMR90_U6"])

c2.data$Norm_Basal_Respiration  <- c2.data$Basal_Respiration/median(c2.data$Basal_Respiration[c2.data$Group == "IMR90_U6"])
c2.data$Norm_ATP_Linked_Respiration  <- c2.data$ATP_Linked_Respiration/median(c2.data$ATP_Linked_Respiration[c2.data$Group == "IMR90_U6"])
c2.data$Norm_Proton_Leak  <- c2.data$Proton_Leak/median(c2.data$Proton_Leak[c2.data$Group == "IMR90_U6"])
c2.data$Norm_Maximal_Respiration  <- c2.data$Maximal_Respiration/median(c2.data$Maximal_Respiration[c2.data$Group == "IMR90_U6"])
c2.data$Norm_Spare_Respiratory_Capacity  <- c2.data$Spare_Respiratory_Capacity/median(c2.data$Spare_Respiratory_Capacity[c2.data$Group == "IMR90_U6"])
c2.data$Norm_Non_mitochondrial_OCR  <- c2.data$Non_mitochondrial_OCR/median(c2.data$Non_mitochondrial_OCR[c2.data$Group == "IMR90_U6"])
c2.data$Norm_Basal_ECAR  <- c2.data$Basal_ECAR/median(c2.data$Basal_ECAR[c2.data$Group == "IMR90_U6"])

# merge the data, if you have multiple cohorts
my.merged <- rbind(c1.data, c2.data)

# define factors for boxplots
my.merged$Group <- factor(my.merged$Group, levels = c("IMR90_U6", "IMR90_Alu", "WI38_U6", "WI38_Alu"))




# BASAL RESPIRATION

# get p-values
Alu.IMR90 <- wilcox.test(my.merged$Norm_Basal_Respiration[my.merged$Group == "IMR90_U6"],my.merged$Norm_Basal_Respiration[my.merged$Group == "IMR90_Alu"])

# Make a list of the data for this metric
basal_resp.plot <- list("IMR90_U6" = my.merged$Norm_Basal_Respiration[my.merged$Group == "IMR90_U6"],
                        "IMR90_Alu" = my.merged$Norm_Basal_Respiration[my.merged$Group == "IMR90_Alu"]
                        )

# Generate the plot
pdf(paste0(Sys.Date(),"_IMR90_Run_1+2_Parsed_Seahorse_Basal_Respiration_Norm-IMR90-U6.pdf"), height = 5, width = 2.5)

    boxplot(basal_resp.plot, 
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(-0.01, 2.5),
            ylab    = "Relative Basal Respiration (A.U.)",
            outline = F ,
            main    = "Basal Respiration 1&2")
    
    beeswarm(basal_resp.plot, add = T, pch = 16)
    
    # add pvalue
    text(1.5,  2.5, signif(Alu.IMR90$p.value,3))

    # Add bars for pvalue comparisons
    segments(x0 = 1, y0 = 2.4, x1 = 2, y1 = 2.4)

    # Add line at y = 1
    abline(h = 1, lty = 2)

dev.off()





# ATP Linked Respiration

# get p-values
Alu.IMR90 <- wilcox.test(my.merged$Norm_ATP_Linked_Respiration[my.merged$Group == "IMR90_U6"],my.merged$Norm_ATP_Linked_Respiration[my.merged$Group == "IMR90_Alu"])

# Make a list of the data for this metric
ATP_resp.plot <- list("IMR90_U6" = my.merged$Norm_ATP_Linked_Respiration[my.merged$Group == "IMR90_U6"],
                      "IMR90_Alu" = my.merged$Norm_ATP_Linked_Respiration[my.merged$Group == "IMR90_Alu"])

# Generate the plot
pdf(paste0(Sys.Date(),"_IMR90_Run_1+2_Parsed_Seahorse_ATP-Linked_Respiration_Norm-IMR90-U6.pdf"), height = 5, width = 2.5)

    boxplot(ATP_resp.plot, 
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(-0.01, 2.5),
            ylab    = "Relative ATP-Linked Respiration (A.U.)",
            outline = F ,
            main    = "ATP-Linked Respiration 1&2")
    
    beeswarm(ATP_resp.plot, add = T, pch = 16)
    
    # add pvalue
    text(1.5,  2.5, signif(Alu.IMR90$p.value,3))

    # Add bars for pvalue comparisons
    segments(x0 = 1, y0 = 2.4, x1 = 2, y1 = 2.4)

    # Add line at y = 1
    abline(h = 1, lty = 2)

dev.off()





# Proton Leak

# get p-values
Alu.IMR90 <- wilcox.test(my.merged$Norm_Proton_Leak[my.merged$Group == "IMR90_U6"],my.merged$Norm_Proton_Leak[my.merged$Group == "IMR90_Alu"])

# Make a list of the data for this metric
Proton_leak.plot <- list("IMR90_U6" = my.merged$Norm_Proton_Leak[my.merged$Group == "IMR90_U6"],
                      "IMR90_Alu" = my.merged$Norm_Proton_Leak[my.merged$Group == "IMR90_Alu"])

# Generate the plot
pdf(paste0(Sys.Date(),"_IMR90_Run_1+2_Parsed_Seahorse_Proton_Leak_Norm-IMR90-U6.pdf"), height = 5, width = 2.5)

    boxplot(Proton_leak.plot, 
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(-0.01, 2.5),
            ylab    = "Relative Proton Leak (A.U.)",
            outline = F ,
            main    = "Proton Leak 1&2")
    
    beeswarm(Proton_leak.plot, add = T, pch = 16)
    
    # add pvalue
    text(1.5,  2.5, signif(Alu.IMR90$p.value,3))

    # Add bars for pvalue comparisons
    segments(x0 = 1, y0 = 2.4, x1 = 2, y1 = 2.4)

    # Add line at y = 1
    abline(h = 1, lty = 2)

dev.off()





# Maximal Respiration

# get p-values
Alu.IMR90 <- wilcox.test(my.merged$Norm_Maximal_Respiration[my.merged$Group == "IMR90_U6"],my.merged$Norm_Maximal_Respiration[my.merged$Group == "IMR90_Alu"])

# Make a list of the data for this metric
Max_resp.plot <- list("IMR90_U6" = my.merged$Norm_Maximal_Respiration[my.merged$Group == "IMR90_U6"],
                      "IMR90_Alu" = my.merged$Norm_Maximal_Respiration[my.merged$Group == "IMR90_Alu"])

# Generate the plot
pdf(paste0(Sys.Date(),"_IMR90_Run_1+2_Parsed_Seahorse_Maximal_Respiration_Norm-IMR90-U6.pdf"), height = 5, width = 2.5)

    boxplot(Max_resp.plot, 
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(-0.01, 2.5),
            ylab    = "Relative Max Respiration (A.U.)",
            outline = F ,
            main    = "Max Respiration 1&2")
    
    beeswarm(Max_resp.plot, add = T, pch = 16)
    
    # add pvalue
    text(1.5,  2.5, signif(Alu.IMR90$p.value,3))

    # Add bars for pvalue comparisons
    segments(x0 = 1, y0 = 2.4, x1 = 2, y1 = 2.4)

    # Add line at y = 1
    abline(h = 1, lty = 2)

dev.off()





# Spare Respiratory Capacity

# get p-values
Alu.IMR90 <- wilcox.test(my.merged$Norm_Spare_Respiratory_Capacity[my.merged$Group == "IMR90_U6"],my.merged$Norm_Spare_Respiratory_Capacity[my.merged$Group == "IMR90_Alu"])

# Make a list of the data for this metric
Spare_capacity.plot <- list("IMR90_U6" = my.merged$Norm_Spare_Respiratory_Capacity[my.merged$Group == "IMR90_U6"],
                            "IMR90_Alu" = my.merged$Norm_Spare_Respiratory_Capacity[my.merged$Group == "IMR90_Alu"])

# Generate the plot
pdf(paste0(Sys.Date(),"_IMR90_Run_1+2_Parsed_Seahorse_Spare_Respiratory_Capacity_Norm-IMR90-U6.pdf"), height = 5, width = 2.5)

    boxplot(Spare_capacity.plot, 
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(-0.01, 2.5),
            ylab    = "Relative Spare Respiratory Capacity (A.U.)",
            outline = F ,
            main    = "Spare Respiratory Capacity 1&2")
    
    beeswarm(Spare_capacity.plot, add = T, pch = 16)
    
    # add pvalue
    text(1.5,  2.5, signif(Alu.IMR90$p.value,3))

    # Add bars for pvalue comparisons
    segments(x0 = 1, y0 = 2.4, x1 = 2, y1 = 2.4)

    # Add line at y = 1
    abline(h = 1, lty = 2)

dev.off()





# Non-mitochondrial OCR

# get p-values
Alu.IMR90 <- wilcox.test(my.merged$Norm_Non_mitochondrial_OCR[my.merged$Group == "IMR90_U6"],my.merged$Norm_Non_mitochondrial_OCR[my.merged$Group == "IMR90_Alu"])

# Make a list of the data for this metric
non_mito.plot <- list("IMR90_U6" = my.merged$Norm_Non_mitochondrial_OCR[my.merged$Group == "IMR90_U6"],
                      "IMR90_Alu" = my.merged$Norm_Non_mitochondrial_OCR[my.merged$Group == "IMR90_Alu"])

# Generate the plot
pdf(paste0(Sys.Date(),"_IMR90_Run_1+2_Parsed_Seahorse_Non-mito-OCR_Norm-IMR90-U6.pdf"), height = 5, width = 2.5)

    boxplot(non_mito.plot, 
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(-0.01, 2.5),
            ylab    = "Relative Non-mitochondrial OCR (A.U.)",
            outline = F ,
            main    = "Non-mitochondrial OCR 1&2")
    
    beeswarm(non_mito.plot, add = T, pch = 16)
    
    # add pvalue
    text(1.5,  2.5, signif(Alu.IMR90$p.value,3))

    # Add bars for pvalue comparisons
    segments(x0 = 1, y0 = 2.4, x1 = 2, y1 = 2.4)

    # Add line at y = 1
    abline(h = 1, lty = 2)

dev.off()





# BASAL ECAR

# get p-values
Alu.IMR90 <- wilcox.test(my.merged$Norm_Basal_ECAR[my.merged$Group == "IMR90_U6"],my.merged$Norm_Basal_ECAR[my.merged$Group == "IMR90_Alu"])

# Make a list of the data for this metric
basal_ecar.plot <- list("IMR90_U6" = my.merged$Norm_Basal_ECAR[my.merged$Group == "IMR90_U6"],
                      "IMR90_Alu" = my.merged$Norm_Basal_ECAR[my.merged$Group == "IMR90_Alu"])

# Generate the plot
pdf(paste0(Sys.Date(),"_IMR90_Run_1+2_Parsed_Seahorse_Basal-ECAR_Norm-IMR90-U6.pdf"), height = 5, width = 2.5)

    boxplot(basal_ecar.plot, 
            col     = c("dodgerblue", "dodgerblue4"), 
            ylim    = c(-0.01, 2.5),
            ylab    = "Relative Basal ECAR (A.U.)",
            outline = F ,
            main    = "Basal ECAR 1&2")
    
    beeswarm(basal_ecar.plot, add = T, pch = 16)
    
    # add pvalue
    text(1.5,  2.5, signif(Alu.IMR90$p.value,3))

    # Add bars for pvalue comparisons
    segments(x0 = 1, y0 = 2.4, x1 = 2, y1 = 2.4)

    # Add line at y = 1
    abline(h = 1, lty = 2)

dev.off()

