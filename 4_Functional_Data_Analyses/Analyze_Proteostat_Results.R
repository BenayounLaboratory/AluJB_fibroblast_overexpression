# ANALYSIS OF PROTEOSTAT DATA

# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(reshape) # Needed if data is melted
library(EnvStats) # for outlier detection

# Specify output directory
dir.output <- c('/Users/juanb/Desktop/Proteostat_Aggresome/Statistical_Analysis/')

# Load the data
exp_1 <- read.csv(file = '/Users/juanb/Desktop/Proteostat_Aggresome/Processed_with_Flowlogic/Experiment_1_MFI_Sample_Processing_Delayed.csv', header = TRUE, row.names = NULL)
exp_5 <- read.csv(file = '/Users/juanb/Desktop/Proteostat_Aggresome/Processed_with_Flowlogic/Experiment_5_MFI.csv', header = TRUE, row.names = NULL)
exp_6 <- read.csv(file = '/Users/juanb/Desktop/Proteostat_Aggresome/Processed_with_Flowlogic/Experiment_6_MFI.csv', header = TRUE, row.names = NULL)
exp_7 <- read.csv(file = '/Users/juanb/Desktop/Proteostat_Aggresome/Processed_with_Flowlogic/Experiment_7_MFI_Different_Cell_Scatter_Pattern.csv', header = TRUE, row.names = NULL)
exp_8 <- read.csv(file = '/Users/juanb/Desktop/Proteostat_Aggresome/Processed_with_Flowlogic/Experiment_8_MFI.csv', header = TRUE, row.names = NULL)
exp_9 <- read.csv(file = '/Users/juanb/Desktop/Proteostat_Aggresome/Processed_with_Flowlogic/Experiment_9_MFI.csv', header = TRUE, row.names = NULL)

    # Update colnames
    colnames(exp_1) <- c('Sample', 'Events', 'B2_MFI')
    colnames(exp_5) <- c('Sample', 'Events', 'B2_MFI')
    colnames(exp_6) <- c('Sample', 'Events', 'B2_MFI')
    colnames(exp_7) <- c('Sample', 'Events', 'B2_MFI')
    colnames(exp_8) <- c('Sample', 'Events', 'B2_MFI')
    colnames(exp_9) <- c('Sample', 'Events', 'B2_MFI')

    # Add a column for group assignment
    exp_1$Group <- c(rep('Empty', 6), rep('Alu', 6))
    exp_5$Group <- c(rep('Empty', 6), rep('Alu', 6))
    exp_6$Group <- c(rep('Empty', 6), rep('Alu', 6))
    exp_7$Group <- c(rep('Empty', 6), rep('Alu', 6))
    exp_8$Group <- c(rep('Empty', 6), rep('Alu', 6))
    exp_9$Group <- c(rep('Empty', 6), rep('Alu', 6))

# Generate a second combinated df with each experiment normalized to the media of the control group

    # Define df to hold median normalized data
    exp1.normalized <- exp_1
    exp5.normalized <- exp_5
    exp6.normalized <- exp_6
    exp7.normalized <- exp_7
    exp8.normalized <- exp_8
    exp9.normalized <- exp_9
    
    # Define the median for each experiment
    exp_1.med <- median(exp_1[which(exp_1$Group == 'Empty'), 'B2_MFI'])
    exp_5.med <- median(exp_5[which(exp_5$Group == 'Empty'), 'B2_MFI'])
    exp_6.med <- median(exp_6[which(exp_6$Group == 'Empty'), 'B2_MFI'])
    exp_7.med <- median(exp_7[which(exp_7$Group == 'Empty'), 'B2_MFI'])
    exp_8.med <- median(exp_8[which(exp_8$Group == 'Empty'), 'B2_MFI'])
    exp_9.med <- median(exp_9[which(exp_9$Group == 'Empty'), 'B2_MFI'])
    
    # Median-normalize each experiment
    exp1.normalized$B2_MFI <- exp1.normalized$B2_MFI/exp_1.med
    exp5.normalized$B2_MFI <- exp5.normalized$B2_MFI/exp_5.med
    exp6.normalized$B2_MFI <- exp6.normalized$B2_MFI/exp_6.med
    exp7.normalized$B2_MFI <- exp7.normalized$B2_MFI/exp_7.med
    exp8.normalized$B2_MFI <- exp8.normalized$B2_MFI/exp_8.med
    exp9.normalized$B2_MFI <- exp9.normalized$B2_MFI/exp_9.med
    
    # Combine normalized experiments
    exps_clean.normalized <-rbind(exp5.normalized, exp6.normalized, exp8.normalized, exp9.normalized)
    
    # Check for outliers
    my.outliers <- rosnerTest(exps_clean.normalized[which(exps_clean.normalized$Group == 'Empty'), 'B2_MFI'], k = 2, alpha = 0.05, warn = TRUE)
  
# Reorder the group levels in the dataframe to match the order I want in the plot
exp_1$Group <- factor(exp_1$Group , levels = c('Empty', 'Alu'))
exp_5$Group <- factor(exp_5$Group , levels = c('Empty', 'Alu'))
exp_6$Group <- factor(exp_6$Group , levels = c('Empty', 'Alu'))
exp_7$Group <- factor(exp_7$Group , levels = c('Empty', 'Alu'))
exp_8$Group <- factor(exp_8$Group , levels = c('Empty', 'Alu'))
exp_9$Group <- factor(exp_9$Group , levels = c('Empty', 'Alu'))

exp1.normalized$Group <- factor(exp1.normalized$Group , levels = c('Empty', 'Alu'))
exp5.normalized$Group <- factor(exp5.normalized$Group , levels = c('Empty', 'Alu'))
exp6.normalized$Group <- factor(exp6.normalized$Group , levels = c('Empty', 'Alu'))
exp7.normalized$Group <- factor(exp7.normalized$Group , levels = c('Empty', 'Alu'))
exp8.normalized$Group <- factor(exp8.normalized$Group , levels = c('Empty', 'Alu'))
exp9.normalized$Group <- factor(exp9.normalized$Group , levels = c('Empty', 'Alu'))


exps_clean.normalized$Group <- factor(exps_clean.normalized$Group , levels = c('Empty', 'Alu'))


# Run statistical analyses

    Stats.exp_1 <- wilcox.test(exp_1[which(exp_1$Group == 'Empty'), 'B2_MFI'],
                               exp_1[which(exp_1$Group == 'Alu'), 'B2_MFI'])
    
    Stats.exp_5 <- wilcox.test(exp_5[which(exp_5$Group == 'Empty'), 'B2_MFI'],
                               exp_5[which(exp_5$Group == 'Alu'), 'B2_MFI'])
    
    Stats.exp_6 <- wilcox.test(exp_6[which(exp_6$Group == 'Empty'), 'B2_MFI'],
                               exp_6[which(exp_6$Group == 'Alu'), 'B2_MFI'])
    
    Stats.exp_7 <- wilcox.test(exp_7[which(exp_7$Group == 'Empty'), 'B2_MFI'],
                               exp_7[which(exp_7$Group == 'Alu'), 'B2_MFI'])
    
    Stats.exp_8 <- wilcox.test(exp_8[which(exp_8$Group == 'Empty'), 'B2_MFI'],
                               exp_8[which(exp_8$Group == 'Alu'), 'B2_MFI'])
    
    Stats.exp_9 <- wilcox.test(exp_9[which(exp_9$Group == 'Empty'), 'B2_MFI'],
                               exp_9[which(exp_9$Group == 'Alu'), 'B2_MFI'])
    
    Stats.exps_clean.norm <- wilcox.test(exps_clean.normalized[which(exps_clean.normalized$Group == 'Empty'), 'B2_MFI'],
                                         exps_clean.normalized[which(exps_clean.normalized$Group == 'Alu'), 'B2_MFI'])


    Stats.exp_1$p.value
    Stats.exp_5$p.value
    Stats.exp_6$p.value
    Stats.exp_7$p.value
    Stats.exp_8$p.value
    Stats.exp_9$p.value
    Stats.exps_clean.norm$p.value





# GENERATE PLOTS

# Define general plots parameters
    
    # Specify xlabels tick positions
    plot.other.ticks <- c(1, 2)
    
    # Specify xlabels 
    plot.other.label <- c('Empty', 'AluJb')

    # Define point colors
    my.main.colors <- c('black', 'black')
    other.point.color <- rep(my.main.colors, 6)


  
# Start PDF to plot EXPERIMENT 1
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_1_Sample_Processing_Delayed", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp_1,
            ylim = c(0, 80),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 1', sep = ''),
            ylab = c('MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp_1, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 78, paste('p = ', signif(Stats.exp_1$p.value, 3), sep = ''), cex = 1)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 1 NORMALIZED
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_1_Sample_Processing_Delayed_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp1.normalized,
            ylim = c(0, 2),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 5', sep = ''),
            ylab = c('Relative MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp1.normalized, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 2, paste('p = ', signif(Stats.exp_1$p.value, 3), sep = ''), cex = 1)
    
    # Add line at y = 1
    abline(h = 1, lty = 2)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 5
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_5", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp_5,
            ylim = c(0, 80),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 5', sep = ''),
            ylab = c('MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp_5, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 78, paste('p = ', signif(Stats.exp_5$p.value, 3), sep = ''), cex = 1)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 5 NORMALIZED
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_5_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp5.normalized,
            ylim = c(0.5, 1.5),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 5', sep = ''),
            ylab = c('Relative MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp5.normalized, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 1.5, paste('p = ', signif(Stats.exp_5$p.value, 3), sep = ''), cex = 1)
    
    # Add line at y = 1
    abline(h = 1, lty = 2)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 6
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_6", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp_6,
            ylim = c(0, 80),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 6', sep = ''),
            ylab = c('MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp_6, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 78, paste('p = ', signif(Stats.exp_6$p.value, 3), sep = ''), cex = 1)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 6 NORMALIZED
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_6_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp6.normalized,
            ylim = c(0.5, 1.5),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 6', sep = ''),
            ylab = c('Relative MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp6.normalized, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 1.5, paste('p = ', signif(Stats.exp_6$p.value, 3), sep = ''), cex = 1)
    
    # Add line at y = 1
    abline(h = 1, lty = 2)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 7
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_7_Different_Cell_Scatter_Pattern", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp_7,
            ylim = c(0, 80),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 7', sep = ''),
            ylab = c('MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp_7, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 78, paste('p = ', signif(Stats.exp_7$p.value, 3), sep = ''), cex = 1)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 7
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_7_Different_Cell_Scatter_Pattern_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp7.normalized,
            ylim = c(0.5, 1.5),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 7', sep = ''),
            ylab = c('Relative MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp7.normalized, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 1.5, paste('p = ', signif(Stats.exp_7$p.value, 3), sep = ''), cex = 1)
    
    # Add line at y = 1
    abline(h = 1, lty = 2)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 8
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_8", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp_8,
            ylim = c(0, 80),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 8', sep = ''),
            ylab = c('MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp_8, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 78, paste('p = ', signif(Stats.exp_8$p.value, 3), sep = ''), cex = 1)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 8 NORMALIZED
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_8_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp8.normalized,
            ylim = c(0.5, 1.5),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 8', sep = ''),
            ylab = c('Relative MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp8.normalized, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 1.5, paste('p = ', signif(Stats.exp_8$p.value, 3), sep = ''), cex = 1)
    
    # Add line at y = 1
    abline(h = 1, lty = 2)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 9
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_9", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp_9,
            ylim = c(0, 80),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 9', sep = ''),
            ylab = c('MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp_9, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 78, paste('p = ', signif(Stats.exp_9$p.value, 3), sep = ''), cex = 1)
    
# End pdf
dev.off()





# Start PDF to plot EXPERIMENT 9 NORMALIZED
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_Exp_9_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exp9.normalized,
            ylim = c(0.5, 1.5),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 9', sep = ''),
            ylab = c('Relative MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exp9.normalized, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 1.5, paste('p = ', signif(Stats.exp_9$p.value, 3), sep = ''), cex = 1)
    
    # Add line at y = 1
    abline(h = 1, lty = 2)
    
# End pdf
dev.off()





# Start PDF to plot NORMALIZED COMBINED data
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Proteostat_5+6+8+9_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # make boxplots 
    boxplot(B2_MFI ~ Group, 
            col = c("dodgerblue", "dodgerblue4"), 
            data = exps_clean.normalized,
            ylim = c(0.5, 1.5),
            outline = F, 
            #col = boxplot.col,
            main = paste('Proteostat Exp 5+6+8+9', sep = ''),
            ylab = c('Relative MFI'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(B2_MFI ~ Group, data = exps_clean.normalized, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 1.5, paste('p = ', signif(Stats.exps_clean.norm$p.value, 3), sep = ''), cex = 1)
    
    # Add line at y = 1
    abline(h = 1, lty = 2)
    
# End pdf
dev.off()





# Clean the environment
rm(list=ls())           
