# ANALYSIS OF CELL CYCLE DATA

# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(reshape) # Needed if data is melted

# Specify output directory
dir.output <- c('/Users/juanb/Desktop/Propidium_Iodide_Cell_Cycle/Statistical_Analysis/')

# Load the data
exp_2 <- read.csv(file = '/Users/juanb/Desktop/Propidium_Iodide_Cell_Cycle/Processed_with_Flowlogic/Exp2.csv', header = TRUE, row.names = NULL)
exp_5 <- read.csv(file = '/Users/juanb/Desktop/Propidium_Iodide_Cell_Cycle/Processed_with_Flowlogic/Exp5.csv', header = TRUE, row.names = NULL)
exp_6 <- read.csv(file = '/Users/juanb/Desktop/Propidium_Iodide_Cell_Cycle/Processed_with_Flowlogic/Exp6.csv', header = TRUE, row.names = NULL)
    
    # Update colnames
    colnames(exp_2) <- c('Sample', 'Events', 'G2_M', 'G0_G1', 'S')
    colnames(exp_5) <- c('Sample', 'Events', 'G2_M', 'G0_G1', 'S')
    colnames(exp_6) <- c('Sample', 'Events', 'G2_M', 'G0_G1', 'S')
    
    # Add a column for group assignment
    exp_2$Group <- c(rep('Empty', 6), rep('Alu', 6))
    exp_5$Group <- c(rep('Empty', 6), rep('Alu', 6))
    exp_6$Group <- c(rep('Empty', 6), rep('Alu', 6))
    
    # Make columns for Treatment x Cell Cycle Phase
    exp_2$G0G1_Group <- paste('G0/G1_', exp_2$Group, sep = '')
    exp_2$S_Group <- paste('S_', exp_2$Group, sep = '')
    exp_2$G2M_Group <- paste('G2/M_', exp_2$Group, sep = '')
    
    exp_5$G0G1_Group <- paste('G0/G1_', exp_5$Group, sep = '')
    exp_5$S_Group <- paste('S_', exp_5$Group, sep = '')
    exp_5$G2M_Group <- paste('G2/M_', exp_5$Group, sep = '')
    
    exp_6$G0G1_Group <- paste('G0/G1_', exp_6$Group, sep = '')
    exp_6$S_Group <- paste('S_', exp_6$Group, sep = '')
    exp_6$G2M_Group <- paste('G2/M_', exp_6$Group, sep = '')
    
    # Reorganize to two column format
    exp_2 <- data.frame(Condition = c(exp_2$G0G1_Group, exp_2$S_Group, exp_2$G2M_Group),
                        Percent = c(exp_2$G0_G1, exp_2$S, exp_2$G2_M)
                        )
    
    exp_5 <- data.frame(Condition = c(exp_5$G0G1_Group, exp_5$S_Group, exp_5$G2M_Group),
                        Percent = c(exp_5$G0_G1, exp_5$S, exp_5$G2_M)
                        )
   
    exp_6 <- data.frame(Condition = c(exp_6$G0G1_Group, exp_6$S_Group, exp_6$G2M_Group),
                        Percent = c(exp_6$G0_G1, exp_6$S, exp_6$G2_M)
                        )

# Define general plots parameters
    
    # Specify xlabels tick positions
    plot.other.ticks <- c(1, 2)
    
    # Specify xlabels 
    plot.other.label <- c('Empty', 'AluJb')

    # Define point colors
    my.main.colors <- c('black')
    other.point.color <- rep(my.main.colors, 12)

    # boxplot colors
    box.col <- c("dodgerblue", "dodgerblue4")
  
    
    
    
    
# PLOT EXPERIMENT 2
    
# Run statistical analyses
Stats.G0G1 <- wilcox.test(exp_2[which(exp_2$Condition == 'G0/G1_Empty'), 'Percent'],
                          exp_2[which(exp_2$Condition == 'G0/G1_Alu'), 'Percent'])

Stats.S <- wilcox.test(exp_2[which(exp_2$Condition == 'S_Empty'), 'Percent'],
                       exp_2[which(exp_2$Condition == 'S_Alu'), 'Percent'])

Stats.G2M <- wilcox.test(exp_2[which(exp_2$Condition == 'G2/M_Empty'), 'Percent'],
                         exp_2[which(exp_2$Condition == 'G2/M_Alu'), 'Percent'])



Stats.G0G1$p.value
Stats.S$p.value
Stats.G2M$p.value
    
    
# Start PDF to for G0/G1
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_2_G0G1", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_2[1:12, ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('G0/G1_Empty', 'G0/G1_Alu'))

    # make boxplots 
    boxplot(Percent ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(65, 80),
            outline = F, 
            #col = boxplot.col,
            main = paste('G0/G1 Phase Abundances', sep = ''),
            ylab = c('% of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Percent ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 80, paste('p = ', signif(Stats.G0G1$p.value, 3), sep = ''), cex = 1)

    # Add line at control median
    abline(h = median(plot.data$Percent[1:6]), lty = 2)
    
# End pdf
dev.off()


# Start PDF to for S Phase
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_2_S", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_2[13:24, ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('S_Empty', 'S_Alu'))
    
    # make boxplots 
    boxplot(Percent ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(5, 10),
            outline = F, 
            #col = boxplot.col,
            main = paste('S Phase Abundance', sep = ''),
            ylab = c('% of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Percent ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 10, paste('p = ', signif(Stats.S$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = median(plot.data$Percent[1:6]), lty = 2)
    
# End pdf
dev.off()


# Start PDF to for G2/M
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_2_G2M", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_2[25:36, ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('G2/M_Empty', 'G2/M_Alu'))
    
    # make boxplots 
    boxplot(Percent ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(15, 20),
            outline = F, 
            #col = boxplot.col,
            main = paste('G2/M Phase Abundance', sep = ''),
            ylab = c('% of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Percent ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 20, paste('p = ', signif(Stats.G2M$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = median(plot.data$Percent[1:6]), lty = 2)
    
# End pdf
dev.off()





# PLOT EXPERIMENT 5
    
# Run statistical analyses
Stats.G0G1 <- wilcox.test(exp_5[which(exp_5$Condition == 'G0/G1_Empty'), 'Percent'],
                          exp_5[which(exp_5$Condition == 'G0/G1_Alu'), 'Percent'])

Stats.S <- wilcox.test(exp_5[which(exp_5$Condition == 'S_Empty'), 'Percent'],
                       exp_5[which(exp_5$Condition == 'S_Alu'), 'Percent'])

Stats.G2M <- wilcox.test(exp_5[which(exp_5$Condition == 'G2/M_Empty'), 'Percent'],
                         exp_5[which(exp_5$Condition == 'G2/M_Alu'), 'Percent'])



Stats.G0G1$p.value
Stats.S$p.value
Stats.G2M$p.value
    
    
# Start PDF to for G0/G1
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_5_G0G1", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_5[1:12, ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('G0/G1_Empty', 'G0/G1_Alu'))
    
    # make boxplots 
    boxplot(Percent ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(55, 70),
            outline = F, 
            #col = boxplot.col,
            main = paste('G0/G1 Phase Abundances', sep = ''),
            ylab = c('% of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Percent ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 70, paste('p = ', signif(Stats.G0G1$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = median(plot.data$Percent[1:6]), lty = 2)
    
# End pdf
dev.off()


# Start PDF to for S Phase
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_5_S", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_5[13:24, ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('S_Empty', 'S_Alu'))
    
    # make boxplots 
    boxplot(Percent ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(10, 20),
            outline = F, 
            #col = boxplot.col,
            main = paste('S Phase Abundance', sep = ''),
            ylab = c('% of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Percent ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 20, paste('p = ', signif(Stats.S$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = median(plot.data$Percent[1:6]), lty = 2)
    
# End pdf
dev.off()


# Start PDF to for G2/M
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_5_G2M", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_5[25:36, ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('G2/M_Empty', 'G2/M_Alu'))
    
    # make boxplots 
    boxplot(Percent ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(15, 25),
            outline = F, 
            #col = boxplot.col,
            main = paste('G2/M Phase Abundance', sep = ''),
            ylab = c('% of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Percent ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 25, paste('p = ', signif(Stats.G2M$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = median(plot.data$Percent[1:6]), lty = 2)
    
# End pdf
dev.off()





# PLOT EXPERIMENT 6
    
# Run statistical analyses
Stats.G0G1 <- wilcox.test(exp_6[which(exp_6$Condition == 'G0/G1_Empty'), 'Percent'],
                          exp_6[which(exp_6$Condition == 'G0/G1_Alu'), 'Percent'])

Stats.S <- wilcox.test(exp_6[which(exp_6$Condition == 'S_Empty'), 'Percent'],
                       exp_6[which(exp_6$Condition == 'S_Alu'), 'Percent'])

Stats.G2M <- wilcox.test(exp_6[which(exp_6$Condition == 'G2/M_Empty'), 'Percent'],
                         exp_6[which(exp_6$Condition == 'G2/M_Alu'), 'Percent'])

Stats.G0G1$p.value
Stats.S$p.value
Stats.G2M$p.value
    
    
# Start PDF to for G0/G1
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_6_G0G1", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_6[1:12, ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('G0/G1_Empty', 'G0/G1_Alu'))
    
    # make boxplots 
    boxplot(Percent ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(50, 83),
            outline = F, 
            #col = boxplot.col,
            main = paste('G0/G1 Phase Abundances', sep = ''),
            ylab = c('% of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Percent ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 83, paste('p = ', signif(Stats.G0G1$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = median(plot.data$Percent[1:6]), lty = 2)
    
# End pdf
dev.off()


# Start PDF to for S Phase
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_6_S", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_6[13:24, ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('S_Empty', 'S_Alu'))
    
    # make boxplots 
    boxplot(Percent ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(10, 30),
            outline = F, 
            #col = boxplot.col,
            main = paste('S Phase Abundance', sep = ''),
            ylab = c('% of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Percent ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 30, paste('p = ', signif(Stats.S$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = median(plot.data$Percent[1:6]), lty = 2)
    
# End pdf
dev.off()


# Start PDF to for G2/M
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_6_G2M", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_6[25:36, ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('G2/M_Empty', 'G2/M_Alu'))
    
    # make boxplots 
    boxplot(Percent ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(7.5, 20),
            outline = F, 
            #col = boxplot.col,
            main = paste('G2/M Phase Abundance', sep = ''),
            ylab = c('% of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Percent ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 20, paste('p = ', signif(Stats.G2M$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = median(plot.data$Percent[1:6]), lty = 2)
    
# End pdf
dev.off()






# PLOT EXPERIMENT 2+5+6 NORMALIZED
    
# Write a function to do a within-experiment and within-phase normalization
normalize.exp <- function(exp_data){
  
  # Make a column to hold the normalized data
  exp_data$Normalized_value <- NA
  
  # Find the median within the control group for each phase
  median.G0G1 <- median(exp_data[which(exp_data$Condition == "G0/G1_Empty"), 'Percent'])
  median.S <- median(exp_data[which(exp_data$Condition == "S_Empty"), 'Percent'])
  median.G2M <- median(exp_data[which(exp_data$Condition == "G2/M_Empty"), 'Percent'])
  
  # Normalize each phase by dividing by the median
  exp_data[which(exp_data$Condition == "G0/G1_Empty" | exp_data$Condition == "G0/G1_Alu"), 'Normalized_value'] <- exp_data[which(exp_data$Condition == "G0/G1_Empty" | exp_data$Condition == "G0/G1_Alu"), 'Percent']/median.G0G1
  exp_data[which(exp_data$Condition == "S_Empty" | exp_data$Condition == "S_Alu"), 'Normalized_value'] <- exp_data[which(exp_data$Condition == "S_Empty" | exp_data$Condition == "S_Alu"), 'Percent']/median.S
  exp_data[which(exp_data$Condition == "G2/M_Empty" | exp_data$Condition == "G2/M_Alu"), 'Normalized_value'] <- exp_data[which(exp_data$Condition == "G2/M_Empty" | exp_data$Condition == "G2/M_Alu"), 'Percent']/median.G2M

  # Return the updated dataframe
  return(exp_data)
  
}

# Conduct normalization
exp_2 <- normalize.exp(exp_data = exp_2)
exp_5 <- normalize.exp(exp_data = exp_5)
exp_6 <- normalize.exp(exp_data = exp_6)

# Combine normalized experiments
exp_combined <- rbind(exp_2, exp_5, exp_6)

# Run statistical analyses
Stats.G0G1 <- wilcox.test(exp_combined[which(exp_combined$Condition == 'G0/G1_Empty'), 'Normalized_value'],
                          exp_combined[which(exp_combined$Condition == 'G0/G1_Alu'), 'Normalized_value'])

Stats.S <- wilcox.test(exp_combined[which(exp_combined$Condition == 'S_Empty'), 'Normalized_value'],
                       exp_combined[which(exp_combined$Condition == 'S_Alu'), 'Normalized_value'])

Stats.G2M <- wilcox.test(exp_combined[which(exp_combined$Condition == 'G2/M_Empty'), 'Normalized_value'],
                         exp_combined[which(exp_combined$Condition == 'G2/M_Alu'), 'Normalized_value'])

Stats.G0G1$p.value
Stats.S$p.value
Stats.G2M$p.value
    
    
# Start PDF to for G0/G1
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_2+5+6_G0G1-Phase_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_combined[which(exp_combined$Condition == "G0/G1_Empty" | exp_combined$Condition == "G0/G1_Alu"), ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('G0/G1_Empty', 'G0/G1_Alu'))
    
    # Find the median within the control group
    plot.median <- median(plot.data[which(plot.data$Condition == "G0/G1_Empty"), 'Normalized_value'])
    
    # make boxplots 
    boxplot(Normalized_value ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(0.7, 1.3),
            outline = F, 
            #col = boxplot.col,
            main = paste('G0/G1 Phase Abundances', sep = ''),
            ylab = c('Relative # of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Normalized_value ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 1.3, paste('p = ', signif(Stats.G0G1$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = plot.median, lty = 2)
    
# End pdf
dev.off()


# Start PDF to for S Phase
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_2+5+6_S-Phase_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_combined[which(exp_combined$Condition == "S_Empty" | exp_combined$Condition == "S_Alu"), ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('S_Empty', 'S_Alu'))
    
    # Find the median within the control group
    plot.median <- median(plot.data[which(plot.data$Condition == "S_Empty"), 'Normalized_value'])
    
    # make boxplots 
    boxplot(Normalized_value ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(0.5, 1.75),
            outline = F, 
            #col = boxplot.col,
            main = paste('S Phase Abundance', sep = ''),
            ylab = c('Relative # of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Normalized_value ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 1.75, paste('p = ', signif(Stats.S$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = plot.median, lty = 2)
    
# End pdf
dev.off()


# Start PDF to for G2/M
pdf(paste(dir.output, "Boxplot_IMR90_Alu_Cell_Cycle_Exp_2+5+6_G2M_Normalized", ".pdf", sep=""), width = 2.5, height = 5)

    # subset the data
    plot.data <- exp_combined[which(exp_combined$Condition == "G2/M_Empty" | exp_combined$Condition == "G2/M_Alu"), ]
    plot.data$Condition <- factor(plot.data$Condition , levels = c('G2/M_Empty', 'G2/M_Alu'))
    
    # Find the median within the control group
    plot.median <- median(plot.data[which(plot.data$Condition == "G2/M_Empty"), 'Normalized_value'])
    
    # make boxplots 
    boxplot(Normalized_value ~ Condition, 
            col = box.col, 
            data = plot.data,
            ylim = c(0.6, 1.3),
            outline = F, 
            #col = boxplot.col,
            main = paste('G2/M Phase Abundance', sep = ''),
            ylab = c('Relative # of cells'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(Normalized_value ~ Condition, data = plot.data, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 1.3, paste('p = ', signif(Stats.G2M$p.value, 3), sep = ''), cex = 1)
    
    # Add line at control median
    abline(h = plot.median, lty = 2)
    
# End pdf
dev.off()





# Clean the environment
rm(list=ls())           
