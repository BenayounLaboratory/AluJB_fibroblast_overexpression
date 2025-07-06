# Define a function for error bars
error.bar <- function(x, y, error, my.length=0.05, my.color){
  arrows(x, y+error, x, y-error, angle=90, code = 3, length = my.length, col = my.color)
}


# Define another function for error bars
error.bar.2d <- function(x, y, error, my.length=0.05, my.color){
  arrows(x+error, y, x-error, y, angle=90, code = 3, length = my.length, col = my.color)
}


# Define a function to parse Seahorse excel files
parse_seahorse <- function (seahorse.data.file, seahorse.plate.file, cohort.name) {
  
  ############## parse Seahorse data and extract functional metrics

  # Specify sheets by their name for the normalized rate data and the plate-sample map
  fibroblast.age.data  <- read_excel(seahorse.data.file , sheet = "Normalized Rate")
  fibroblast.age.plate <- read_excel(seahorse.plate.file, sheet = "Well_Sample",  col_names = F)
  colnames(fibroblast.age.plate) <- c("Well","Sample")
  
  # merge with annotation
  fibroblast.age.annot <- merge(fibroblast.age.data, fibroblast.age.plate, by = "Well")
  
  # remove background wells and empty wells
  fibroblast.age.annot <- fibroblast.age.annot[fibroblast.age.annot$Sample != "Background",]
  fibroblast.age.annot <- fibroblast.age.annot[fibroblast.age.annot$Sample != "Empty",]
  fibroblast.age.annot <- fibroblast.age.annot[fibroblast.age.annot$Sample != "Unassigned",]
  
  # remove negative data points
  my.neg.data  <- bitOr(fibroblast.age.annot$OCR < 0, fibroblast.age.annot$ECAR < 0) >0
  my.neg.wells <- unique(fibroblast.age.annot$Well[my.neg.data])
  
  fibroblast.age.annot <- fibroblast.age.annot[!(fibroblast.age.annot$Well %in% my.neg.wells),]

  # Parse data
  my.samples <- unique(fibroblast.age.annot$Sample)
  
  parsed.metab           <- data.frame(matrix(0, length(my.samples), 9))
  rownames(parsed.metab) <- my.samples
  colnames(parsed.metab) <- c("Group"                        ,
                              "Cohort"                       ,
                              "Basal_Respiration"            ,
                              "ATP_Linked_Respiration"       ,
                              "Proton_Leak"                  ,
                              "Maximal_Respiration"          ,
                              "Spare_Respiratory_Capacity"   ,
                              "Non_mitochondrial_OCR"        ,
                              "Basal_ECAR")
  
  parsed.metab$Cohort <- cohort.name
  parsed.metab$Group  <- "IMR90_U6"
  parsed.metab$Group[grep("IMR90_Alu",my.samples)] <- "IMR90_Alu"
  parsed.metab$Group[grep("WI38_U6",my.samples)] <- "WI38_U6"
  parsed.metab$Group[grep("WI38_Alu",my.samples)] <- "WI38_Alu"

  # use formulas from Stress Test product page
  for (i in 1:length(my.samples)) {
    
    my.sample <- my.samples[i]
    
    # extract all data of technical replicate wells AND time replicate measurements
    samp.data <- fibroblast.age.annot[fibroblast.age.annot$Sample == my.sample, ]
    
    # Calculate functional parameters, using the median value of the technical repliates for all relevent timepoints
    parsed.metab[i,]$Basal_Respiration              <- median(samp.data[samp.data$Measurement %in% c(1,2,3),]$OCR) - median(samp.data[samp.data$Measurement %in% c(10,11,12),]$OCR)
    parsed.metab[i,]$ATP_Linked_Respiration         <- median(samp.data[samp.data$Measurement %in% c(1,2,3),]$OCR)    - median(samp.data[samp.data$Measurement %in% c(4,5,6),]$OCR)
    parsed.metab[i,]$Proton_Leak                    <- median(samp.data[samp.data$Measurement %in% c(4,5,6),]$OCR)    - median(samp.data[samp.data$Measurement %in% c(10,11,12),]$OCR)
    parsed.metab[i,]$Maximal_Respiration            <- median(samp.data[samp.data$Measurement %in% c(7,8,9),]$OCR)    - median(samp.data[samp.data$Measurement %in% c(10,11,12),]$OCR)
    parsed.metab[i,]$Spare_Respiratory_Capacity     <- median(samp.data[samp.data$Measurement %in% c(7,8,9),]$OCR) - median(samp.data[samp.data$Measurement %in% c(1,2,3),]$OCR)
    parsed.metab[i,]$Non_mitochondrial_OCR          <- median(samp.data[samp.data$Measurement %in% c(10,11,12),]$OCR) 
    
    parsed.metab[i,]$Basal_ECAR                     <- median(samp.data[samp.data$Measurement %in% c(1,2,3),]$ECAR)
    
  }

  write.table(parsed.metab, file = paste0(Sys.Date(), '_', cohort.name, "_Parsed_Seahorse_Data.txt"), sep = "\t", quote = F, col.names = NA)
  
  
  
  
  
  ############## parse and plot ECAR and OCR curves per experiment
  
  # need to summarize each sample at each time point
  time.pts <- sort(unique(fibroblast.age.annot$Time))
  ECAR.samp           <- data.frame(matrix(0, length(time.pts), 1+ length(my.samples) ))
  OCR.samp            <- data.frame(matrix(0, length(time.pts), 1+ length(my.samples) ))
  colnames(ECAR.samp) <- c("Time", my.samples)
  colnames(OCR.samp)  <- c("Time", my.samples)
  ECAR.samp$Time      <- time.pts
  OCR.samp$Time       <- time.pts
  
  for (i in 1:length(my.samples)) { # loop over bio samples
    for (j in 1:length(time.pts)) { # loop over time points
      
      tmp.data <- fibroblast.age.annot[bitAnd(fibroblast.age.annot$Sample == my.samples[i], fibroblast.age.annot$Time == time.pts[j]) >0,]  # obtain data (including technical replicates) at time point j for sample i
      
      OCR.samp [j, 1+i] <- median(tmp.data[,"OCR"]) # take the median OCR at time point j for sample i
      ECAR.samp[j, 1+i] <- median(tmp.data[,"ECAR"]) # take the median ECAR at time point j for sample i
      
    }
  }
  
  # OCR: get mean and SEM per biological group for the cohort
  OCR.cohort  <- data.frame(matrix(0, length(time.pts), 9 ))
  colnames(OCR.cohort  ) <- c("Time", paste0(c("IMR90_U6", "IMR90_Alu", "WI38_U6", "WI38_Alu"),"_Av"), paste0(c("IMR90_U6", "IMR90_Alu", "WI38_U6", "WI38_Alu"),"_SEM"))
  OCR.cohort $Time      <- time.pts
  
  OCR.cohort$IMR90_U6_Av <- apply(OCR.samp[,grep("IMR90_U6", colnames(OCR.samp))],1,mean)
  OCR.cohort$IMR90_Alu_Av <- apply(OCR.samp[,grep("IMR90_Alu", colnames(OCR.samp))],1,mean)
  OCR.cohort$WI38_U6_Av <- apply(OCR.samp[,grep("WI38_U6", colnames(OCR.samp))],1,mean)
  OCR.cohort$WI38_Alu_Av <- apply(OCR.samp[,grep("WI38_Alu", colnames(OCR.samp))],1,mean)
  
  OCR.cohort$IMR90_U6_SEM <- apply(OCR.samp[,grep("IMR90_U6", colnames(OCR.samp))],1,sem)
  OCR.cohort$IMR90_Alu_SEM <- apply(OCR.samp[,grep("IMR90_Alu", colnames(OCR.samp))],1,sem)
  OCR.cohort$WI38_U6_SEM <- apply(OCR.samp[,grep("WI38_U6", colnames(OCR.samp))],1,sem)
  OCR.cohort$WI38_Alu_SEM <- apply(OCR.samp[,grep("WI38_Alu", colnames(OCR.samp))],1,sem)
  
  # ECAR: get mean and SEM per biological group for the cohort
  ECAR.cohort <- data.frame(matrix(0, length(time.pts), 9 ))
  colnames(ECAR.cohort ) <- c("Time", paste0(c("IMR90_U6", "IMR90_Alu", "WI38_U6", "WI38_Alu"),"_Av"), paste0(c("IMR90_U6", "IMR90_Alu", "WI38_U6", "WI38_Alu"),"_SEM"))
  ECAR.cohort$Time      <- time.pts
  
  ECAR.cohort$IMR90_U6_Av <- apply(ECAR.samp[,grep("IMR90_U6", colnames(ECAR.samp))],1,mean)
  ECAR.cohort$IMR90_Alu_Av <- apply(ECAR.samp[,grep("IMR90_Alu", colnames(ECAR.samp))],1,mean)
  ECAR.cohort$WI38_U6_Av <- apply(ECAR.samp[,grep("WI38_U6", colnames(ECAR.samp))],1,mean)
  ECAR.cohort$WI38_Alu_Av <- apply(ECAR.samp[,grep("WI38_Alu", colnames(ECAR.samp))],1,mean)
  
  ECAR.cohort$IMR90_U6_SEM <- apply(ECAR.samp[,grep("IMR90_U6", colnames(ECAR.samp))],1,sem)
  ECAR.cohort$IMR90_Alu_SEM <- apply(ECAR.samp[,grep("IMR90_Alu", colnames(ECAR.samp))],1,sem)
  ECAR.cohort$WI38_U6_SEM <- apply(ECAR.samp[,grep("WI38_U6", colnames(ECAR.samp))],1,sem)
  ECAR.cohort$WI38_Alu_SEM <- apply(ECAR.samp[,grep("WI38_Alu", colnames(ECAR.samp))],1,sem)
  
  
  # Define the y-axis limit for the OCR plot
  my.ylim.ocr <- round(2*max(OCR.cohort[,grep("Av", colnames(OCR.cohort))]))
  
  # Generate the OCR plot
  pdf(paste0(Sys.Date(), "_OCR_", cohort.name, "_per_Group.pdf"), width = 6, height = 5)
  
      # Plot the points and error bars for each group
      plot(OCR.cohort$Time,      
           OCR.cohort$IMR90_U6_Av, 
           type = 'b', 
           ylim = c(0, my.ylim.ocr), 
           pch = 16, 
           col = "dodgerblue",
           ylab = "Normalized OCR (pmol/min/DNA)", 
           xlab = "Time (minutes)", 
           main = paste0("OCR (", cohort.name, ")") )
      
      error.bar(OCR.cohort$Time, OCR.cohort$IMR90_U6_Av, OCR.cohort$IMR90_U6_SEM, my.color = "dodgerblue")
      
      points(OCR.cohort$Time,    OCR.cohort$IMR90_Alu_Av, type = 'b', pch = 16, col = "dodgerblue4")
      error.bar(OCR.cohort$Time, OCR.cohort$IMR90_Alu_Av, OCR.cohort$IMR90_Alu_SEM, my.color = "dodgerblue4")
      
      # Plot vertical lines corresponding to the injection times for each compount
      abline(v = 21.2  , lty = "dashed", col = "grey55", lwd = 0.5)
      abline(v = 45.67, lty = "dashed", col = "grey55", lwd = 0.5)
      abline(v = 70.16  , lty = "dashed", col = "grey55", lwd = 0.5)
      
      # Add labels for the injected compounds
      text(21.2  , 0.95*my.ylim.ocr, "Oligomycin", col = "grey55")
      text(45.67, 0.95*my.ylim.ocr, "FCCP"      , col = "grey55")
      text(70.16  , 0.95*my.ylim.ocr, "Rot/AA"    , col = "grey55")
      
  dev.off()
  
  # Define the y-axis limit for the ECAR plot
  my.ylim.ecar <- round(2*max(ECAR.cohort[,grep("Av", colnames(ECAR.cohort))]), digits = 1)
  
  # Generate the ECAR plot
  pdf(paste0(Sys.Date(),"_ECAR_", cohort.name, "_per_Group.pdf"), width = 6, height = 5)
  
      # Plot the points and error bars for each group
      plot(ECAR.cohort$Time,      
           ECAR.cohort$IMR90_U6_Av, type = 'b', 
           ylim = c(0, 0.3), 
           pch = 16, 
           col = "dodgerblue",
           ylab = "Normalized ECAR (mpH/min/DNA)", 
           xlab = "Time (minutes)", 
           main = paste0("ECAR (", cohort.name, ")") )
      
      error.bar(ECAR.cohort$Time, ECAR.cohort$IMR90_U6_Av, ECAR.cohort$IMR90_U6_SEM, my.color = "dodgerblue")
      
      points(ECAR.cohort$Time,    ECAR.cohort$IMR90_Alu_Av, type = 'b', pch = 16, col = "dodgerblue4")
      error.bar(ECAR.cohort$Time, ECAR.cohort$IMR90_Alu_Av, ECAR.cohort$IMR90_Alu_SEM, my.color = "dodgerblue4")
      
      # Plot vertical lines corresponding to the injection times for each compount
      abline(v = 21.2  , lty = "dashed", col = "grey55", lwd = 0.5)
      abline(v = 45.67, lty = "dashed", col = "grey55", lwd = 0.5)
      abline(v = 70.16  , lty = "dashed", col = "grey55", lwd = 0.5)
      
      # Add labels for the injected compounds
      text(21.2  , 0.95*my.ylim.ecar, "Oligomycin", col = "grey55")
      text(45.67, 0.95*my.ylim.ecar, "FCCP"      , col = "grey55")
      text(70.16  , 0.95*my.ylim.ecar, "Rot/AA"    , col = "grey55")
      
  dev.off()


  
}