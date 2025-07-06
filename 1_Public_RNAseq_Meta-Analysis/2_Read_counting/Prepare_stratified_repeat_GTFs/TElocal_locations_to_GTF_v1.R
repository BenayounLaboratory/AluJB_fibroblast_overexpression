# Load libraries
library(splitstackshape)

# Define the path of the repeat location file
repeat.locations.path <- '/Users/juanb/Library/CloudStorage/Dropbox/Research/Shared_Bioinformatic_Resources/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE.gtf.locInd.locations'

# Define the output path
output.dir <- '/Users/juanb/Library/CloudStorage/Dropbox/Research/Shared_Bioinformatic_Resources/STAR_genome_indices/TElocal_rmsk_files/'

# Load the TElocal locations file
repeat.locations <- read.csv(repeat.locations.path, header = FALSE, sep = '\t')

# Split #1 for ':'
repeat.locations <- as.data.frame(splitstackshape::cSplit(repeat.locations, splitCols = "V2", sep = ":", direction = "wide", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = TRUE))

# Split #2 for '-'
repeat.locations <- as.data.frame(splitstackshape::cSplit(repeat.locations, splitCols = "V2_2", sep = "-", direction = "wide", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = TRUE))

# Update the colnames
colnames(repeat.locations) <- repeat.locations[1, ]

# Remove the row with column labels
repeat.locations <- repeat.locations[-c(1), ]

# Add GTF features
  
    # Add source
    repeat.locations$source <- 'Repeat_masker'
    
    # Add feature
    repeat.locations$feature <- 'repeat'
    
    # Add score
    repeat.locations$score <- '.'
    
    # Add frame
    repeat.locations$frame <- '.'

# Organize columns in GTF format
repeat.locations <- repeat.locations[, c('chromsome', 'source', 'feature', 'start', 'stop', 'score', 'strand', 'frame', 'TE')]

# Save
write.table(repeat.locations, file = paste(output.dir, 'GRCh38_GENCODE_rmsk_TE.GTF', sep = ''), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t', )