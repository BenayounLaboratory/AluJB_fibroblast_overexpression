# Set strings as factors
options(stringsAsFactors = F)

# Load needed libraries
library(bitops)
library(goeveg)
library("readxl")

# Set a working directory
setwd('/Users/juanb/Library/CloudStorage/Dropbox/Research/2018-2024_Benayoun_Lab/Research/Data/Seahorse/2024_12_17_AluJb_OE_IMR90_WI38_Run_1_and_2_combined/Plot_IMR90_only/Seahorse_Results_AluJb_OE/')

# Load functions
source('/Users/juanb/Library/CloudStorage/Dropbox/Research/2018-2024_Benayoun_Lab/Research/Data/Seahorse/2024_12_17_AluJb_OE_IMR90_WI38_Run_1_and_2_combined/Plot_IMR90_only/Parse_Seahorse_Excel_IMR90_WI38_AluJb-OE_functions_vJB01.R')

# 2024-12-19
# Parse IMR90/WI38 seahorse data to get measurements
# 2 cohorts prepared by Juan and Eyael


####################################################################################
# AluJb OE IMR90 and WI38 Run #1
parse_seahorse(seahorse.data.file = "/Users/juanb/Library/CloudStorage/Dropbox/Research/2018-2024_Benayoun_Lab/Research/Data/Seahorse/2024_12_17_AluJb_OE_IMR90_WI38_Run_1_and_2_combined/Raw_Data/2024_12_05_AluJb_IMR90_WI38_XF_Cell_Mito_Stress_Test_1.xlsx", 
               seahorse.plate.file = "/Users/juanb/Library/CloudStorage/Dropbox/Research/2018-2024_Benayoun_Lab/Research/Data/Seahorse/2024_12_17_AluJb_OE_IMR90_WI38_Run_1_and_2_combined/Raw_Data/2024_12_05_AluJb_IMR90_WI38_XF_Cell_Mito_Stress_Test_1_Plate_Format.xlsx", 
               cohort.name = "IMR90_AluJb_OE_Run_1")

####################################################################################
# AluJb OE IMR90 and WI38 Run #2
parse_seahorse(seahorse.data.file = "/Users/juanb/Library/CloudStorage/Dropbox/Research/2018-2024_Benayoun_Lab/Research/Data/Seahorse/2024_12_17_AluJb_OE_IMR90_WI38_Run_1_and_2_combined/Raw_Data/2024_12_17_AluJb_IMR90_WI38_XF_Cell_Mito_Stress_Test_2.xlsx", 
               seahorse.plate.file = "/Users/juanb/Library/CloudStorage/Dropbox/Research/2018-2024_Benayoun_Lab/Research/Data/Seahorse/2024_12_17_AluJb_OE_IMR90_WI38_Run_1_and_2_combined/Raw_Data/2024_12_17_AluJb_IMR90_WI38_XF_Cell_Mito_Stress_Test_2_Plate_Format.xlsx", 
               cohort.name = "IMR90_AluJb_OE_Run_2")

####################################################################################
sink(paste0(Sys.Date(),"_IMR90_Seahorse_1+2_Parsing_sessionInfo.txt"))
sessionInfo()
sink()







  