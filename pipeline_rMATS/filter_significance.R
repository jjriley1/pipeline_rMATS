#Load libraries
library(data.table)
library(dplyr)

#Take file input on command line to allow integration into pipelines (Rscript summarize_CTs.R 'file.csv' 'output_file.csv')
file_input = commandArgs(trailingOnly = T)

#Read as a data.table
file = fread(file_input[1])

file = file %>% subset(., select=which(!duplicated(names(.))))
file = file %>% dplyr::filter(PValue<0.05)
file = file %>% dplyr::filter(FDR<0.05)

#Write out
fwrite(file, file_input[2])


