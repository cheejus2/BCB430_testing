#!/usr/bin/env Rscript

# script that combines a directory of csv files produced by "count_matrix.R" or "count_matrix_dir.R"
# into larger combined count matrix.

# Argument should be the dir containing the csv files (individual count matrices)

# to run script in command line type:
# "Rscript combine_csv.R arg1"
# where arg1 is the directory containing csv files

# reading in and checking arguments from the command line
args = commandArgs(trailingOnly = TRUE)

if (length(args) > 1 | length(args) < 1) {
  stop("1 argument must be supplied: 1. path to directory containing csv files", call.=FALSE)
} else if (!dir.exists(args[1])) {
  stop("Invalid csv directory", call.False)
}

# list of csv files in given directory
csv_files <- list.files(path = args[1], pattern = "\\.csv$")

# ignore csv produced by script combine_csv.R
csv_files <- grep(csv_files, pattern='combined_counts.csv', inv=T, value=T)

# initializing a dataframe to hold the combined csv files (count matrices)
combined <- data.frame(sample = character(),
                       Gene_ID = character(),
                       exons = integer(), 
                       introns = integer(),
                       exon_per = integer(),
                       intron_per = integer())

for (i in csv_files) {
  
  count_matrix <- read.csv(file = paste(args[1], "/", i, sep = ""))
  
  holder <- data.frame(sample = tools::file_path_sans_ext(i),
                       Gene_ID = count_matrix$Gene_ID,
                       exons = count_matrix$exons, 
                       introns = count_matrix$introns,
                       exon_per = (count_matrix$exons) / 
                         (count_matrix$exons + count_matrix$introns),
                       intron_per = (count_matrix$introns) /
                         (count_matrix$exons + count_matrix$introns))
  
  combined <- rbind(combined, holder)
}

write.csv(combined, file = paste(args[1], "/", "combined_counts.csv", sep = ""))
