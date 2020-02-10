#!/usr/bin/env Rscript

# script that produces count matrix from directory of bamfiles
# arguments should be: 
# 1. directory of bamfiles 2. gtf file

# reading in and checking arguments from the command line
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("One argument must be supplied (directory of bamfiles)", call.=FALSE)
} else if (length(args) > 1) {
  stop("One argument must be supplied (directory of bamfiles)", call.=FALSE)
} else if (!dir.exists(args[1])) {
  stop("Invalid directory", call.False)
}

# list of bam files in given directory
folder_files <- list.files(path = args[1], pattern = "\\.bam$")

for (i in folder_files) {
  
  fn = paste0(args[1], "/", i)
  
  
  # add code here
  
  
  # writing temp dataframe (count matrix) as csv file
  write.csv(temp_df, file = paste(args[1], "/", tools::file_path_sans_ext(i), "_count.csv", sep = ""))
  
}
