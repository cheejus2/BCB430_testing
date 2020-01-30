#!/usr/bin/env Rscript

# script that produces count matrix from directory of bamfiles
# arguments should be: 
# 1. directory of bamfiles 2. gtf file

library("Rsubread")
library("dplyr")

# function to create dataframe with exon, gene and intron data
combine_counts <- function(e_count, g_count) {
  
  # merge the exon and gene data
  temp <- merge(e_count, g_count, by = "row.names")
  
  # dataframe to store count data 
  df <- data.frame(Gene_ID = character(),
                   exons = integer(), 
                   genes = integer(), 
                   introns = integer())
  
  for (i in 1:length(temp[,1])) {
    if (temp[i, 2] != 0 && temp[i, 3] != 0) {
      if (temp[i, 3] - temp[i, 2] >= 0) {
        df <- rbind(df, data.frame("Gene_ID" = temp[i, 1], "exons" = temp[i, 2], 
                                   "genes" = temp[i, 3], "introns" = temp[i, 3] - temp[i, 2]))
      } else {
        df <- rbind(df, data.frame("Gene_ID" = temp[i, 1], "exons" = temp[i, 2], 
                                   "genes" = temp[i, 3], "introns" = 0))
      }
    }
  }
  
  return(df)
}

# reading in and checking arguments from the command line
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("two argument must be supplied (directory of bamfiles)", call.=FALSE)
} else if (length(args) > 2) {
  stop("2 arguments must be supplied: 1. directory of bamfiles 2. gtf file", call.=FALSE)
} else if (!dir.exists(args[1])) {
  stop("Invalid directory", call.False)
}

# list of bam files in given directory
folder_files <- list.files(path = args[1], pattern = "\\.bam$")

for (i in folder_files) {
  
  # variable holding exon counts
  temp_e <- featureCounts(files = paste(args[1], "/", i, sep = ""), 
                          annot.ext = args[2],
                          isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType="gene_id", 
                          isPairedEnd = TRUE)
  # variable holding gene counts
  temp_g <- featureCounts(files = paste(args[1], "/", i, sep = ""), 
                          annot.ext = args[2],
                          isGTFAnnotationFile = TRUE, GTF.featureType="gene", GTF.attrType="gene_id", 
                          isPairedEnd = TRUE)
  
  # temp dataframe holding exon, gene and intron counts
  temp_df <- combine_counts(temp_e$counts, temp_g$counts)
  
  # writing temp dataframe (count matrix) as csv file
  write.csv(temp_df, file = paste(args[1], "/", tools::file_path_sans_ext(i), "_count.csv", sep = ""))

}
