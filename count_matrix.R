#!/usr/bin/env Rscript

# script that produces count matrix from directory of bamfiles

library("Rsubread")

# function to find genes with successful alignments
nonzero_counts <- function(counts, cat) {
  
  # dataframe to store count data 
  df <- data.frame(Gene_ID = character(),
                   exons = integer(), 
                   genes = integer(), 
                   introns = integer())
  
  for (i in 1:length(counts)) {
    if (counts[i] != 0) {
      if (cat == "exons") {
        df <- rbind(df, data.frame("Gene_ID" = rownames(counts)[i], "exons" = counts[i], 
                                   "genes" = 0, "introns" = 0))
      } else if (cat == "genes") {
        df <- rbind(df, data.frame("Gene_ID" = rownames(counts)[i], "exons" = 0, 
                                   "genes" = counts[i], "introns" = 0))
      } else if (cat == "introns") {
        df <- rbind(df, data.frame("Gene_ID" = rownames(counts)[i], "exons" = 0, 
                                   "genes" = 0, "introns" = counts[i]))
      }
    }
  }
  
  return(df)
}

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
  stop("One argument must be supplied (directory of bamfiles)", call.=FALSE)
} else if (length(args) > 1) {
  stop("No more than one argument must be supplied (directory of bamfiles)", call.=FALSE)
} else if (!dir.exists(args[1])) {
  stop("Invalid directory", call.False)
}

# list of bam files in given directory
folder_files <- list.files(path = args[1], pattern = "\\.bam$")

print(folder_files)
