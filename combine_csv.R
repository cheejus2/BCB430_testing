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
csv_files <- grep(csv_files, pattern='^combined', inv=T, value=T)

# initializing a dataframe to hold the combined csv files (count matrices)
combined_exons <- data.frame(Gene_ID = character())

combined_introns <- data.frame(Gene_ID = character())

combined_genes <- data.frame(Gene_ID = character())

for (i in csv_files) {
  
  # temp variable to hold count matrix data for sample
  count_matrix <- read.csv(file = paste(args[1], "/", i, sep = ""))
  sample_name <- tools::file_path_sans_ext(i)
  
  # temp df to hold exon data for sample
  exon_holder <- data.frame(Gene_ID = count_matrix$Gene_ID,
                            temp_name = count_matrix$exons)
  names(exon_holder)[names(exon_holder) == "temp_name"] <- sample_name
  
  # temp df to hold intron data for sample
  intron_holder <- data.frame(Gene_ID = count_matrix$Gene_ID,
                              temp_name = count_matrix$introns)
  names(intron_holder)[names(intron_holder) == "temp_name"] <- sample_name
  
  # temp df to hold gene data for sample
  gene_holder <- data.frame(Gene_ID = count_matrix$Gene_ID,
                              temp_name = count_matrix$genes)
  names(gene_holder)[names(gene_holder) == "temp_name"] <- sample_name
  
  # adding individual sample data to combined count matrices
  combined_exons <- merge(combined_exons, exon_holder, all = TRUE)
  combined_introns <- merge(combined_introns, intron_holder, all = TRUE)
  combined_genes <- merge(combined_genes, gene_holder, all = TRUE)
}

write.csv(combined_exons, file = paste(args[1], "/", "combined_exons.csv", sep = ""))
write.csv(combined_introns, file = paste(args[1], "/", "combined_introns.csv", sep = ""))
write.csv(combined_genes, file = paste(args[1], "/", "combined_genes.csv", sep = ""))
