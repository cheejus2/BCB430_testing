# resource: http://bioconductor.org/help/course-materials/2015/BioC2015/bioc2015rnaseq.html
# resource: https://bioinformatics-core-shared-training.github.io/RNAseq-R/rna-seq-annotation-visualisation.nb.html
# resource: https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf
# resource: https://www.bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html#7_composite_plots_for_multiple_chromosomes

library("GenomicFeatures")
library("Gviz")
library("GenomicRanges")
library("GenomicAlignments")
library("rtracklayer")
library("Rsubread")
library("erer")
library("plyr")
library("ggplot2")

# making txdb
txdb <- makeTxDbFromGFF("datafiles/Mus_musculus.GRCm38.98.gtf", format="gtf")

# setting genomic range, retreved from: https://www.ncbi.nlm.nih.gov/gene/14810
total.range <- GRanges("2", IRanges(25291177, 25319253))

# genome alignments
ga_wt_f <- readGAlignmentPairs("datafiles/wildtype_female_Grin1.bam", param=ScanBamParam(which=total.range))
ga_wt_f
ga_wt_m <- readGAlignmentPairs("datafiles/wildtype_male_Grin1.bam", param=ScanBamParam(which=total.range))
ga_wt_m
ga_mt_f <- readGAlignmentPairs("datafiles/mutant_female_Grin1.bam", param=ScanBamParam(which=total.range))
ga_mt_f

# computing coverage
cov_wt_f <- coverage(ga_wt_f)
cov_wt_m <- coverage(ga_wt_m)
cov_mt_f <- coverage(ga_mt_f)

zoom <- GRanges("2", IRanges(25291177, 25319253))

cov_wt_f[zoom]
cov_wt_m[zoom]
cov_mt_f[zoom]

cov_wt_f.numeric <- as.numeric(cov_wt_f[zoom][[1]])
cov_wt_m.numeric <- as.numeric(cov_wt_m[zoom][[1]])
cov_mt_f.numeric <- as.numeric(cov_mt_f[zoom][[1]])

median(cov_mt_f.numeric, cov_wt_f.numeric, cov_wt_m.numeric)

# plotting
plot(cov_wt_f.numeric, type="h")
plot(cov_wt_m.numeric, type="h")
plot(cov_mt_f.numeric, type="h")

# visualize read pileups
options(ucscChromosomeNames=FALSE)
grtrack <- GeneRegionTrack(txdb)
gtrack <- GenomeAxisTrack()
altrack_wt_f <- AlignmentsTrack("datafiles/wildtype_female_Grin1.bam")
altrack_wt_m <- AlignmentsTrack("datafiles/wildtype_male_Grin1.bam")
altrack_mt_f <- AlignmentsTrack("datafiles/mutant_female_Grin1.bam")
plotTracks(list(gtrack, grtrack, altrack_wt_f), from=25291177, to=25319253, chromosome="2")
plotTracks(list(gtrack, grtrack, altrack_wt_m), from=25291177, to=25319253, chromosome="2")
plotTracks(list(gtrack, grtrack, altrack_mt_f), from=25291177, to=25319253, chromosome="2")

plotTracks(list(gtrack, grtrack, altrack_wt_f, altrack_mt_f, altrack_wt_m), 
           from=25291177, to=25319253, chromosome="2")

ot <- OverlayTrack(trackList = list(altrack_wt_f, altrack_mt_f, altrack_wt_m))

plotTracks(list(gtrack, grtrack, ot), from=25291177, to=25319253, chromosome="2")

# just coverage
altrack_wt_f_c <- AlignmentsTrack("datafiles/wildtype_female_Grin1.bam", type=c("coverage"), legend = TRUE)
altrack_wt_m_c <- AlignmentsTrack("datafiles/wildtype_male_Grin1.bam", type=c("coverage"), legend = TRUE)
altrack_mt_f_c <- AlignmentsTrack("datafiles/mutant_female_Grin1.bam", type=c("coverage"), legend = TRUE)

displayPars(altrack_wt_f_c) <- list(col = "red", fill = "transparent", legend = TRUE)
displayPars(altrack_mt_f_c) <- list(col = "blue", fill = "transparent", legend = TRUE)
displayPars(altrack_wt_m_c) <- list(col = "orange", fill = "transparent", legend = TRUE)

plotTracks(list(gtrack, grtrack, altrack_wt_f_c, altrack_mt_f_c, altrack_wt_m_c), 
           from=25291177, to=25319253, chromosome="2", legend = TRUE)

ot_c <- OverlayTrack(trackList = list(altrack_wt_f_c, altrack_mt_f_c, altrack_wt_m_c))

plotTracks(list(gtrack, grtrack, ot_c), from=25291177, to=25319253, chromosome="2", legend = TRUE)

# visualize sashimi plot
altrack_wt_f_s <- AlignmentsTrack("datafiles/wildtype_female_Grin1.bam", type=c("coverage","sashimi"))
altrack_wt_m_s <- AlignmentsTrack("datafiles/wildtype_male_Grin1.bam", type=c("coverage","sashimi"))
altrack_mt_f_s <- AlignmentsTrack("datafiles/mutant_female_Grin1.bam", type=c("coverage","sashimi"))
plotTracks(list(gtrack, grtrack, altrack_wt_f_s), from=25291177, to=25319253, chromosome="2")
plotTracks(list(gtrack, grtrack, altrack_wt_m_s), from=25291177, to=25319253, chromosome="2")
plotTracks(list(gtrack, grtrack, altrack_mt_f_s), from=25291177, to=25319253, chromosome="2")

displayPars(altrack_wt_f_s) <- list(col = "red", fill = "transparent", legend = TRUE)
displayPars(altrack_mt_f_s) <- list(col = "blue", fill = "transparent", legend = TRUE)
displayPars(altrack_wt_m_s) <- list(col = "orange", fill = "transparent", legend = TRUE)

ot_s <- OverlayTrack(trackList = list(altrack_wt_f_s, altrack_mt_f_s, altrack_wt_m_s))

plotTracks(list(gtrack, grtrack, ot_s), from=25291177, to=25319253, chromosome="2")

# determining features contained within the gtf file
gtf <- rtracklayer::import('datafiles/Mus_musculus.GRCm38.98.gtf')
gtf_df = as.data.frame(gtf)
feature_list <- unique(gtf_df[7])
feature_list

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

# quantifying data with Rsubread
temp1 <- featureCounts(files="datafiles/wildtype_female_Grin1.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE)
temp2 <- featureCounts(files="datafiles/wildtype_male_Grin1.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE)
temp3 <- featureCounts(files="datafiles/mutant_female_Grin1.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE)
temp4 <- featureCounts(files="datafiles/wildtype_female_Grin1.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="gene", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE, allowMultiOverlap = FALSE)
temp5 <- featureCounts(files="datafiles/wildtype_male_Grin1.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="gene", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE, allowMultiOverlap = TRUE)
temp6 <- featureCounts(files="datafiles/mutant_female_Grin1.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="gene", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE, allowMultiOverlap = TRUE)
temp7 <- featureCounts(files="datafiles/wildtype_female_Grin1.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="CDS", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE)
temp8 <- featureCounts(files="datafiles/wildtype_male_Grin1.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="CDS", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE)
temp9 <- featureCounts(files="datafiles/mutant_female_Grin1.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="CDS", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE)

# comparison count table (mutant female chr2)
temp_mt_f_e <- featureCounts(files="datafiles/mutant_female_2.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE)
temp_mt_f_g <- featureCounts(files="datafiles/mutant_female_2.bam",
                       annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                       isGTFAnnotationFile = TRUE, GTF.featureType="gene", GTF.attrType="gene_id", 
                       isPairedEnd = TRUE)

result_e <- nonzero_counts(temp_mt_f_e$counts, "exons")
result_g <- nonzero_counts(temp_mt_f_g$counts, "genes")
result <- merge(result_e, result_g, by = "Gene_ID", all = TRUE)

# comparison count table (mutant female full genome)
mt_f_e <- featureCounts(files="datafiles/mutant_female.bam",
                             annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                             isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType="gene_id", 
                             isPairedEnd = TRUE)
mt_f_g <- featureCounts(files="datafiles/mutant_female.bam",
                             annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                             isGTFAnnotationFile = TRUE, GTF.featureType="gene", GTF.attrType="gene_id", 
                             isPairedEnd = TRUE)

# comparison count table (wildtype female full genome)
wt_f_e <- featureCounts(files="datafiles/wildtype_female.bam",
                        annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                        isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType="gene_id", 
                        isPairedEnd = TRUE)
wt_f_g <- featureCounts(files="datafiles/wildtype_female.bam",
                        annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                        isGTFAnnotationFile = TRUE, GTF.featureType="gene", GTF.attrType="gene_id", 
                        isPairedEnd = TRUE)

# comparison count table (wildtype male full genome)
wt_m_e <- featureCounts(files="datafiles/wildtype_male.bam",
                        annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                        isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType="gene_id", 
                        isPairedEnd = TRUE)
wt_m_g <- featureCounts(files="datafiles/wildtype_male.bam",
                        annot.ext="datafiles/Mus_musculus.GRCm38.98.gtf",
                        isGTFAnnotationFile = TRUE, GTF.featureType="gene", GTF.attrType="gene_id", 
                        isPairedEnd = TRUE)

# combine count matrices
mt_f_df <- combine_counts(mt_f_e$counts, mt_f_g$counts)
wt_f_df <- combine_counts(wt_f_e$counts, wt_f_g$counts)
wt_m_df <- combine_counts(wt_m_e$counts, wt_m_g$counts)

# build df for grouped bar plot
group_df <- data.frame(Gene_ID = character(), condition = character(), value = integer())

group_df <- rbind(group_df, data.frame("sample" = "mt female", "condition" = "exons", 
                                       "value" = sum(mt_f_df$exons)))
group_df <- rbind(group_df, data.frame("sample" = "mt female", "condition" = "genes", 
                                       "value" = sum(mt_f_df$genes)))
group_df <- rbind(group_df, data.frame("sample" = "mt female", "condition" = "introns", 
                                       "value" = sum(mt_f_df$introns)))

group_df <- rbind(group_df, data.frame("sample" = "wt female", "condition" = "exons", 
                                       "value" = sum(wt_f_df$exons)))
group_df <- rbind(group_df, data.frame("sample" = "wt female", "condition" = "genes", 
                                       "value" = sum(wt_f_df$genes)))
group_df <- rbind(group_df, data.frame("sample" = "wt female", "condition" = "introns", 
                                       "value" = sum(wt_f_df$introns)))

group_df <- rbind(group_df, data.frame("sample" = "wt male", "condition" = "exons", 
                                       "value" = sum(wt_m_df$exons)))
group_df <- rbind(group_df, data.frame("sample" = "wt male", "condition" = "genes", 
                                       "value" = sum(wt_m_df$genes)))
group_df <- rbind(group_df, data.frame("sample" = "wt male", "condition" = "introns", 
                                       "value" = sum(wt_m_df$introns)))
# grouped plotting
ggplot(group_df, aes(fill=condition, y=value, x=sample)) + geom_bar(position="dodge", stat="identity")

# stacked plotting 
ggplot(group_df, aes(fill=condition, y=value, x=sample)) + geom_bar(position="stack", stat="identity")

# exon vs intron scatter plotting - mt female
ggplot(mt_f_df, aes(x=exons, y=introns)) + geom_point()
ggplot(mt_f_df, aes(x=exons, y=introns)) + geom_point() + geom_text(label=mt_f_df$Gene_ID) 

# exon vs intron scatter plotting - wt female
ggplot(wt_f_df, aes(x=exons, y=introns)) + geom_point()
ggplot(wt_f_df, aes(x=exons, y=introns)) + geom_point() + geom_text(label=wt_f_df$Gene_ID)

# exon vs intron scatter plotting - wt male
ggplot(wt_m_df, aes(x=exons, y=introns)) + geom_point()
ggplot(wt_m_df, aes(x=exons, y=introns)) + geom_point() + geom_text(label=wt_m_df$Gene_ID)
