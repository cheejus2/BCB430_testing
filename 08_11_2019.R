# resource: http://bioconductor.org/help/course-materials/2015/BioC2015/bioc2015rnaseq.html
# resource: https://bioinformatics-core-shared-training.github.io/RNAseq-R/rna-seq-annotation-visualisation.nb.html
# resource: https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf

library("GenomicFeatures")
library("Gviz")
library("GenomicRanges")
library("GenomicAlignments")
library("rtracklayer")

# making txdb
txdb <- makeTxDbFromGFF("datafiles/Mus_musculus.GRCm38.98.gtf", format="gtf")

# setting genomic range
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

# visualize sashimi plot
altrack_wt_f_s <- AlignmentsTrack("datafiles/wildtype_female_Grin1.bam", type=c("coverage","sashimi"))
altrack_wt_m_s <- AlignmentsTrack("datafiles/wildtype_male_Grin1.bam", type=c("coverage","sashimi"))
altrack_mt_f_s <- AlignmentsTrack("datafiles/mutant_female_Grin1.bam", type=c("coverage","sashimi"))
plotTracks(list(gtrack, grtrack, altrack_wt_f_s), from=25291177, to=25319253, chromosome="2")
plotTracks(list(gtrack, grtrack, altrack_wt_m_s), from=25291177, to=25319253, chromosome="2")
plotTracks(list(gtrack, grtrack, altrack_mt_f_s), from=25291177, to=25319253, chromosome="2")

# determining features contained within the gtf file
gtf <- rtracklayer::import('datafiles/Mus_musculus.GRCm38.98.gtf')
gtf_df = as.data.frame(gtf)
feature_list <- unique(gtf_df[7])

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
                       isPairedEnd = TRUE, allowMultiOverlap = TRUE)
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
