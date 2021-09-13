#!/usr/bin/Rscript
library("optparse")

options = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Input BAM file"),
  make_option(c("-t", "--type"), type="character", default="occ", 
              help="Types of profiles to be computed, separated by commas only [options: occ, dyads; default = %default]"),
  make_option(c("-l", "--minLength"), type="integer", default=50, 
              help="The smallest DNA fragment to be considered [default = %default]"),
  make_option(c("-L", "--maxLength"), type="integer", default=200, 
              help="The largest DNA fragment to be considered [default = %default]"),
  make_option(c("-O", "--outputDir"), type="character", default="./", 
              help="The largest DNA fragment to be considered [default = %default]")  
) 

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least the dataset file name must be supplied.", call.=FALSE)
}

# Load the necessary R packages
suppressPackageStartupMessages({
  library(Rsamtools)
  library(rtracklayer)
  library(tools)
})

##################
# Initialization #
##################
inputFilename <- opt$file # input filename
plotType <- toupper(opt$type) # type of profile to construct (occ or dyads)
Lmin <- opt$minLength  # the smallest DNA fragment to be considered
Lmax <- opt$maxLength  # the largest DNA fragment to be considered
outputDir <- opt$outputDir

# Import paired-end reads
#sample.name <- sub(".bam", "", inputFilename)
sample.name <- sub(".bam", "", basename(inputFilename))
all_fields <- c("rname", "pos", "isize")
param <- ScanBamParam(what = all_fields, 
                     flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, 
                                        isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                                        isMinusStrand = FALSE, isMateMinusStrand = TRUE,
                                        isNotPassingQualityControls = FALSE))
bam = scanBam(inputFilename, param=param)

# Keep only the proper reads, with the length > 0
posStrandReads = (bam[[1]]$isize > 0)

reads = GRanges(seqnames=Rle(bam[[1]]$rname[posStrandReads]),
                ranges = IRanges(start=bam[[1]]$pos[posStrandReads], width=bam[[1]]$isize[posStrandReads]),
                strand = "*")
rm(bam)

# Discard the reads from the yeast mitochondrial DNA
# reads = reads[seqnames(reads) != 'chrM']

# Size selection
readLength = width(reads)
goodReadsInd = ((readLength >= Lmin) & (readLength <= Lmax))
reads = reads[goodReadsInd]

TotalNoReads = length(reads)

#################
# Normalization #
#################
# Compute rescaling coefficients, such that after rescaling the average occupancy is 1, for each chromosome

Occ = coverage(reads)
chrLabel = seqlevels(Occ)
noChr = length(chrLabel)
coverageWeight = list()
for(chr in chrLabel)
{
  coverageWeight[[chr]] = 1/mean(Occ[[chr]])
}

switch(plotType, 
       OCC={
         # Compute occupancy
         Occ <- coverage(reads, weight=coverageWeight)
         
         outFileName <- paste("Normalized_Occ.", sample.name, ".", Lmin, "-", Lmax, ".bw",sep="")
         outFilePath <- file.path(outputDir, outFileName)
         
         # Save bigwig files
         export.bw(Occ, outFilePath)
       },
       DYADS={
         # Compute Dyads
         Dyads <- coverage(resize(reads, width = 1, fix = "center"), weight=coverageWeight)
         
         outFileName <- paste("Normalized_Dyads.", sample.name, ".", Lmin, "-", Lmax, ".bw",sep="")
         outFilePath <- file.path(outputDir, outFileName)
         
         # Save bigwig files
         export.bw(Dyads, outFilePath)
       }
)