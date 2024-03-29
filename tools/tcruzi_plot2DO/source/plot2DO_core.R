suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(pracma)
  library(doParallel)
  library(foreach)
})


ComputeNormalizationFactors <- function(reads)  {

  # Compute rescaling coefficients, such that after rescaling the average occupancy is 1, for each chromosome 
  occ <- coverage(reads)
  chrLabel <- seqlevels(occ)
  coverageWeight <- lapply(chrLabel, function(chr) 1/mean(occ[[chr]]))
  names(coverageWeight) <- chrLabel
  return(coverageWeight)
}


noCores <- switch(Sys.info()[['sysname']],
                  Windows = 1, # the parallel functions do not work in Windows...
                  Linux   = floor(detectCores(logical = FALSE) / 2), # logical = FALSE does not work in linux..., so we divide by 2
                  Darwin  = detectCores(logical = FALSE)) # Mac

if(noCores > 1) {
  noCores <- min(noCores - 1, 8) # Do not use more than 8 cores (to prevent memory issues)
}

# reads <- resized_reads
# coverageWeight
# referenceGRanges, 
# referenceGenome <- params$genome

ComputeCoverageMatrixGenes <- function(reads, coverageWeight, 
                                       referenceGRanges, referenceGenome)
{
  # For human and mouse data, limit the number of cores to 2 (increase this # according to the available memory that is available on your system)
  if (referenceGenome %in% c('mm9', 'mm10', 'hg18', 'hg19', 'hg38')){
    noCores = min(c(2, noCores)) 
  }
  
  occ <- coverage(reads)
  chrLabel <- seqlevels(occ)
  noChr <- length(chrLabel)
  
  registerDoParallel(cores=noCores)
  
  
  #referenceGRanges_subset <- referenceGRanges[seq(1, 15)]
  #referenceGRanges <- referenceGRanges_subset
  
  # For each fragment size, compute the corresponding coverage/occupancy
  #occMatrixTranspose <- foreach(l=lMin:lMax, .combine='cbind') %dopar% {
  occMatrixTranspose <- foreach(region_index = seq(1, length(referenceGRanges)), 
                                .export=c("AlignRegionsTranspose"), 
                                .packages=c('GenomicRanges', 'pracma'), 
                                .combine='cbind') %dopar% {
                                  
                                  # Keep only the reads with the specific length l
                                  #goodReadsInd <- (readLength == l)
                                  #goodReads <- reads[goodReadsInd]
                                  
                                  gene_region <- referenceGRanges[region_index]   
                                  # TODO: ver si esta bien considerar que un gen es un registro de los utrme
                                  # Keep only the reads within the specified region
                                  goodReadsInd <- findOverlaps(reads, gene_region, 
                                                               type = "within", 
                                                               ignore.strand = FALSE) 
                                  goodReads <- reads[queryHits(goodReadsInd)]
                                  
                                  #chr_name <- unique(seqnames(goodReads))
                                  
                                  # Compute average occupancy
                                  occ <- coverage(goodReads, weight = coverageWeight)
                                  
                                  #occ_region <- occ[[chr_name]]
                                  
                                  alignedRegionsTranspose <- AlignRegionsTranspose(occ, referenceGRanges)
                                  
                                  # mean by gene (or utr):
                                  # sort_by <- colMeans(alignedRegionsTranspose) 
                                  
                                  # NO: mean by position:
                                  # rowMeans(alignedRegionsTranspose)
                                  alignedRegionsTranspose[, region_index]
                                }
  occMatrix <- t(occMatrixTranspose)
  rownames(occMatrix) <- c()
  
  stopImplicitCluster()
  #registerDoSEQ()
  
  # sort by mean gene occupancy:
  gene_mean_occ <- rowMeans(occMatrix)
  order_indexes <- order(gene_mean_occ)
  result <- occMatrix[order_indexes, ]
  
  return(result)
}



# Parallelized version, much faster (~100x faster on my macbook pro)
ComputeCoverageMatrix <- function(lMin, lMax, beforeRef, afterRef, reads, 
                                  coverageWeight, referenceGRanges, readLength, referenceGenome)
{
  # For human and mouse data, limit the number of cores to 2 (increase this # according to the available memory that is available on your system)
  if (referenceGenome %in% c('mm9', 'mm10', 'hg18', 'hg19', 'hg38')){
    noCores = min(c(2, noCores)) 
  }
  
  occ <- coverage(reads)
  chrLabel <- seqlevels(occ)
  noChr <- length(chrLabel)
  
  registerDoParallel(cores=noCores)
  
  # For each fragment size, compute the corresponding coverage/occupancy
  #occMatrixTranspose <- foreach(l=lMin:lMax, .combine='cbind') %dopar% {
  occMatrixTranspose <- foreach(l=lMin:lMax, 
                                .export=c("AlignRegionsTranspose"), 
                                .packages=c('GenomicRanges', 'pracma'), 
                                .combine='cbind') %dopar% {
    # Keep only the reads with the specific length l
    goodReadsInd <- (readLength == l)
    goodReads <- reads[goodReadsInd]
    
    # Compute average occupancy
    occ <- coverage(goodReads, weight=coverageWeight)
    alignedRegionsTranspose <- AlignRegionsTranspose(occ, referenceGRanges)
    rowMeans(alignedRegionsTranspose)
  }
  occMatrix <- t(occMatrixTranspose)
  rownames(occMatrix) <- c()
  
  stopImplicitCluster()
  #registerDoSEQ()
  
  return(occMatrix)
}

AlignRegionsTranspose <- function(profile, referenceGRanges)
{
  # Create Views with all the referenceGRanges
  chrName <- unique(as.character(seqnames(referenceGRanges)))
  myViews <- Views(profile[chrName], as(referenceGRanges, "IntegerRangesList")[chrName])
  
  alignedProfilesList <- lapply(myViews, function(gr) viewApply(gr, as.vector))
  alignedProfiles <- do.call("cbind", alignedProfilesList)
  
  ## Get the index of referenceGRanges, which were reorganized by as(referenceGRanges, "IntegerRangesList")
  listInd <- split(1:length(referenceGRanges), as.factor(seqnames(referenceGRanges)))
  idx <- do.call("c", listInd)
  
  alignedProfiles <- alignedProfiles[, order(idx)]
  
  ## Flip regions from the Crick strand
  crickInd <- which(strand(referenceGRanges) == "-")
  alignedProfiles[ , crickInd] <- flipud(alignedProfiles[ , crickInd])
  
  return(alignedProfiles)
}

ConstructReferenceGRanges <- function(optSites, annotat, selectedReference, beforeRef, afterRef, genome, optAlign) {
  
  chrLen <- annotat$chrLen

  annotations <- annotat$annotations
  
  if (is.null(optSites)) {
    
    if(selectedReference == "Plus1" & genome != 'sacCer3') {
      print_help(opt_parser)
      stop("Plus1 annotations are available only for the sacCer3 genome. For other genomes, the available options for --reference are: TSS (default), TTS", call.=FALSE)
    } 
    
    referencePos <- annotations[[selectedReference]] 
    referenceChr <- annotations$Chr
    refStrand <- annotations$Strand
    
    validValues <- (! is.nan(referencePos))
    referencePos <- referencePos[validValues]
    referenceChr <- referenceChr[validValues]
    refStrand <- refStrand[validValues]
    
    # TODO: fix for strand + or - instead of +1 -1 ??
    watson <- (refStrand == 1)
    
  } else {
    
    referenceFilename <- optSites
    referenceFilePath <- file.path(annotationsBasePath, referenceFilename)
    #referenceFilePath <- optSites
    sites <- import.bed(referenceFilePath)
    
  if(genome == 'tair10') {
    seqlevelsStyle(sites) <- 'Ensembl'
  } else if (genome == 'TcruziCLBrenerEsmeraldo'){ 
    # do nothing or check tritryp style
  } else if (genome == 'Sylvio'){ 
    # do nothing or check tritryp style
  } else {
    # Make sure the chromosome names are following the 'UCSC' convention
    seqlevelsStyle(sites) <- 'UCSC'
  }
    
    switch(optAlign, 
           "center"    ={ fixAlign = "center" },
           "fivePrime" ={ fixAlign = "start" },
           "threePrime"={ fixAlign = "end" })
    
    referenceGRanges <- resize(sites, width = 1, fix = fixAlign)
    
    referencePos <- start(referenceGRanges)
    referenceChr <- seqnames(referenceGRanges)
    refStrand <- strand(referenceGRanges)
    
    # TODO: fix for strand + or - instead of +1 -1 ??
    watson <- as.vector(strand(refStrand) != '-')      # ('+' or '*')
  }
  
  leftEdge <- referencePos - beforeRef
  rightEdge <- referencePos + afterRef
  
  leftEdgeCrick <- referencePos - afterRef
  rightEdgeCrick <- referencePos + beforeRef
  
  leftEdge[!watson] <- leftEdgeCrick[!watson]
  rightEdge[!watson] <- rightEdgeCrick[!watson]
  
  wholeChr <- GRanges(seqnames=names(chrLen), IRanges(start=rep(1, length(chrLen)), end=chrLen))  
  
  referenceGRanges <- GRanges(seqnames=referenceChr, IRanges(start=leftEdge, end=rightEdge), strand=refStrand)
  result <- subsetByOverlaps(referenceGRanges, wholeChr, type="within", ignore.strand=TRUE)
  
  return(result)
  
}

CalculatePlotData <- function(params, reads, referenceGRanges, genes_analysis = FALSE) {
  
  # Compute the histogram
  readLength <- width(reads) 
  h <- hist(readLength, breaks=seq(from = 0.5, to = max(readLength, 1000, params$lMax)+0.5, by = 1), plot=FALSE)
  lengthHist <- 100 * h$density[params$lMin:params$lMax]
  
  totalNoReads <- length(reads)
  
  resized_reads <- switch(params$plotType, 
                      "OCC"            = reads,
                      "DYADS"          = resize(reads, width = 1, fix = "center"),
                      "FIVEPRIME_ENDS" = resize(reads, width = 1, fix = "start"),
                      "THREEPRIME_ENDS"= resize(reads, width = 1, fix = "end"))
  
  coverageWeight <- ComputeNormalizationFactors(resized_reads)
  
  outputFilePath <- GetOutputMatrixFilePath(params$plotType, params$referencePointsBed, 
                                            params$reference, params$siteLabel, 
                                            params$lMin, params$lMax, params$sampleName)
  
  if (genes_analysis) {
    occMatrix <- ComputeCoverageMatrixGenes(resized_reads, 
                                            coverageWeight, referenceGRanges, 
                                            params$genome)
    
    genes_out_base_path <- dirname(outputFilePath)
    genes_file_name <- basename(outputFilePath)
    genes_file_name <- paste0('genes_', genes_file_name)
    outputFilePath <- file.path(genes_out_base_path, genes_file_name)
    
  } else {
    occMatrix <- ComputeCoverageMatrix(params$lMin, params$lMax, 
                                       params$beforeRef, params$afterRef, 
                                       resized_reads, coverageWeight, 
                                       referenceGRanges, readLength, 
                                       params$genome)
    
  }
  
  print(outputFilePath)
  
  # Parameters to be saved:
  lMin <- params$lMin
  lMax <- params$lMax
  beforeRef <- params$beforeRef
  afterRef <- params$afterRef
  sampleName <- params$sampleName
  siteLabel <- params$siteLabel
  save(occMatrix, lMin, lMax, beforeRef, afterRef, totalNoReads, lengthHist, sampleName, siteLabel, file=outputFilePath)
  
}

