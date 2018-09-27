coverage <- function(sample, targetFile, covThr, mapQThr) {
    library(GenomicFeatures)
    library(GenomicRanges)
    library(GenomicAlignments)
    library(Rsamtools)
    
    get.start <- function(pos, targetPos) {
        p <- pmax(pos, targetPos)
        return(p)
    }
    get.end <- function(pos, targetPos) {
        p <- pmin(pos, targetPos)
        return(p)
    }
    
    
    bamdir <- "/Bam"
    covdir <- "/Coverage"
    
    file <- paste(bamdir, sample, ".recal.bam", sep = "")
    covFile <- paste(covdir, sample, ".Low.Coverage.bed", sep = "")
    
    dir.create(file.path(covdir), showWarnings = FALSE)
    
    targetData <- read.table(targetFile, header = F)
    colnames(targetData) <- c("chr", "start", "end", "id", "score", "strand")
    targetGrObject <- with(targetData, GRanges(chr, IRanges(start + 1, end)))
    seqlevelsStyle(targetGrObject) <- "UCSC"
    
    # make vector with ranges from bed file in order to compare with ranges of
    # undercovered regions
    ranges <- cbind(as.character(targetData$chr), rep(":", length(targetData$chr)), 
        as.character(targetData$start + 1), rep("-", length(targetData$chr)), as.character(targetData$end))
    
    ranges <- apply(ranges, 1, paste, collapse = "")
    
    flag <- scanBamFlag(isNotPassingQualityControls = FALSE)
    
    param <- ScanBamParam(which = targetGrObject, flag = flag, mapqFilter = mapQThr)
    
    
    reads <- readGAlignments(bamFile, param = param)
    
    cvg <- coverage(reads)  #total no of reads
    cvg.i <- slice(cvg, upper = covThr, lower = 0)
    
    cvg.ranges <- reduce(as(cvg.i, "GRanges"))
    
    fo <- findOverlaps(cvg.ranges, targetGrObject)
    
    qAreas <- cvg.ranges[queryHits(fo)]
    sAreas <- targetData[subjectHits(fo), ]
    
    df <- data.frame(seqnames = seqnames(qAreas), 
					starts = get.start(start(qAreas) - 1, sAreas$start), 
					ends = get.end(end(qAreas), sAreas$end), 
					names = c(rep("LowCoverage", length(qAreas))), 
					strands = strand(qAreas), 
					targetseqnames = sAreas$chr, 
        			targetstarts = sAreas$start, 
					targetends = sAreas$end, 
					targetname = sAreas$id)
    
    df$scores <- df$ends - df$starts
    
    
    write.table(df, file = covFile, sep = "\t", col.names = FALSE, row.names = FALSE, 
        quote = FALSE)
    
    detach("package:GenomicFeatures")
    detach("package:GenomicRanges")
    detach("package:GenomicAlignments")
    detach("package:Rsamtools")
}
