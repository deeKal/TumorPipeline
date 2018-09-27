wd <- "/media/data1/TumorPipeline/tumorPipeline"
setwd(wd)

source("config.R")
source("fun.R")
source("coverageFun.R")
source("annotation.R")

setwd(rundir)

library("fastqcr")
metricsdir <- "/Metrics"
dir.create(file.path(metricsdir), showWarnings = FALSE)

fastqc(fqdir, metricsdir)

detach("package:fastqcr")

for (file in files) {
    sample <- substr(file, 1, nchar(file) - 21)
    
    picardFastqToSam(sample, pathToPicard, RG, lab, singlEnd, umi)
    
    picardMarkIlluminaAdapters(sample, pathToPicard)
    
    bwaAlignment(sample, reference, pathToPicard, pathToBWA)
    
    gatkBaseRecalibration(sample, reference, pathToGatk, targets, knownSites)
    
    coverage(sample, targets, covThr, mapQThr)
    
    mutect2VariantCalling(sample, reference, pathToGatk, targets, dbSNPSites, cosmicSites)
    
    annotateVCF(sample, pathToVep, tabDelimitScript, vepReference)
}
