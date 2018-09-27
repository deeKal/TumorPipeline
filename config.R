rundir <- "/home/dnalab/MolecularDiagnosticsLab/DataForAnalysis/QiaSeq4_Test"  #Need to be changed each time
umi <- TRUE  #Run a UMI analysis??? ALWAYS FALSE on Devyser panel
singlEnd <- FALSE  #Paired end or single end sequencing? by default paired end
RG <- "G13GV.1"  #Need to be changed each time (Flowcell code.1)
lab <- "MLCDiag"  #Sequencing Center

reference <- "/media/data1/ReferenceGenomes/hg19.fa"  #Path to reference file
QiaSeqTargets <- "/home/dnalab/MolecularDiagnosticsLab/QiaSeqBEDs/DHS-102Z.primers-150bp.bed"
dvsrTargets <- "/home/dnalab/MolecularDiagnosticsLab/DevyserBEDs/BED_BRCA_hg19_2018_05_08.bed"
vepReference <- "~/.vep/homo_sapiens/89_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"

targets <- QiaSeqTargets
if (singlEnd) {
    targets <- dvsrTargets
}

pathToPicard <- "/home/dnalab/Tools/picard/build/libs/picard.jar"
pathToBWA <- "/home/dnalab/Tools/BWA/bwa"
pathToGatk <- "/home/dnalab/Tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"
pathToVep <- "/home/dnalab/MolecularDiagnosticsLab/Annotation/ensembl-vep/vep"
tabDelimitScript <- "/media/data1/TumorPipeline/tumorPipeline/tabDelimit.pl"

dbSNPSites <- "/media/data1/KnownSites/dbsnp_138.hg19.vcf.gz"
ThGSites <- "/media/data1/KnownSites/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
cosmicSites <- "/media/data1/KnownSites/CosmicCodingMuts.vcf.gz"
knownSites <- paste(" -knownSites ", c(dbSNPSites, ThGSites), sep = " ", collapse = " ")

covThr <- 50  #Threshold for coverage
mapQThr <- 10  #Threshold for mapping quality

fqdir <- paste(rundir, "/Fastq", sep = "")

files <- list.files(fqdir, full.names = FALSE, pattern = "1_001.fastq.gz")

