picardFastqToSam <- function(sample, pathToPicard, RG, lab, singlEnd, umi){
	#Convert Fastq files to unaligned bam files (ubam) in order to store metadata
	
	fqdir <- "/Fastq"
	otherbamdir <- "/Bam/Other"
	sample <- substr(file,1,nchar(file)-21)
	ubam  <-  paste(otherbamdir,"/", sample, ".ubam.bam", sep="")	

	dir.create(file.path(otherbamdir), showWarnings = FALSE)

	if(umi){
		file1out <- paste(substr(file,1,nchar(file)-8), "UMI.fastq.gz", sep="")
		file2out <- paste(substr(file2,1,nchar(file)-8), "UMI.fastq.gz", sep="")

		command <- paste("umi_tools extract -I " ,file, 
		" --bc-pattern=NNNNNNNNNNNN --read2-in=" ,file2, 
		" --stdout=" ,file1out, " --read2-out=", file2out, sep="")
		system(command,intern = TRUE, wait=TRUE)

		file1 <- file1out
		file2 <- file2out
	}
	
	if(singlEnd){
		command <- paste("java -Xmx8G -jar ", pathToPicard, " FastqToSam FASTQ=", 
		file, " OUTPUT=", ubam, " READ_GROUP_NAME=", RG, " SAMPLE_NAME=", sample,
		" LIBRARY_NAME=A PLATFORM_UNIT=", RG,
		" PLATFORM=illumina SEQUENCING_CENTER=", lab, sep="")
	} else{
		file2 <- paste(fqdir,"/",sample,"_L001_R2_001.fastq.gz", sep="")
		command <- paste("java -Xmx8G -jar ", pathToPicard, " FastqToSam FASTQ=", 
		file, " FASTQ2=", file2, " OUTPUT=", ubam, " READ_GROUP_NAME=", RG, 
		" SAMPLE_NAME=", sample," LIBRARY_NAME=A PLATFORM_UNIT=", RG,
		" PLATFORM=illumina SEQUENCING_CENTER=", lab, sep="")
	}
	system(command,intern = TRUE, wait=TRUE)

}

picardMarkIlluminaAdapters <- function(sample, pathToPicard){
	#Mark the adapters used by Illumina so that they can be clipped
	#Generate metrics file
	
	otherbamdir <- "/Bam/Other"
	ubam <- paste(otherbamdir,"/", sample, ".ubam.bam", sep="")	
	MAbam <- paste(otherbamdir,"/", sample, ".MA.bam", sep="")
	metricsdir <- "/Metrics"

	command <- paste("java -Xmx8G -jar ", pathToPicard, " MarkIlluminaAdapters I=", 
	ubam," O=", MAbam, " M=", metricsdir, sample,".MarkAdapt.txt", sep="")
	system(command,intern = TRUE, wait=TRUE)
}

bwaAlignment <- function(sample, reference, pathToPicard, pathToBWA) {
	
	otherbamdir <- "/Bam/Other"
	ubam  <-  paste(otherbamdir,"/", sample, ".ubam.bam", sep="")	
	MAbam <- paste(otherbamdir, sample, ".MA.bam", sep="")
	clbam <- paste(otherbamdir, sample, ".clean.bam", sep="")

	command <- "set -o pipefail"
	system(command,intern = TRUE, wait=TRUE)
	#Revert back to Fastq, while clipping adapters|Align to reference genome|Convert SAM file to BAM
    command <- paste("java -Xmx8G -jar ", pathToPicard, " SamToFastq I=", 		
        MAbam, " FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 \
		INTERLEAVE=true NON_PF=true TMP_DIR=temp \
		| /home/dnalab/Tools/BWA/bwa mem -M -t 2 -p ", 							 
        reference, " /dev/stdin \
		| java -Xmx16G -jar ", pathToPicard, 									
		" MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=", 				
        ubam, " OUTPUT=", clbam, " R=", reference, 
		" CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		ATTRIBUTES_TO_RETAIN=XS TMP_DIR=temp", sep = "")
    system(command, intern = TRUE, wait = TRUE)

}


gatkBaseRecalibration <- function(sample, reference, pathToGatk, targets, knownSites){

	otherbamdir <- "/Bam/Other"
	bamdir <- "/Bam"
	metricsdir <- "/Metrics"

	clbam <- paste(otherbamdir, sample, ".clean.bam", sep="")	
	recalBam <- paste(bamdir, sample, ".recal.bam", sep="")
	recalTable<- paste(metricsdir,sample, ".recal_data.table", sep="")
	postRecalTable<- paste(metricsdir,sample, ".post.recal_data.table", sep="")
	recalPlots <- paste(metricsdir, sample,".recalibration_plots.pdf", sep="")
	
	dir.create(file.path(bamdir), showWarnings = FALSE)

	command <- paste("java -jar " ,pathToGatk, " -T BaseRecalibrator -R ",
	reference," -L " ,targets ," -I ", clbam, knownSites," -o ", recalTable, sep="")
	system(command,intern = TRUE, wait=TRUE)

	command <- paste("java -jar " ,pathToGatk, " -T BaseRecalibrator -R ",
	reference," -L " ,targets ," -I ",clbam, knownSites," -BQSR ", recalTable, 
	"  -o ", postRecalTable, sep="")
	system(command,intern = TRUE, wait=TRUE) 

	command <- paste("java -jar " ,pathToGatk, " -T AnalyzeCovariates -R ",reference,
	" -before  ",metricsdir, sample,".recal_data.table -after ", 
	postRecalTable," -plots ", recalPlots, sep="")
	system(command,intern = TRUE, wait=TRUE)

	command <- paste("java -jar " ,pathToGatk, " -T PrintReads -R ",reference,
	" -I ",file," -L " ,targets ," -BQSR ", recalTable, " -o ", recalBam, sep="")
	system(command,intern = TRUE, wait=TRUE)

}



mutect2VariantCalling <- function(sample, reference, pathToGatk, targets, dbSNPSites, cosmicSites){
	
	bamdir <- "/Bam"
	vcfdir <- "/Vcf"
	
	recalBam <- paste(bamdir, sample, ".recal.bam", sep="")
	vcfFile <- paste(vcfdir, sample, ".Mutect2.vcf", sep="")

	dir.create(file.path(vcfdir), showWarnings = FALSE)

	command <- paste("java -jar ", pathToGatk, " -T MuTect2 -R ", reference, 
	" -I:tumor ", recalBam," --dbsnp ", dbSNPSites, " --cosmic ", cosmicSites, 
	" -L ", targets, " -o ", vcfFile, sep="")
	system(command,intern = TRUE, wait=TRUE)
}


