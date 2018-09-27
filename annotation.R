annotateVCF <- function(sample, pathToVep, tabDelimitScript, vepReference) {
    
    library(vcfR)
    library(tidyr)
    library(GenomeInfoDb)
    library(VariantAnnotation)
    library(xlsx)
    
    vcfdir <- "/Vcf"
    vcfFileM <- paste(vcfdir, sample, ".Mutect2.vcf", sep = "")
    
    vcfFile <- paste(vcfdir, sample, ".vcf", sep = "")
    
    file <- paste(vcfdir, sample, ".annotated.vcf", sep = "")
    
    myFile <- paste(vcfdir, sample, ".annotated.tsv", sep = "")
    finalFile <- paste(vcfdir, sample, ".annotated.final.tsv", sep = "")
    excelFile <- paste(vcfdir, sample, ".xlsx", sep = "")
    
    vcf <- readVcf(vcfFileM)
    seqlevelsStyle(vcf) <- "Ensembl"
    writeVcf(vcf, vcfFile)
    
    # Connect to VEP Server and annotate file
    command <- paste("perl ", pathToVep, " --cache -i ", vcfFile, " -o ", file, 
	" --vcf --everything --fork 6 --assembly GRCh37 --format vcf  --cache \
	--variant_class --humdiv --gene_phenotype --regulatory  --numbers --domains \
	--hgvs --protein --symbol --ccds --uniprot --tsl --appris --canonical \
	--biotype --af --af_1kg --af_exac   --pubmed --merged --port 3337 \
	-check_existing -sift b -polyphen p --offline --fasta ", 
        vepReference, " --no_stats --chr 13,17 --flag_pick")
    system(command, intern = TRUE, wait = TRUE)
    
    vcf.frame <- vcfR2tidy(vcf, info_only = FALSE, single_frame = TRUE, toss_INFO_column = TRUE)
    vcf.dat <- vcf.frame$dat
    
    
    vcf.dat <- subset(vcf.dat, select = -c(DB, NLOD, PON, RPA, RU, Indiv, gt_GQ, 
        gt_PGT, gt_PID, gt_PL))
    vcf.dat <- vcf.dat %>% separate_rows(CSQ, sep = ",")
    vcf.split <- strsplit(vcf.dat$gt_AD, split = ",")
    
    for (i in 1:length(vcf.split)) {
        if (length(vcf.split[[i]]) == 3) {
            vcf.split[[i]] <- vcf.split[[i]][-3]
        }
    }
    
    l <- as.data.frame(vcf.split, row.names = c("Ref", "Mut"))
    l <- t(l)
    l <- as.data.frame(l)
    l$Mut <- as.numeric(as.character(l$Mut))
    l$Ref <- as.numeric(as.character(l$Ref))
    vcf.dat$gt_DP <- l$Ref + l$Mut
    
    vcf.dat <- vcf.dat[c(1:7, 14:16, 19, 8:13, 17, 18, 20:25)]
    
    
    write.table(vcf.dat, file = myFile, sep = "\t", col.names = TRUE, 
		row.names = FALSE, quote = FALSE)
    
    
    command <- paste("perl ", tabDelimitScript, myFile)
    system(command, intern = TRUE, wait = TRUE)
    
    
    newdata <- read.table(finalFile, header = TRUE, sep = "\t")
    newdata <- subset(newdata, Feature == "NM_000059.3" | Feature == "NM_007294.3")
    write.xlsx(newdata, excelFile, row.names = FALSE, showNA = FALSE)
    
    file.remove(vcfFile)
    file.remove(file)
    file.remove(myFile)
    file.remove(finalFile)
    
    detach("package:vcfR")
    detach("package:tidyr")
    detach("package:GenomeInfoDb")
    detach("package:VariantAnnotation")
    detach("package:xlsx")
    
}

