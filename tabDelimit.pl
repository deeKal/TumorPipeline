#Perl script for parsing vcf files to tab delimited files
#!/bin/perl

use warnings;
use strict;
use POSIX;

my $line;
my $file = $ARGV[0];
my $outFile = substr($file, 0, -3)."final.tsv";

open(FILE,"<",$file) or die "$!\n";

open (OUTFILE, ">>",$outFile) or die $!; 	


while ($line = <FILE>) {
	$line =~ tr/\|/\t/;

	$line =~ s/CSQ/Allele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tPICK\tVARIANT_CLASS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tTSL\tAPPRIS\tCCDS\tENSP\tSWISSPROT\tTREMBL\tUNIPARC\tREFSEQ_MATCH\tSOURCE\tGENE_PHENO\tSIFT\tPolyPhen\tDOMAINS\tHGVS_OFFSET\tAF\tAFR_AF\tAMR_AF\tEAS_AF\tEUR_AF\tSAS_AF\tAA_AF\tEA_AF\tExAC_AF\tExAC_Adj_AF\tExAC_AFR_AF\tExAC_AMR_AF\tExAC_EAS_AF\tExAC_FIN_AF\tExAC_NFE_AF\tExAC_OTH_AF\tExAC_SAS_AF\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE/;

	print OUTFILE $line;

}

