## Author: Dario Galanti
## INPUT: ewas output from GEMMA
## OUTPUT: Manhattan and qqplot in png
## RUN: Rscript --vanilla qqman_GEMMA.R SupFamTIPs_log_1out_RIR/SupFamTIPs_log_1out_RIR.assoc.txt

library(qqman)
library(fastman)  # For plot with larger lables
library(data.table)

## OPTIONAL: Define number of independent tests to caldulate the significance threshold.
## In GWAS this is done by SNP pruning with Plink according to the study
## "Addressing population-specific multiple testing burdens in genetic association studies"
## In EWAS we can use the number of DMRs or use something similar to SNP pruning to count independent cytosines.
## If "tests" is not defined, the total number of cytosines will be used
#tests <- 548   # 548 for Erysi DMRs or 162 for Mpersicae DMRs

args <- commandArgs(trailingOnly = TRUE)

Results <- fread(args[1], header=T, data.table = F, sep="\t")
file <- basename(args[1])
file <- substr(file, 1, nchar(file)-10)
outdir <- dirname(args[1])

## FORMAT FOR QQMAN
Results <- Results[c(2,1,3,12)]
colnames(Results) <- c("SNP","CHR","BP","P")
Results$CHR <- as.integer(gsub("Scaffold_","",Results$CHR))

# MANHATTAN PLOT
topval <- -log10(min(Results$P)) #save top value for manhattan ylim scaling
"##CALCULATE BONFERRONI THRESHOLD"
SNPs <- nrow(Results)
if (exists("tests")==FALSE) {tests <- SNPs}
Bonferroni <- -log10(0.05/tests)

"##CALCULATE FDR TRESHOLD"
pvec <- sort(Results$P)
FDRtr <- 6.0
for (i in 1:length(pvec)) {
  P <- pvec[i]
  cv <- (i/SNPs)*0.2
  if(P < cv) { FDRtr <- -log10(P)
  }
}
print(paste("Bonferroni treshold for",tests,"indipendent tests: ", Bonferroni))
print(paste("FDR 0.20 treshold for",SNPs,"SNPs:",FDRtr))

##MAKING MANHATTAN
ytop <- max(Bonferroni, topval)
## "cex" regulates dot size; "cex.axes" axes numbers; "cex.lab" axes lables
man.pal <- c("darkblue","darkcyan") # NEW

## FASTMAN WITH LARGER LABLES
png(paste(outdir,"/FDR_0.2_",file,"_manhattan_GEMMA.png",sep=""), width = 1400, height = 600)
par(mar = c(5, 5.4, 4, 2) + 0.1) # Add space for the y axes label. c(bottom, left, top, right) margins with default c(5,4,4,2)+0.1
par(mgp=c(3.1,1.1,0)) #Add space between axes lables and plot. c(axis.title, axis.label, axis.line) with default c(3,1,0)
fastman(Results, cex=1.6, cex.axis=2.1, cex.lab=2.1, col = man.pal, main=file, maxP=NULL,
          genomewideline = Bonferroni,
          #suggestiveline = FDRtr,
          ylim = c(0, ceiling(ytop)))
dev.off()

## MAKING QQPLOT
png(paste(outdir,"/",file,"_qqplot_GEMMA.png",sep=""), width = 500, height = 500)
par(mar = c(5, 5.2, 4, 2) + 0.1) # Add space for the y axes lable
qq(Results$P, cex = 1.2, cex.axis = 1.7, cex.lab=1.8, col="darkblue")
dev.off()

print("#########")
print("## END QQMAN ##")
print
