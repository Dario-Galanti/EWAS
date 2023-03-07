## Author: Dario Galanti
## INPUT 1: EWAS results in format EWAS_outdir/${context}/${context}.bedGraph_EWAS_IGVresults.gwas
## INPUT 2: Genome index file: genome.index.fa.fai
## OUTPUT 1: Manhattan plots with different colours for the different contexts and chromosomes separated by space
## OUTPUT 2: qqplots with different colours for the different contexts
## RUN: Rscript --vanilla EWAS_manhattan_contColours.R EWAS_output_phenotype_4download genome.fa.fai n°tests

library(ggplot2)
library(tidyverse)
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

## DEFINE INPUT AND OUTPUTS
EWAS_dir <- args[1] # EWAS output directory for download with results that can be visualized in IGV
index_file <- args[2]
tests <- as.numeric(args[3]) # Define number of independent tests 548 for Erysi DMRs or 162 for Mpersicae DMRs

file <- basename(EWAS_dir)
outdir <- EWAS_dir

intChr_space <- 10000000 # Define space between chromosomes (in bp). 10 MB should work, even less with dense manhattan

CpG <- read.table(paste(EWAS_dir,"/CpG/",file,"_IGVresults.gwas", sep=""), header=TRUE, sep="\t")
CHG <- read.table(paste(EWAS_dir,"/CHG/",file,"_IGVresults.gwas", sep=""), header=TRUE, sep="\t")
CHH <- read.table(paste(EWAS_dir,"/CHH/",file,"_IGVresults.gwas", sep=""), header=TRUE, sep="\t")


## Add context and merge them
CpG$context <- "CG"
CHG$context <- "CHG"
CHH$context <- "CHH"


Results <- rbind(CpG,CHG,CHH)
names(Results) <- c("CHR", "BP", "SNP", "P", "context")  ## Fix column names to be the same as EWAS pipe output
Results$context <- factor(Results$context, levels=c("CG","CHG","CHH"))
Results$Scaff <- as.integer(gsub("Scaffold_","",Results$CHR))
Results <- Results[with(Results, order(Scaff, BP)),]
Results$logP <- -log10(Results$P)

### 1) OPTIONAL: Remove low p.val positions for lighter and faster plotting
Results <- Results[Results$logP>0.10,] # Not good for qqplot


### 2B) Retrieve scaffold length info from genome index and calculate cumulative bp
index <- fread(index_file, header=F)
all_scaff <- unique(Results$CHR)
vec <- index$V1 %in% all_scaff
index <- index[vec,c(1,2)]
# Add inter chromosomal space to chromosomes length
index$V2 <- index$V2 + intChr_space
## This command gives a warning, but it should not be a problem
ind_cum <- as_tibble(cbind(index$V1,c(0,cumsum(index$V2))[1:(length(index$V2))]))
colnames(ind_cum) <- c("CHR","bp_add")
ind_cum$bp_add <- as.numeric(ind_cum$bp_add)


### 3a) Add cumulative bp to SNPs
#colnames(ind_cum) <- c("snp_scaff","snp_bp_add")
Results <- Results %>% 
  inner_join(ind_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)

### 4) Get center cum_bp for each chromosome (to position chr names)
index[1,2] <- index[1,2] - (intChr_space/2) # We remove half of the inter chromosomal space from the first scaffold
chrlimits <- c(0,cumsum(index$V2))
axis_set <- as.data.frame(cbind(index$V1,(index$V2/2),c(0,cumsum(index$V2))[1:(length(index$V2))]))
axis_set$center <- (as.numeric(axis_set$V2) + as.numeric(axis_set$V3))
axis_set <- axis_set[c(1,4)]
colnames(axis_set) <- c("scaff","center")
# If willing to use chr number
axis_set$scaff <- gsub('[a-zA-Z]*_', "", axis_set$scaff)


# CALCULATE BONFERRONI
topval <- max(Results$logP) #save top value for manhattan ylim scaling
#Results$P=10^((-Results$P)) #rrblup results already contain -log(p), convert to p
"##CALCULATE BONFERRONI THRESHOLD"
SNPs <- nrow(Results)
if (exists("tests")==FALSE) {tests <- SNPs}
Bonferroni <- -log10(0.05/tests)


## MAKING MANHATTAN
ytop <- ceiling(max(Bonferroni, topval))
#cont.pantone.pal=c("#CD5B45", "#EEC900", "#1E90FF") # Rust (no real name), Gold2, Dodgerblue
cont.pantone.pal=c("#CD5B45", "#D6B500", "#1E90FF") # Rust (no real name), Gold between 2 and 3, Dodgerblue

chrlimits <- chrlimits[2:(length(chrlimits)-1)]
### MANHATTAN PLOT
pdf(paste(outdir,"/",file,"_manhattan_contCol.pdf",sep=""), width = 10.5, height = 4)  # In pdf (heavier)
ggplot(Results, aes(x = bp_cum, y = logP, color = context, size = logP)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_hline(yintercept = Bonferroni, color = "orangered3", linetype = "dashed") + 
  geom_vline(xintercept = chrlimits, size=0.2, color="grey12", linetype = "dotted") + # Add chromosome vertical limits
  #scale_x_continuous(label = axis_set$scaff, breaks=axis_set$center) +
  scale_x_continuous(label = axis_set$scaff, breaks=axis_set$center, expand = c(0.03,0.03)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ytop)) +
  #scale_fill_manual(values = cont.pantone.pal) +
  scale_color_manual(values = cont.pantone.pal) +
  #scale_size_continuous(range = c(1,1.5)) + # if willing to have bigger dots at the top
  labs(x = NULL, y = "-log(p)") + 
  theme_classic() +
  guides(shape = guide_legend(override.aes = list(size = 2.2))) + 
  theme(
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=15),
    axis.title = element_text(size=16),
    axis.line = element_line(colour = 'grey30', size = 0.1),
    axis.ticks = element_line(colour = "grey30", size = 0.1),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12)
  )
dev.off()



## CONTEXT COLOURED QQPLOT
n <- nrow(Results)
Results <- Results[order(Results$P),]
Results$expected <- -log10(ppoints(n)) # Add expected -logP
Results$clower <- -log10(qbeta(p = 0.05 / 2, shape1 = 1:n, shape2 = n:1)) # Add lower conf int
Results$cupper <- -log10(qbeta(p = 1.95 / 2, shape1 = 1:n, shape2 = n:1)) # Add upper conf int

png(paste(outdir,"/",file,"_qqplot_contCol.png",sep=""), width = 400, height = 360)
ggplot(Results) +
  #geom_ribbon(mapping = aes(x=expected, ymin=clower, ymax=cupper),alpha = 0.1) + # Add ribbon to the expected curve
  geom_point(aes(expected, logP, color = context), alpha = 0.7, size = 2) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5, color="red") +
  xlab("Expected -log(p)") +
  ylab("Observed -log(p)") +
  theme_bw() +
  theme(aspect.ratio=1,
    panel.grid = element_blank(),
    axis.text = element_text(size=16),
    axis.title = element_text(size=16),
    axis.line = element_line(colour = 'grey30', size = 0.1),
    axis.ticks = element_line(colour = "grey30", size = 0.1),
    legend.title = element_text(size=15),
    legend.text = element_text(size=14)
  )
dev.off()





