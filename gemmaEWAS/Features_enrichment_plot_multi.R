## Author: Dario Galanti April 2021
## Aim: Plot enrichment of variants in different genomic features for increasing -log(p) similarly to Atwell et al. 2010.
## Input: txt files with enrichment analysis for increasing -log(p) values for several phenotypes

library(ggplot2)
library(dplyr)

Features_colors <- c(
  #"global" = "grey40",
  #"genes" = ,
  "CDS" = "#D06465",
  "introns" = "#BF6B92",
  "TEs" = "#698838",
  "prom" = "#CB6E34",
  "interg" = "#7BA5D6")


file <- basename("Exoreads_enrichment_features_EWAS")
file <- substr(file, 1, nchar(file)-5)
outdir <- dirname(".")

## 1) READ DATA AND TIDY UP
enrichment1 <- read.table("ErysiRes_noout_noConv1Out/Cs_enrichment_features_EWAS.txt", header=TRUE, sep="\t")
enrichment2 <- read.table("MpersicaeRes_noout2_noConv1Out/Cs_enrichment_features_EWAS.txt", header=TRUE, sep="\t")

enrichment1$ReadGroup <- "ErysiRes"
enrichment2$ReadGroup <- "MpersicaeRes"

enrichment <- rbind(enrichment1, enrichment2)
enrichment$ReadGroup <- factor(enrichment$ReadGroup, levels=c("MpersicaeRes","ErysiRes"))

## OPTIONAL: Remove lines with observed frequency = 0
#enrichment <- enrichment[!(enrichment$Observed_freq == 0),]

## Rename features
enrichment$Feature <- gsub("intergenic","interg",enrichment$Feature)
enrichment$Feature <- gsub("promoters","prom",enrichment$Feature)

## REMOVE GENES AND GLOBAL
enrichment <- enrichment[!enrichment$Feature == "global",]
enrichment <- enrichment[!enrichment$Feature == "genes",]

## 3a)  PRINT ENRICHMENT AND FDR PLOT
max_enrich <- ceiling(max(enrichment$Enrichment, 8))


pdf(paste(outdir,"/",file,".pdf",sep=""), width = 9, height = 4.3)
ggplot(enrichment, aes(x=logP, y=Enrichment, color=Feature)) +
  geom_line(size=1.3) +
  scale_x_continuous(name = "-log(p)") +
  scale_y_continuous(name = "Enrichment", limits = c(0,max_enrich) ) + 
  scale_color_manual(values=Features_colors) +
  facet_wrap(~ ReadGroup) +
  theme_light() +
  #guides(color = guide_legend(override.aes = list(size = 1.2))) + 
  theme(
    aspect.ratio=1.05,
    axis.title.x = element_text(hjust = 0.5, size=16),
    axis.text.x = element_text(size=14),
    axis.title.y = element_text(hjust = 0.5, size=16),
    axis.text.y = element_text(size=14),
    panel.grid = element_line(color="#E6E6E6"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
    strip.background = element_rect(fill="white", color="grey20"),
    panel.background = element_rect(fill="white"), strip.text = element_text(color="grey20", size = 13)
  )
dev.off()







