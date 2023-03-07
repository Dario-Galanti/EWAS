Code for conducting Epigenome Wide Association Studies (EWAS) for multiple phenotypes with [GEMMA](https://github.com/genetics-statistics/GEMMA), starting from methylation unionbed files. <br/>

IMPORTANT:<br/>
- Input unionbed files can be generated and filtered for 1) NAs and 2) Minor Epiallele Frequency (MEF) with code stored [here](https://github.com/Dario-Galanti/WGBS_downstream/tree/main/WGBS_simpleworkflow)<br/>
- Unionbed files with region methylation averages instead of positions can be obtained with code stored [here](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_simpleworkflow/region_meth/average_over_bed.py)<br/>
- To reduce the multiple testing burden of millions of cytosines, one could 1) use region averages instead of positions, 2) use positions but only within predifined regions or 3) Use a strict MEF filtering to only test highly variable positions. The code mentioned above allows all these options. <br/>
- We use [GEMMA](https://github.com/genetics-statistics/GEMMA) to run mixed models correcting for population structure with a distance matrix <br/>
- We use the BIMBAM input files for [GEMMA](https://github.com/genetics-statistics/GEMMA) to provide quantitative methylation values <br/>
- This repository contains scripts for additional plotting and powerful downstream analysis to interpret the results.<br/>
<br/>

SCRIPTS DESCRIPTION:<br/>
<br/>
[gemmaEWA_multi_Epi.sh](https://github.com/Dario-Galanti/EWAS/new/main/gemmaEWAS/gemmaEWA_multi_Epi.sh)<br/>
Run EWAS for multiple phenotypes in [GEMMA](https://github.com/genetics-statistics/GEMMA) from context-specific unionbed files. <br/>
Additionally it uses [qqman_GEMMA.R](https://github.com/Dario-Galanti/EWAS/new/main/gemmaEWAS/qqman_GEMMA.R) to draw manhattan and qqplots and [enrichment_plot.R](https://github.com/Dario-Galanti/EWAS/new/main/gemmaEWAS/enrichment_plot.R) to plot enrichment of a priori candidate genes (if present).

[EWAS_manhattan_contColours.R](https://github.com/Dario-Galanti/EWAS/new/main/gemmaEWAS/EWAS_manhattan_contColours.R)<br/>
Draw manhattan and qqplots for all sequence contexts combined and for a single phenotype, producing graphs similar to these:
![image](https://user-images.githubusercontent.com/58292612/223488263-18394163-3e9a-4ff7-95db-eb877895257a.png)


[Features_enrichment_gemmaEWA.sh](https://github.com/Dario-Galanti/EWAS/new/main/gemmaEWAS/Features_enrichment_gemmaEWA.sh)<br/>
Run an enrichment analysis of variants in different genomic features (CDS, TEs, promoters ...) to test whether high -log(P) associations are enriched in a specific genomic feature.

[Features_enrichment_plot_multi.R](https://github.com/Dario-Galanti/EWAS/new/main/gemmaEWAS/Features_enrichment_plot_multi.R)<br/>
Plot results from [Features_enrichment_gemmaEWA.sh](https://github.com/Dario-Galanti/EWAS/new/main/gemmaEWAS/Features_enrichment_gemmaEWA.sh) for multiple phenotypes.
![image](https://user-images.githubusercontent.com/58292612/223489111-93669311-a691-47bf-b37d-b75af8a3cb1d.png)



