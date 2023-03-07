#!/bin/bash
## Author: Dario Galanti (January 2023)
## Aim: Calculate enrichment of associations in different genomic features for increasing -log(P) thresholds.
## Aim: Similar to what Atwell et al. 2010 does with a-priori candidates, but using whole genomic features instead
## Input: EWAS output directory of a specific phenotype with structure ${phenotype}/${context}/*.assoc.txt
## Output: txt file with headers -logP, feature, Cs_enrichment
## Run: bash Features_enrichment_gemmaEWA.sh /scr/episan/RP07/EWAS/gemmaEWAS/Exoreads/MpersicaeRes_noout2
## Run: sbatch --partition test --cpus-per-task 2 --mem 10G --time 8:00:00 --wrap "bash Features_enrichment_EWAS.sh /scr/episan/RP07/EWAS/EWAS_output"

# DEFINE STUFF
outdir=$1			## NB: Give the full path or the symlinks won't work
downDir=${outdir}/$(basename $outdir)_4download
ref_index=/scr/episan/RP07/Thlaspi_genome/Ta_v5_genome/final.fasta.fai
plot_script=/scr/episan/RP07/EWAS/Features_enrichment_plot.R
Rscript=~/conda/R4/bin/Rscript		#NB: Use R4 or a part of the script will crush
combined_results=${outdir}/All.bedGraph.txt
fout=${outdir}/Cs_enrichment_features_EWAS.txt


## Combine all contexts into single file and add -log(p)
echo -e ID beta stats pvalue FDR logP context | tr " " "\t" > $combined_results
awk 'FNR>1{OFS="\t";match(FILENAME, /C[pH][HG]/);cont=substr(FILENAME, RSTART, RLENGTH);logp=-log($12)/log(10);print $1":"$5"-"$6,$8,$9,$12,"*",logp,cont}' ${outdir}/C*/*.assoc.txt | sort -k1,1 -k2,2n >> $combined_results


## DEFINE FEATURE BED FILES
# NB: Feature filenames must contain feature name just before file extention
featureDir=/scr/episan/RP07/region_meth/Features_meth_v5/feature_beds
genes=${featureDir}/Ta_v5_genes.bed
CDS=${featureDir}/Ta_v5_CDS.bed
introns=${featureDir}/Ta_v5_introns.bed
TEs=${featureDir}/Ta_v5_TEs.bed
promoters=${featureDir}/Ta_v5_promoters.bed
intergenic=${featureDir}/Ta_v5_intergenic.bed

# MAKE FILE ARRAYS
reg_arr=( global $genes $CDS $introns $TEs $promoters $intergenic )


### CALCULATE ENRICHMENT IN DIFFERENT FEATURES FOR DIFFERENT -LOGP

# Print headers to fout
echo -e logP"\t"Feature"\t"Enrichment'\t'Observed_freq'\t'Expected_freq > $fout


## 1) ITERATE  THROUGH GENOMIC FEATURES (REGIONS)
 # First calculate for whole genome (no intersection)
for reg in ${reg_arr[@]};
do
 if [ $reg == "global" ];
 then
  feature=global
 else
  feature=$(basename $reg .bed | cut -d"_" -f3)
 fi

 ## 2) ITERATE THROUGH -LOGP (we multiply and then divide by 10 for simplicity)
 enrich_bed=${outdir}/Enrichment_${feature}.bed
 enrichment=${outdir}/Enrichment_${feature}.txt
 echo -e -logP'\t'Cand_sig_vrts'\t'Sig_vrts'\t'Observed_freq'\t'Expected_freq'\t'Enrichment'\t'Feature > $enrichment
 
 if [ $reg == "global" ];
 then
  awk 'NR>1{OFS="\t";split($1,a,":");split(a[2],b,"-");print a[1],b[1],b[2],$6,"1"}' $combined_results > $enrich_bed
 else
  awk 'NR>1{OFS="\t";split($1,a,":");split(a[2],b,"-");print a[1],b[1],b[2],$6}' $combined_results | bedtools intersect -a stdin -b $reg -wa -c > $enrich_bed
 fi
 exp_freq=$(awk '{if($5>0){cand++}} END {print cand/NR}' $enrich_bed)
 maxi=$(tail -n+2  $combined_results | awk 'BEFORE{m=0} {if($6>m){m=$6}} END {print int(m*10)}')
 echo -e -logP'\t'Cand_sig_vrts'\t'Sig_vrts'\t'Observed_freq'\t'Expected_freq'\t'Enrichment'\t'Feature > $enrichment
 
 for ((p=0;p<=maxi;p+=5));										# Iterate through -log(p) by 0.5 increasing steps and calculate enrichment
 do
  # For high -log(p) there can be 0 significant Cs, leading to division by 0 error with the old command
  awk -v logp=$p -v exp_fq=$exp_freq -v reg=$feature 'BEGIN{logp=(logp/10);cand=0;sig=0}{if($4>=logp){sig++;if($5>0){cand++}}} END{OFS="\t";if(sig>0){obs_fq=(cand/sig)}else{obs_fq=0};print logp,cand,sig,obs_fq,exp_fq,(obs_fq/exp_fq),reg}' $enrich_bed >> $enrichment
 done
 rm $enrich_bed
 
 ## COMBINE RESULTS FROM DIFFERENT FEATURES
 tail -n+2 $enrichment | awk 'OFS="\t"{print $1,$7,$6,$4,$5}' >> $fout
 #rm $enrichment
done

