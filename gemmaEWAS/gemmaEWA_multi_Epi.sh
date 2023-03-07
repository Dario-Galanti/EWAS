#!/bin/bash
## Author: Dario Galanti (Jan 2023)
## Aim: Run EWAS for multiple phenotypes with GEMMA with BIMBAM input files, correcting for population structure with a distance matrix.
## Aim: Enrichment of a-priori candidate variants can also be performed as in Atwell et al. 2010.
## Input: 1) Unionbed file --> Will be transformed to BIMBAM "Mean genoype file". It can contain 5% NAs. Filename must include sequence context (CpG, CHG or CHH)
## Input: 2) Phenotype file with headers "id" "empty" "pheno1" "pheno2" ...
## Input: 3) kinship matrix
## Input: 4) OPTIONAL: Flanked (20kb) a priori candidate genes in bed format for the enrichment analysis
## Documentation: https://www.xzlab.org/software/GEMMAmanual.pdf
## Run: bash gemmaEWA_multi_Epi.sh ${cont}_unionbed.bed /path/to/myphenotypes.txt outdir
## Run: bash gemmaEWA_multi_Epi.sh /scr/episan/RP07/EWAS/gemmaEWAS/input_v5/CpG_aphidDMRpos_v5_5NAs_2MEC_5ED.bed /scr/episan/RP07/EWAS/gemmaEWAS/pheno_v5/AphidRes.txt Exoreads

## Output: ${outdir}/${phenotype}/${cont}/* and ${outdir}/${phenotype}/${outdir}_4download/${cont}/*

##NB: GEMMA runs with two kinds of input files:
## - PLINK -> Easy to obtain and more symplistic
## - BIMBAM -> It stores dosage information or posterior mean genotypes (Is a value between 0 to 2 that can be interpreted as the
## minor allele dosage: 0 is homozygous for the major allele, 1 is a heterozygote, and 2 is a homozygote for the minor allele.)
## GEMMA is a GWAS software, but we use BIMBAM files and the “-notsnp” option, to feed it methylation values instead

## PROCESS AND RESOURCES:
## GEMMA cannot be run with multiple threads, so even though they have a model that runs multiple phenotypes, it is very slow
## So we run the single phenotype lmm and we parallelise submitting each phenotype as a different job to the queue

## DEFINE TOOLS
Rscript=~/conda/R4/bin/Rscript
gemma=~/conda/gemma/bin/gemma
bedtools=/usr/local/bin/bedtools
work=/scr/episan/RP07/EWAS/gemmaEWAS
qqman=${work}/qqman_GEMMA.R
enrich_script=/scr/episan/RP07/EWAS/gemmaEWAS/enrichment_plot.R
cd ${work}					# This allows to use relative paths

#flanked_candidates=/scr/episan/RP07/GWAS/candidates/flanked_DefResponse_TaAnno_v5_genes.bed		# Defence response (from Ta annotation with ancestor terms)
flanked_candidates=/scr/episan/RP07/GWAS/candidates/DefResp_AtOrtho_genes_BSMT1_20kbFlk.bed		# Defence response (from Ta anno with ancestor terms and BSMT1)

## DEFINE INPUT
union=$1
pheno_file=$2													# /beegfs/work/bbmdg01/GWAS/pheno_v5/SupFamTIPs_withUS_GWA_log_1out.txt
geno_BIM=${work}/BIM_geno/$(basename $1 .bed).bimbam			# Will be created and stored in the "BIM_geno" dir
pheno_BIM=${work}/BIM_pheno/$(basename $2 .txt).bimbam			# Will be created and stored in the "BIM_pheno" dir
cont=$(basename $union | grep -o 'C[pH][HG]')					# Extract context

## DEFINE K
kinship=${work}/input_v5/Ta_v3_vrts_1MAF_imp_8prun_METspls_withref.kinship
k_nohead=${work}/input_v5/Ta_v3_vrts_1MAF_imp_8prun_METspls_withref_nohead.kinship
## DEFINE OUTPUT
outdir=${work}/$3
mkdir -p $outdir BIM_geno BIM_pheno
report=${outdir}/gemmaEWA_report.txt

## 0) FIRST CHECK IF UNIONBED AND KINSHIP HAVE THE SAME SAMPLE NAMES
union_spls=$(head -1 $union | cut -f4- | tr -d '\r')
k_spls=$(head -1 $kinship | tr -d '\r')
if [ "$union_spls" != "$k_spls" ] ; then echo "ERROR: Sample names in unionbed and K matrix are not the same"; exit 1 ;fi
spls_arr=($union_spls)

### 1) CHECK AND FORMAT INPUT
pheno_arr=($(head -1 $pheno_file | cut -f3- | tr -d '\r'))
num_pheno=${#pheno_arr[@]}

## 1a) FORMAT PHENO (Order in the same way as geno and kinship files and add NAs for spls with missing phenotypes)
if [ ! -f "$pheno_BIM" ]; then
 for spl in "${spls_arr[@]}";
 do
  # grep values from phenotype file and replace empty fields with NAs (unix newline is \n, while DOS/windows is \r\n)
  val=$(grep -P $spl"\t" $pheno_file | cut -f3- | tr -d '\r' | awk -F"\t" '{OFS=" ";for(i=1;i<=NF;i++){if($i==""){$i="NA"}};print}')
  if [[ ${#val} == 0 ]];then val=$(seq $num_pheno | sed "c NA" | tr "\n" " ");fi
  echo $val >> $pheno_BIM
 done
 echo "Pheno BIM file created"
else
 echo "Pheno BIM file already present"
fi

## 1b) FORMAT KINSHIP
if [ ! -f "$k_nohead" ]; then tail -n+2 $kinship > $k_nohead;fi

## 1c) FORMAT GENO
## NB: 2nd and 3rd column are usually allele, but they are ignored by GEMMA, so we will store start and end, without a need for the SNP annotation file
if [ ! -f "$geno_BIM" ]; then
 awk 'NR>1{OFS="\t";mid=int((($2+$3)/2)+0.6);$1=$1":"mid; print}' $union > $geno_BIM
 echo "Epigeno BIM file created"
else
 echo "Epigeno BIM file already present"
fi


### PRINT SOME INFO
echo epigenotype file: $union > $report
echo phenotype file: $pheno_file >> $report
echo phenotypes: ${pheno_arr[@]} >> $report
echo kinship: $kinship >> $report


### LOOP THROUGH PHENOTYPES AND SUBMIT ONE JOB FOR EACH
mkdir -p ${work}/work/${3}
mkdir -p ${work}/logs/${3}
n=0
for phe in ${pheno_arr[@]};
do
	n=$((n+1))
	jobName=${work}/work/${3}/gemma.${cont}_${phe}.sh
	(
	echo "#!/bin/bash"
	echo "#SBATCH --cpus-per-task 2"	#Nodes and cores
	echo "#SBATCH --time 10:00:00"
	echo "#SBATCH --mem 8G"
	echo "#SBATCH --job-name gemma.${cont}_${phe}"
	echo "#SBATCH --partition test"
	echo "#SBATCH --output ${work}/logs/${3}/gemma.${cont}_${phe}.%A.out"
	echo "#SBATCH --error ${work}/logs/${3}/gemma.${cont}_${phe}.%A.err"
	echo ""
	
	## Define and prepare input and output
	echo "ph_dir=${outdir}/${phe}/${cont}"
	echo "results=\${ph_dir}/${phe}.assoc.txt"
	echo "mkdir -p \$ph_dir"
	echo "cd \$ph_dir"
	echo ""

	## 2) RUN GWAS
	echo "$gemma -g $geno_BIM -p $pheno_BIM -k $k_nohead -n $n -notsnp -lmm 1 -o $phe"
	echo "awk '{OFS=\"\t\";if(NR<2){print}else{split(\$2,p,\":\");\$1=p[1];\$3=p[2];print}}' output/*.assoc.txt > \$results"
	echo "rm -r output"
	
	## 3) VISUALIZATION (manhattan and qqplot)
	echo "$Rscript --vanilla $qqman \$results"
	echo ""
	
	## 4) ENRICHMENT OF A-PRIORI CANDIDATES
	echo "enrich_bed=Enrichment_${phe}.bed"
	echo "enrichment=Enrichment_${phe}.txt"
	echo "awk 'NR>1{OFS=\"\t\";logp=-log(\$12)/log(10);print \$1,\$5,\$6,logp}' \$results | $bedtools intersect -a stdin -b $flanked_candidates -wa -c > \$enrich_bed"
	echo "maxi=\$(awk 'BEFORE{m=0} {if(\$4>m){m=\$4}} END {print int(m*10)}' \$enrich_bed)"
	echo "exp_freq=\$(awk '{if(\$5>0){cand++}} END {print cand/NR}' \$enrich_bed)"
	# When running few positions/DMRs is possible that nothing is close to candidate genes.
	echo "if [ \$exp_freq -gt 0 ]"
	echo "then"
	echo " echo -e -logP'\t'Cand_sig_vrts'\t'Sig_vrts'\t'Observed_freq'\t'Expected_freq'\t'Enrichment > \$enrichment"
	echo " for ((p=0;p<=maxi;p+=5));"										# Iterate through -log(p) by 0.5 increasing steps and calculate enrichment
	echo " do"
	echo "  awk -v logp=\$p -v exp_fq=\$exp_freq 'BEGIN{logp=(logp/10)}{if(\$4>=logp){sig++;if(\$5>0){cand++}}} END{OFS=\"\t\";obs_fq=(cand/sig);print logp,cand,sig,obs_fq,exp_fq,(obs_fq/exp_fq)}' \$enrich_bed >> \$enrichment"
	echo " done"
	echo " $Rscript --vanilla $enrich_script \$enrichment"
	echo "else"
	echo "echo No marker is an a-priori candidate. Expected frequency is 0. We skip enrichment analysis"
	echo "fi"
	echo "rm \$enrich_bed"
	echo ""
	
	## 5): FORMAT RESULTS FOR IGV (https://software.broadinstitute.org/software/igv/GWAS) AND PREPARE FOLDER FOR DOWNLOAD
	## I remove all SNPs with -log(p)<1 to make the results lighter and easier to download and I format them in the correct way to use them in IGV
	echo "newdir=${outdir}/${phe}/${phe}_4download/${cont}"
	echo "mkdir -p \$newdir"
	echo "IGVfout=\${newdir}/${phe}_IGVresults.gwas"
	echo "Cs=\$(wc -l \$results | cut -d' ' -f1)"
	echo "if [ \$Cs -gt 500000 ];"			# Only remove low p-values in large files
	echo "then"
	echo " awk 'OFS=\"\t\"{if(NR==1){print \"chr\",\"bp\",\"snp\",\"p\"}else{if(\$12<0.1){print \$1,\$3,\$2,\$12}}}' \$results > \$IGVfout"
	echo "else"
	echo " awk 'OFS=\"\t\"{if(NR==1){print \"chr\",\"bp\",\"snp\",\"p\"}else{print \$1,\$3,\$2,\$12}}' \$results > \$IGVfout"
	echo "fi"
	echo "cp ./*.png \${newdir}"
	echo "cp ./*.pdf \${newdir}"
	
	) > $jobName
		chmod +x $jobName
		echo "Bash file $jobName created"
		sbatch ${jobName}
done
exit

