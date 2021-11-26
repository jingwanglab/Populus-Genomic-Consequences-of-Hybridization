#!/bin/sh
#SBATCH -c 1 --mem 4G
#SBATCH -o /UserHome/user260/pipeline/population_structure/populus227_pipeline/populus227.population_structure.out
#SBATCH -e /UserHome/user260/pipeline/population_structure/populus227_pipeline/populus227.population_structure.err

vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
bgzip="/UserHome/user260/tools/bcftools-1.2/htslib-1.2.1/bgzip"
plink="/data/apps/plink/20181202/plink"
structure="/UserHome/user260/tools/fastStructure/fastStructure/structure.py"
choosek="/UserHome/user260/tools/fastStructure/fastStructure/chooseK.py"
smartpca="/data/apps/eigensoft/EIG-7.2.1/bin/smartpca"
twstats="/data/apps/eigensoft/EIG-7.2.1/bin/twstats"
twtable="/data/apps/eigensoft/EIG-7.2.1/POPGEN/twtable"
admixture="/data/apps/admixture/1.3.0/admixture"

InputDir="/w/user260/multiple_aspen/populus259_vcf"
vcf="$InputDir/populus259.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz"

step=$1  ###whether to work on the whole dataset, or those removing clones or remove hybrids
method=$2  ###the methods to perform, e.g. admixture, faststructure, pca, snphylo...

populus227="/w/user260/multiple_aspen/populus227_vcf/samples/populus227.samples.remove_clones.csv"

OutDir="/w/user260/multiple_aspen/populus227_vcf/"

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

plink_OutDir=$OutDir/plink_ld0_2

if [ ! -d "$plink_OutDir" ]; then
mkdir -p $plink_OutDir
fi


##########################################################
##########################################################
###the step1 is to 
if [ $step == "1" ]; then

cut -d ',' -f 4 $populus227 |sed '1d' > populus227.keep.txt
$vcftools --gzvcf $vcf --keep populus227.keep.txt --max-missing 0.8 --maf 0.0000001 --recode --recode-INFO-all --out $OutDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed
$bcftools view $OutDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf -O z -o $OutDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz && rm $OutDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf


####################################################################################################
#####This script was used to prune the SNPs and Indels in order to do population structure analysis
###To make population structure analysis, there is several criteria to make:
#1. MAF>5%
#2. Missing rate<10%
#3. SNPs in LD were filtered with a window size of 50 SNPs (advancing 5 SNPs at a time) and a r2 threshold of 0.2, --indep 50 5 2
#--indep 50 5 2 where the 2 stands for the vif threshold (VIF is 1/(1-R^2)) which means in this case r2 = 0,50.
##r2=0.2
$plink --vcf $OutDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz --maf 0.05 --geno 0.1 --indep-pairwise 50 10 0.2 --allow-extra-chr --allow-no-sex --make-bed --out $plink_OutDir/populus227.snp.maf5geno10.ld_pruned0_2
##remove the SNPs located in scaffolds for downstream analysis
grep -v "scaffold" $plink_OutDir/populus227.snp.maf5geno10.ld_pruned0_2.prune.in > $plink_OutDir/populus227.snp.maf5geno10.ld_pruned0_2.no_scaff.prune.in 
$plink --vcf $OutDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz --extract $plink_OutDir/populus227.snp.maf5geno10.ld_pruned0_2.no_scaff.prune.in --recode --allow-extra-chr --allow-no-sex --make-bed --out $plink_OutDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned0_2


elif [ $step == "2" ]; then
###performing the faststructure and pca analysis

if [ $method == "faststructure" ]; then
bed="$plink_OutDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned0_2"
#
faststructure_OutDir=$plink_OutDir/faststructure/
#
if [ ! -d "$faststructure_OutDir" ]; then
mkdir -p $faststructure_OutDir
fi
#
#k=$3  ###the assumed k value
#python $structure -K $k --input=$bed --output=$faststructure_OutDir/populus227.maf5geno10.ld_pruned0_2.faststurcture.k$k --full --seed=100

#######################################################################3
##choose the best k
python $choosek --input=$faststructure_OutDir/populus227.maf5geno10.ld_pruned0_2.faststurcture



elif [ $method == "pca" ]; then

pca_OutDir=$plink_OutDir/pca

if [ ! -d "$pca_OutDir" ]; then
mkdir -p $pca_OutDir
fi

#1. performing pca analysis
#$smartpca -p parfile_smartpca.populus227.txt > smartpca.populus227.log

#2. using the twstats to select the most significant PCs
$twstats -t $twtable -i $pca_OutDir/populus227.pca.eval -o $pca_OutDir/populus227.twstats.out

#3.transform the evec file to .csv to read by r
sed 's/[ ][ ]*/,/g' $pca_OutDir/populus227.pca.evec | cut -d ',' -f 2- > $pca_OutDir/populus227.pca.evec.csv


fi
fi


