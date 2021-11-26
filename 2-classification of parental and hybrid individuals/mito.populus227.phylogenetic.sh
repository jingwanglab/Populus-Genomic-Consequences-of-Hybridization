#!/bin/sh
#SBATCH -c 1 --mem 8G
#SBATCH -o /UserHome/user260/pipeline/phylo_tree/mito_chro/mito.populus227.phylo.out
#SBATCH -e /UserHome/user260/pipeline/phylo_tree/mito_chro/mito.populus227.phylo.err

vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
bgzip="/UserHome/user260/tools/bcftools-1.2/htslib-1.2.1/bgzip"
export PERL5LIB=/data/apps/vcftools/vcftools-0.1.15/share/perl5:$PERL5LIB
vcf_annotate="/data/apps/vcftools/vcftools-0.1.15/bin/vcf-annotate"
plink="/data/apps/plink/20181202/plink"

InputDir="/w/user261/data/sangyupeng/Populus_Unified/mitochondrion/m_vcf"
m_vcf="$InputDir/m_merge.vcf.gz"
QD_filter="/UserHome/user260/pipeline/phylo_tree/mito_chro/QD_filter.txt"
populus227="/UserHome/user260/pipeline/phylo_tree/mito_chro/populus227/populus227.outgroup.txt"
populus227_reheader="/UserHome/user260/pipeline/phylo_tree/mito_chro/populus227/populus227.outgroup.reheader.txt"
populus227_old="/UserHome/user260/pipeline/phylo_tree/mito_chro/populus227/populus227.old.txt"
populus227_reheader_new="/UserHome/user260/pipeline/phylo_tree/mito_chro/populus227/populus227.txt"


step=$1

OutDir="/w/user260/multiple_aspen/populus227_vcf/mito_chro/mito"
if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

OutDir_phylo=$OutDir/phylo_nj
if [ ! -d "$OutDir_phylo" ]; then
mkdir -p $OutDir_phylo
fi


Out=${m_vcf##*/}
Out_snp_Suffix=${Out%.vcf.gz}.snp
snp_vcf=$Out_snp_Suffix.recode.vcf.gz

Out_biallelic_Suffix=${Out%.vcf.gz}.biallelic
biallelic_vcf=$Out_biallelic_Suffix.recode.vcf.gz

Out_DP_Suffix=${Out%.vcf.gz}.DP100
DP_vcf=$Out_DP_Suffix.recode.vcf.gz

Out_GQ_Suffix=${Out%.vcf.gz}.GQ30
GQ_vcf=$Out_GQ_Suffix.recode.vcf.gz

Out_populus227_Suffix=${Out%.vcf.gz}.populus227
populus227_vcf=$Out_populus227_Suffix.recode.vcf.gz

Out_populus227_reheader_Suffix=${Out%.vcf.gz}.populus227.reheader
populus227_reheader_vcf=$Out_populus227_reheader_Suffix.recode.vcf.gz

Out_missing_Suffix=${Out%.vcf.gz}.populus227.rm_missing
RM_vcf=$Out_missing_Suffix.recode.vcf.gz

Out_QD_Suffix=${Out%.vcf.gz}.populus227.QD
QD_vcf=$Out_QD_Suffix.recode.vcf.gz

Out_QD_filtering_Suffix=${Out%.vcf.gz}.populus227.QD.filter
QD_filtering_vcf=$Out_QD_filtering_Suffix.recode.vcf.gz

Out_non_ref_Suffix=${Out%.vcf.gz}.populus227.maf
Non_ref_vcf=$Out_non_ref_Suffix.recode.vcf.gz

Out_maf_Suffix=${Out%.vcf.gz}.populus227.maf

###########################################################################################################
###########################################################################################################
###step1: filtering on the called SNPs

if [ $step == "1" ]; then 
#remove indels
$vcftools --gzvcf $m_vcf --remove-indels --out $OutDir/$Out_snp_Suffix --recode --recode-INFO-all
$bgzip $OutDir/$Out_snp_Suffix.recode.vcf

#Bi-alleleic filtering
$vcftools --gzvcf $OutDir/$snp_vcf --max-alleles 2 --recode --recode-INFO-all --out $OutDir/$Out_biallelic_Suffix
$bgzip $OutDir/$Out_biallelic_Suffix.recode.vcf

#Depth filtering
$vcftools --gzvcf $OutDir/$biallelic_vcf --minDP 100 --recode --recode-INFO-all --out $OutDir/$Out_DP_Suffix
$bgzip $OutDir/$Out_DP_Suffix.recode.vcf

#Genotype quality filtering
$vcftools --gzvcf $OutDir/$DP_vcf --minGQ 30 --recode --recode-INFO-all --out $OutDir/$Out_GQ_Suffix
$bgzip $OutDir/$Out_GQ_Suffix.recode.vcf

##extracting populus227 samples
$vcftools --gzvcf $OutDir/$Out_GQ_Suffix.recode.vcf.gz --keep $populus227_old --recode --recode-INFO-all --out $OutDir/$Out_populus227_Suffix
$bgzip $OutDir/$Out_populus227_Suffix.recode.vcf

##reheader the populus227 samples
$bcftools reheader -s $populus227_reheader $OutDir/$Out_populus227_Suffix.recode.vcf.gz -o $OutDir/$Out_populus227_reheader_Suffix.recode.vcf.gz 

#Max-missing filtering
$vcftools --gzvcf $OutDir/$populus227_reheader_vcf --max-missing 0.8 --recode --recode-INFO-all --out $OutDir/$Out_missing_Suffix
$bgzip $OutDir/$Out_missing_Suffix.recode.vcf

#Quality by depth filtering
zcat $OutDir/$RM_vcf | $vcf_annotate -f $QD_filter > $OutDir/$Out_QD_Suffix.recode.vcf
$bgzip $OutDir/$Out_QD_Suffix.recode.vcf
$vcftools --gzvcf $OutDir/$QD_vcf --remove-filtered MinQD --recode --recode-INFO-all --out $OutDir/$Out_QD_filtering_Suffix
$bgzip $OutDir/$Out_QD_filtering_Suffix.recode.vcf
#
#remove maf <0.000001 non-polymorhic sites
##3,902 SNPs
$vcftools --gzvcf $OutDir/$QD_filtering_vcf --maf 0.000001 --recode --recode-INFO-all --out $OutDir/$Out_non_ref_Suffix
$bgzip $OutDir/$Out_non_ref_Suffix.recode.vcf

##only extract those populus 259 samples without outgroup samples
#2,915 SNPs
$vcftools --gzvcf $OutDir/$Out_non_ref_Suffix.recode.vcf.gz --keep $populus227_reheader_new --maf 0.000001 --recode --recode-INFO-all --out $OutDir/${Out_non_ref_Suffix}.no_outgroup
$bgzip $OutDir/${Out_non_ref_Suffix}.no_outgroup.recode.vcf

rm $OutDir/$snp_vcf $OutDir/$biallelic_vcf $OutDir/$DP_vcf $OutDir/$GQ_vcf $OutDir/$Out_populus227_Suffix.recode.vcf.gz $OutDir/$Out_populus227_reheader_Suffix.recode.vcf.gz $OutDir/$QD_vcf $OutDir/$RM_vcf $OutDir/$QD_filtering_vcf


###########################################################################################################
###########################################################################################################
###step2: reconstruct the phylogenetic tree of all samples with outgroup
elif [ $step == "2" ]; then 
##3895 Sites


$plink --vcf $OutDir/$Out_non_ref_Suffix.recode.vcf.gz --distance '1-ibs' --allow-extra-chr --allow-no-sex --out $OutDir_phylo/$Out_non_ref_Suffix

cut -f 1 $OutDir_phylo/${Out_non_ref_Suffix}.mdist.id | awk 'BEGIN{a="#"}{printf("%s",a);for(i=1;i<=NF;i++){printf($i);printf("")}printf("%s","\n")}' > $OutDir_phylo/${Out_non_ref_Suffix}.mdist.id.mega.txt


echo -e "#mega" > $OutDir_phylo/${Out_non_ref_Suffix}.mdist.meg
echo -e "!Title: populus227.outgroup.mito.snp.mdist;" >> $OutDir_phylo/${Out_non_ref_Suffix}.mdist.meg
echo -e "" >> $OutDir_phylo/${Out_non_ref_Suffix}.mdist.meg
cat $OutDir_phylo/${Out_non_ref_Suffix}.mdist.id.mega.txt >> $OutDir_phylo/${Out_non_ref_Suffix}.mdist.meg
echo -e "" >> $OutDir_phylo/${Out_non_ref_Suffix}.mdist.meg
cat $OutDir_phylo/${Out_non_ref_Suffix}.mdist >> $OutDir_phylo/${Out_non_ref_Suffix}.mdist.meg

rm $OutDir_phylo/${Out_non_ref_Suffix}.mdist.id.mega.txt 


elif [ $step == "3" ]; then
###only use SNPs with no missing for outgroup
##3387 SNPs


$vcftools --gzvcf $OutDir/$Out_non_ref_Suffix.recode.vcf.gz --indv PeuXBY08 --indv PeuYL033 --indv trichocarpa1 --indv trichocarpa2 --max-missing 1.0 --recode --recode-INFO-all --out $OutDir/outgroup
$bgzip $OutDir/${Out_non_ref_Suffix}.no_outgroup.recode.vcf
$bcftools annotate --set-id +'%CHROM\:%POS' -o $OutDir/outgroup.snp.vcf $OutDir/outgroup.recode.vcf
grep -v "#" $OutDir/outgroup.snp.vcf |cut -f 3 > $OutDir/outgroup.snp

$bcftools annotate --set-id +'%CHROM\:%POS' -O z -o $OutDir/$Out_non_ref_Suffix.recode2.vcf.gz $OutDir/$Out_non_ref_Suffix.recode.vcf.gz
$vcftools --gzvcf $OutDir/$Out_non_ref_Suffix.recode2.vcf.gz --snps $OutDir/outgroup.snp --recode --recode-INFO-all --out $OutDir/populus227.outgroup.snp
$bgzip $OutDir/populus227.outgroup.snp.recode.vcf

rm $OutDir/outgroup.recode.vcf $OutDir/outgroup.snp $OutDir/outgroup.snp.vcf $OutDir/$Out_non_ref_Suffix.recode2.vcf.gz


$plink --vcf $OutDir/populus227.outgroup.snp.recode.vcf.gz --distance '1-ibs' --allow-extra-chr --allow-no-sex --out $OutDir_phylo/populus227.outgroup.snp

cut -f 1 $OutDir_phylo/populus227.outgroup.snp.mdist.id | awk 'BEGIN{a="#"}{printf("%s",a);for(i=1;i<=NF;i++){printf($i);printf("")}printf("%s","\n")}' > $OutDir_phylo/populus227.outgroup.snp.mdist.id.mega.txt


echo -e "#mega" > $OutDir_phylo/populus227.outgroup.snp.mdist.meg
echo -e "!Title: populus227.outgroup.snp.mito.snp.mdist;" >> $OutDir_phylo/populus227.outgroup.snp.mdist.meg
echo -e "" >> $OutDir_phylo/populus227.outgroup.snp.mdist.meg
cat $OutDir_phylo/populus227.outgroup.snp.mdist.id.mega.txt >> $OutDir_phylo/populus227.outgroup.snp.mdist.meg
echo -e "" >> $OutDir_phylo/populus227.outgroup.snp.mdist.meg
cat $OutDir_phylo/populus227.outgroup.snp.mdist >> $OutDir_phylo/populus227.outgroup.snp.mdist.meg

rm $OutDir_phylo/populus227.outgroup.snp.mdist.id.mega.txt

elif [ $step == "4" ]; then
######performing further filtering of missing rate
#no missing for all samples
#2589 out of a possible 3387 Sites
$vcftools --gzvcf $OutDir/populus227.outgroup.snp.recode.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out $OutDir/populus227.outgroup.snp.no_miss
$bgzip $OutDir/populus227.outgroup.snp.no_miss.recode.vcf

$plink --vcf $OutDir/populus227.outgroup.snp.no_miss.recode.vcf.gz --distance '1-ibs' --allow-extra-chr --allow-no-sex --out $OutDir_phylo/populus227.outgroup.snp.no_miss

cut -f 1 $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.id | awk 'BEGIN{a="#"}{printf("%s",a);for(i=1;i<=NF;i++){printf($i);printf("")}printf("%s","\n")}' > $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.id.mega.txt


echo -e "#mega" > $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.meg
echo -e "!Title: populus227.outgroup.snp.mito.snp.no_miss.mdist;" >> $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.meg
echo -e "" >> $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.meg
cat $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.id.mega.txt >> $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.meg
echo -e "" >> $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.meg
cat $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist >> $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.meg

rm $OutDir_phylo/populus227.outgroup.snp.no_miss.mdist.id.mega.txt


elif [ $step == "5" ]; then
###without outgroup


$plink --vcf $OutDir/${Out_non_ref_Suffix}.no_outgroup.recode.vcf.gz --distance '1-ibs' --allow-extra-chr --allow-no-sex --out $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup

cut -f 1 $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.id | awk 'BEGIN{a="#"}{printf("%s",a);for(i=1;i<=NF;i++){printf($i);printf("")}printf("%s","\n")}' > $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.id.mega.txt


echo -e "#mega" > $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.meg
echo -e "!Title: populus227.mito.snp.mdist;" >> $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.meg
echo -e "" >> $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.meg
cat $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.id.mega.txt >> $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.meg
echo -e "" >> $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.meg
cat $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist >> $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.meg

rm $OutDir_phylo/${Out_non_ref_Suffix}.no_outgroup.mdist.id.mega.txt


fi



