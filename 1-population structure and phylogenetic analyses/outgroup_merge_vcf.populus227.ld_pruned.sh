#!/bin/sh
#SBATCH -c 2 --mem 20G
#SBATCH -o /UserHome/user260/pipeline/phylo_tree/outgroup_merge.populus227.out
#SBATCH -e /UserHome/user260/pipeline/phylo_tree/outgroup_merge.populus227.err


##################################################################################
##################################################################################
#the script is used to merge outgroup species vcf into the 227 populus samples

##tools
vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
plink="/data/apps/plink/20181202/plink"
vcf2phylip="/UserHome/user260/tools/vcf2phylip-master/vcf2phylip.py"

##vcf
OutDir="/w/user260/multiple_aspen/outgroups"
tri1_vcf="/w/user261/data/outgroup_peu_ptri/peutri_unified/trichocarpa1.vcf.gz"
tri2_vcf="/w/user261/data/outgroup_peu_ptri/peutri_unified/trichocarpa2.vcf.gz"
peu1_vcf="/w/user261/data/outgroup_peu_ptri/peutri_unified/PeuXBY08.vcf.gz"
peu2_vcf="/w/user261/data/outgroup_peu_ptri/peutri_unified/PeuYL033.vcf.gz"
populus227_vcf="/w/user260/multiple_aspen/populus227_vcf/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz"


###include the SNPs after LD pruning for the 227 Populus samples
populus227_ld_bim="/w/user260/multiple_aspen/populus227_vcf/plink_ld0_2/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned0_2.bim"

InputDir="/w/user260/multiple_aspen/populus227_vcf/plink_ld0_2"

OutDir_vcf=$InputDir/vcf
if [ ! -d "$OutDir_vcf" ]; then
mkdir -p $OutDir_vcf
fi


OutDir_plink=$OutDir_vcf/plink_ibs_dist
if [ ! -d "$OutDir_plink" ]; then
mkdir -p $OutDir_plink
fi


step=$1


##########################################################
###the step1 is to add SNP ID to the four outgroup samples
if [ $step == "1" ]; then

$bcftools annotate --set-id +'%CHROM\:%POS' -O z -o $OutDir/trichocarpa1.snp.vcf.gz $tri1_vcf && rm $tri1_vcf
$bcftools annotate --set-id +'%CHROM\:%POS' -O z -o $OutDir/trichocarpa2.snp.vcf.gz $tri2_vcf && rm $tri2_vcf
$bcftools annotate --set-id +'%CHROM\:%POS' -O z -o $OutDir/PeuXBY08.snp.vcf.gz $peu1_vcf && rm $peu1_vcf
$bcftools annotate --set-id +'%CHROM\:%POS' -O z -o $OutDir/PeuYL033.snp.vcf.gz $peu2_vcf && rm $peu2_vcf


##########################################################
###step2: extract the SNPs-ID after filtering and also extract the vcf files for these SNPs
elif [ $step == "2" ]; then
###include the SNPs after LD pruning for the 227 Populus samples

cut -f 2 $populus227_ld_bim > $InputDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.snp 

$vcftools --gzvcf $populus227_vcf --snps $InputDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.snp --recode --recode-INFO-all --out $OutDir_vcf/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned
$bcftools view $OutDir_vcf/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.recode.vcf -O z -o $OutDir_vcf/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.recode.vcf.gz && rm $OutDir_vcf/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.recode.vcf

##trichocarpa1
$vcftools --gzvcf $OutDir/trichocarpa1.snp.vcf.gz --snps $InputDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.snp --recode --recode-INFO-all --out $OutDir_vcf/trichocarpa1.snp.ld_pruned
$bcftools view $OutDir_vcf/trichocarpa1.snp.ld_pruned.recode.vcf -O z -o $OutDir_vcf/trichocarpa1.snp.ld_pruned.recode.vcf.gz && rm $OutDir_vcf/trichocarpa1.snp.ld_pruned.recode.vcf

##trichocarpa2
$vcftools --gzvcf $OutDir/trichocarpa2.snp.vcf.gz --snps $InputDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.snp --recode --recode-INFO-all --out $OutDir_vcf/trichocarpa2.snp.ld_pruned
$bcftools view $OutDir_vcf/trichocarpa2.snp.ld_pruned.recode.vcf -O z -o $OutDir_vcf/trichocarpa2.snp.ld_pruned.recode.vcf.gz && rm $OutDir_vcf/trichocarpa2.snp.ld_pruned.recode.vcf

##PeuXBY08
$vcftools --gzvcf $OutDir/PeuXBY08.snp.vcf.gz --snps $InputDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.snp --recode --recode-INFO-all --out $OutDir_vcf/PeuXBY08.snp.ld_pruned
$bcftools view $OutDir_vcf/PeuXBY08.snp.ld_pruned.recode.vcf -O z -o $OutDir_vcf/PeuXBY08.snp.ld_pruned.recode.vcf.gz && rm $OutDir_vcf/PeuXBY08.snp.ld_pruned.recode.vcf

##PeuYL033
$vcftools --gzvcf $OutDir/PeuYL033.snp.vcf.gz --snps $InputDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.snp --recode --recode-INFO-all --out $OutDir_vcf/PeuYL033.snp.ld_pruned
$bcftools view $OutDir_vcf/PeuYL033.snp.ld_pruned.recode.vcf -O z -o $OutDir_vcf/PeuYL033.snp.ld_pruned.recode.vcf.gz && rm $OutDir_vcf/PeuYL033.snp.ld_pruned.recode.vcf


###step3: merge the populus227 samples and the four outgroups
##edit the vcf file
$bcftools annotate -x FORMAT -O z -o $OutDir_vcf/populus227.gt.vcf.gz $OutDir_vcf/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.maf5geno10.ld_pruned.recode.vcf.gz
$bcftools index $OutDir_vcf/populus227.gt.vcf.gz
$bcftools annotate -x FORMAT -O z -o $OutDir_vcf/trichocarpa1.gt.vcf.gz $OutDir_vcf/trichocarpa1.snp.ld_pruned.recode.vcf.gz
$bcftools index $OutDir_vcf/trichocarpa1.gt.vcf.gz
$bcftools annotate -x FORMAT -O z -o $OutDir_vcf/trichocarpa2.gt.vcf.gz $OutDir_vcf/trichocarpa2.snp.ld_pruned.recode.vcf.gz
$bcftools index $OutDir_vcf/trichocarpa2.gt.vcf.gz
$bcftools annotate -x FORMAT -O z -o $OutDir_vcf/PeuXBY08.gt.vcf.gz $OutDir_vcf/PeuXBY08.snp.ld_pruned.recode.vcf.gz
$bcftools index $OutDir_vcf/PeuXBY08.gt.vcf.gz
$bcftools annotate -x FORMAT -O z -o $OutDir_vcf/PeuYL033.gt.vcf.gz $OutDir_vcf/PeuYL033.snp.ld_pruned.recode.vcf.gz
$bcftools index $OutDir_vcf/PeuYL033.gt.vcf.gz

$bcftools view -g ^miss -O z -o $OutDir_vcf/trichocarpa1.gt.no_miss.vcf.gz $OutDir_vcf/trichocarpa1.gt.vcf.gz
$bcftools index $OutDir_vcf/trichocarpa1.gt.no_miss.vcf.gz

$bcftools view -g ^miss -O z -o $OutDir_vcf/trichocarpa2.gt.no_miss.vcf.gz $OutDir_vcf/trichocarpa2.gt.vcf.gz
$bcftools index $OutDir_vcf/trichocarpa2.gt.no_miss.vcf.gz

$bcftools view -g ^miss -O z -o $OutDir_vcf/PeuXBY08.gt.no_miss.vcf.gz $OutDir_vcf/PeuXBY08.gt.vcf.gz
$bcftools index $OutDir_vcf/PeuXBY08.gt.no_miss.vcf.gz

$bcftools view -g ^miss -O z -o $OutDir_vcf/PeuYL033.gt.no_miss.vcf.gz $OutDir_vcf/PeuYL033.gt.vcf.gz
$bcftools index $OutDir_vcf/PeuYL033.gt.no_miss.vcf.gz

$bcftools isec -n=5 $OutDir_vcf/populus227.gt.vcf.gz $OutDir_vcf/trichocarpa1.gt.no_miss.vcf.gz $OutDir_vcf/trichocarpa2.gt.no_miss.vcf.gz $OutDir_vcf/PeuXBY08.gt.no_miss.vcf.gz $OutDir_vcf/PeuYL033.gt.no_miss.vcf.gz -w 1 |grep -v "#" | cut -f 3 > $OutDir_vcf/intersect.snp.txt

$vcftools --gzvcf $OutDir_vcf/populus227.gt.vcf.gz --snps $OutDir_vcf/intersect.snp.txt --recode --recode-INFO-all --out $OutDir_vcf/populus227.gt.isec
$bcftools view $OutDir_vcf/populus227.gt.isec.recode.vcf -O z -o $OutDir_vcf/populus227.gt.isec.recode.vcf.gz && rm $OutDir_vcf/populus227.gt.isec.recode.vcf
$bcftools index $OutDir_vcf/populus227.gt.isec.recode.vcf.gz

$vcftools --gzvcf $OutDir_vcf/trichocarpa1.gt.no_miss.vcf.gz --snps $OutDir_vcf/intersect.snp.txt --recode --recode-INFO-all --out $OutDir_vcf/trichocarpa1.gt.no_miss.isec
$bcftools view $OutDir_vcf/trichocarpa1.gt.no_miss.isec.recode.vcf -O z -o $OutDir_vcf/trichocarpa1.gt.no_miss.isec.recode.vcf.gz && rm $OutDir_vcf/trichocarpa1.gt.no_miss.isec.recode.vcf
$bcftools index $OutDir_vcf/trichocarpa1.gt.no_miss.isec.recode.vcf.gz

$vcftools --gzvcf $OutDir_vcf/trichocarpa2.gt.no_miss.vcf.gz --snps $OutDir_vcf/intersect.snp.txt --recode --recode-INFO-all --out $OutDir_vcf/trichocarpa2.gt.no_miss.isec
$bcftools view $OutDir_vcf/trichocarpa2.gt.no_miss.isec.recode.vcf -O z -o $OutDir_vcf/trichocarpa2.gt.no_miss.isec.recode.vcf.gz && rm $OutDir_vcf/trichocarpa2.gt.no_miss.isec.recode.vcf
$bcftools index $OutDir_vcf/trichocarpa2.gt.no_miss.isec.recode.vcf.gz

$vcftools --gzvcf $OutDir_vcf/PeuXBY08.gt.no_miss.vcf.gz --snps $OutDir_vcf/intersect.snp.txt --recode --recode-INFO-all --out $OutDir_vcf/PeuXBY08.gt.no_miss.isec
$bcftools view $OutDir_vcf/PeuXBY08.gt.no_miss.isec.recode.vcf -O z -o $OutDir_vcf/PeuXBY08.gt.no_miss.isec.recode.vcf.gz && rm $OutDir_vcf/PeuXBY08.gt.no_miss.isec.recode.vcf
$bcftools index $OutDir_vcf/PeuXBY08.gt.no_miss.isec.recode.vcf.gz

$vcftools --gzvcf $OutDir_vcf/PeuYL033.gt.no_miss.vcf.gz --snps $OutDir_vcf/intersect.snp.txt --recode --recode-INFO-all --out $OutDir_vcf/PeuYL033.gt.no_miss.isec
$bcftools view $OutDir_vcf/PeuYL033.gt.no_miss.isec.recode.vcf -O z -o $OutDir_vcf/PeuYL033.gt.no_miss.isec.recode.vcf.gz && rm $OutDir_vcf/PeuYL033.gt.no_miss.isec.recode.vcf
$bcftools index $OutDir_vcf/PeuYL033.gt.no_miss.isec.recode.vcf.gz



$bcftools merge $OutDir_vcf/populus227.gt.isec.recode.vcf.gz $OutDir_vcf/trichocarpa1.gt.no_miss.isec.recode.vcf.gz $OutDir_vcf/trichocarpa2.gt.no_miss.isec.recode.vcf.gz $OutDir_vcf/PeuXBY08.gt.no_miss.isec.recode.vcf.gz $OutDir_vcf/PeuYL033.gt.no_miss.isec.recode.vcf.gz -O z -o $OutDir_vcf/populus227.outgroup.ld_pruned.vcf.gz


###basic filtering to only keep bi-allelic SNPs
$bcftools view -m2 -M2 -v snps -O z -o $OutDir_vcf/populus227.outgroup.ld_pruned.biallelic.vcf.gz $OutDir_vcf/populus227.outgroup.ld_pruned.vcf.gz



elif [ $step == "3" ]; then
####performing NJ tree building by 
$plink --vcf $OutDir_vcf/populus227.outgroup.ld_pruned.biallelic.vcf.gz --distance '1-ibs' --allow-extra-chr --allow-no-sex --out $OutDir_plink/populus227.ougroup.ld_pruned.snp

cut -f 1 $OutDir_plink/populus227.ougroup.ld_pruned.snp.mdist.id | awk 'BEGIN{a="#"}{printf("%s",a);for(i=1;i<=NF;i++){printf($i);printf("")}printf("%s","\n")}' > $OutDir_plink/populus227.ougroup.ld_pruned.snp.mdist.id.mega.txt

echo -e "#mega" > $OutDir_plink/populus227.outgroup.ld_pruned.snp.mdist.meg
echo -e "!Title: populus227.outgroup.ld_pruned.snp.mdist;" >>  $OutDir_plink/populus227.outgroup.ld_pruned.snp.mdist.meg
echo -e "" >> $OutDir_plink/populus227.outgroup.ld_pruned.snp.mdist.meg
cat $OutDir_plink/populus227.ougroup.ld_pruned.snp.mdist.id.mega.txt >> $OutDir_plink/populus227.outgroup.ld_pruned.snp.mdist.meg
echo -e "" >> $OutDir_plink/populus227.outgroup.ld_pruned.snp.mdist.meg
cat $OutDir_plink/populus227.ougroup.ld_pruned.snp.mdist >> $OutDir_plink/populus227.outgroup.ld_pruned.snp.mdist.meg

rm $OutDir_plink/populus227.ougroup.ld_pruned.snp.mdist.id.mega.txt


fi





