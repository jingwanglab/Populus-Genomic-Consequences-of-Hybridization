#!/bin/sh
#SBATCH -c 1 --mem 15G
#SBATCH -o /UserHome/user260/pipeline/snp_effects/dfe_mutation_load/sift4g.mutation_load.out
#SBATCH -e /UserHome/user260/pipeline/snp_effects/dfe_mutation_load/sift4g.mutation_load.err


###The main aim of this script is:
##1.extract the CDS mutation sites for the outgroup species, divide the class of the polymorphic sites
##2.extract the species-specific sites, 
##3.use the est-sfs to predict the unfolded sfs and the derived or ancestral state of the sites


vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
est_sfs="/UserHome/user260/tools/est-sfs-release-2.03/est-sfs"
plink2="/UserHome/user260/tools/plink/plink2"


cds_snp="/w/user260/multiple_aspen/populus162_vcf/vcf/sift/output/populus162.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.all_chr.SIFTannotations.cds.snps"
cds_sift="/w/user260/multiple_aspen/populus162_vcf/vcf/sift/output/populus162.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.all_chr.SIFTannotations.cds.xls"
peu_vcf="/w/user260/multiple_aspen/outgroups/PeuXBY08.snp.vcf.gz"
pti_vcf="/w/user260/multiple_aspen/outgroups/trichocarpa1.snp.vcf.gz"

species_vcf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf"
sift_dir="/w/user260/multiple_aspen/populus162_vcf/vcf/sift/output"

OutDir="/w/user260/multiple_aspen/populus162_vcf/sift4g"
if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

OutDir_outgroup=$OutDir/outgroup_vcf
if [ ! -d "$OutDir_outgroup" ]; then
mkdir -p $OutDir_outgroup
fi

step=$1
species=$2
snp_type=$3
species_sift=$species_vcf/sift_${species}

if [ ! -d "$species_sift" ]; then
mkdir -p $species_sift
fi

species_samples="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/species_samples/${species}.txt"

species_est_sfs=$species_sift/est_sfs
if [ ! -d "$species_est_sfs" ]; then
mkdir -p $species_est_sfs
fi

species_derived_hom_het=$species_est_sfs/derived_hom_het
if [ ! -d "$species_derived_hom_het" ]; then
mkdir -p $species_derived_hom_het
fi


if [ $step == "1" ]; then

#extracting the cds polymorphic sites
$vcftools --gzvcf $peu_vcf --snps $cds_snp --recode --recode-INFO-all --out $OutDir_outgroup/peu.cds_polymorphic
$bcftools view $OutDir_outgroup/peu.cds_polymorphic.recode.vcf -O z -o $OutDir_outgroup/peu.cds_polymorphic.recode.vcf.gz && rm $OutDir_outgroup/peu.cds_polymorphic.recode.vcf
$bcftools view -g ^miss -O z -o $OutDir_outgroup/peu.cds_polymorphic.no_miss.vcf.gz $OutDir_outgroup/peu.cds_polymorphic.recode.vcf.gz && rm $OutDir_outgroup/peu.cds_polymorphic.recode.vcf.gz

$vcftools --gzvcf $OutDir_outgroup/peu.cds_polymorphic.no_miss.vcf.gz --remove-indels --recode --recode-INFO-all --out $OutDir_outgroup/peu.cds_polymorphic.biallelic
$bcftools view $OutDir_outgroup/peu.cds_polymorphic.biallelic.recode.vcf -O z -o $OutDir_outgroup/peu.cds_polymorphic.biallelic.recode.vcf.gz && rm $OutDir_outgroup/peu.cds_polymorphic.biallelic.recode.vcf
mv $OutDir_outgroup/peu.cds_polymorphic.biallelic.recode.vcf.gz $OutDir_outgroup/peu.cds_polymorphic.no_miss.vcf.gz
$bcftools index $OutDir_outgroup/peu.cds_polymorphic.no_miss.vcf.gz

$vcftools --gzvcf $pti_vcf --snps $cds_snp --recode --recode-INFO-all --out $OutDir_outgroup/pti.cds_polymorphic
$bcftools view $OutDir_outgroup/pti.cds_polymorphic.recode.vcf -O z -o $OutDir_outgroup/pti.cds_polymorphic.recode.vcf.gz && rm $OutDir_outgroup/pti.cds_polymorphic.recode.vcf
$bcftools view -g ^miss -O z -o $OutDir_outgroup/pti.cds_polymorphic.no_miss.vcf.gz $OutDir_outgroup/pti.cds_polymorphic.recode.vcf.gz && rm $OutDir_outgroup/pti.cds_polymorphic.recode.vcf.gz

$vcftools --gzvcf $OutDir_outgroup/pti.cds_polymorphic.no_miss.vcf.gz --remove-indels --recode --recode-INFO-all --out $OutDir_outgroup/pti.cds_polymorphic.biallelic
$bcftools view $OutDir_outgroup/pti.cds_polymorphic.biallelic.recode.vcf -O z -o $OutDir_outgroup/pti.cds_polymorphic.biallelic.recode.vcf.gz && rm $OutDir_outgroup/pti.cds_polymorphic.biallelic.recode.vcf
mv $OutDir_outgroup/pti.cds_polymorphic.biallelic.recode.vcf.gz $OutDir_outgroup/pti.cds_polymorphic.no_miss.vcf.gz
$bcftools index $OutDir_outgroup/pti.cds_polymorphic.no_miss.vcf.gz


elif [ $step == "2" ]; then
#divide the class of the cds polymorphic sites into :  loss-of-function(START-LOSS, STOP-GAIN, STOP LOSS), deleterious, tolerate, synonymous
grep "START\|STOP" $cds_sift > $sift_dir/populus162.cds.loss_of_function.xls
awk '$9=="NONSYNONYMOUS"' $cds_sift |awk '$17=="DELETERIOUS"' > $sift_dir/populus162.cds.deleterious.xls
awk '$9=="NONSYNONYMOUS"' $cds_sift |grep -v "DELETER" |awk '$13!="NA"'> $sift_dir/populus162.cds.tolerated.xls
awk '$9=="SYNONYMOUS"' $cds_sift |grep -v "DELETERIOUS" |awk '$13!="NA"' > $sift_dir/populus162.cds.synonymous.xls

##extracting snp id
cut -f 1,2 $sift_dir/populus162.cds.loss_of_function.xls |sed 's/\t/:/g'> $sift_dir/populus162.cds.loss_of_function.snps
cut -f 1,2 $sift_dir/populus162.cds.deleterious.xls |sed 's/\t/:/g'> $sift_dir/populus162.cds.deleterious.snps
cut -f 1,2 $sift_dir/populus162.cds.tolerated.xls |sed 's/\t/:/g'> $sift_dir/populus162.cds.tolerated.snps
cut -f 1,2 $sift_dir/populus162.cds.synonymous.xls |sed 's/\t/:/g'> $sift_dir/populus162.cds.synonymous.snps


elif [ $step == "3" ]; then
###extract the various types of snps for each species and merge the species file with the outgroup files

$vcftools --gzvcf $species_vcf/${species}.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz --snps $sift_dir/populus162.cds.${snp_type}.snps --recode --recode-INFO-all --out $species_sift/${species}.snp.${snp_type}
$bcftools view $species_sift/${species}.snp.${snp_type}.recode.vcf -O z -o $species_sift/${species}.snp.${snp_type}.recode.vcf.gz && rm $species_sift/${species}.snp.${snp_type}.recode.vcf

$bcftools annotate -x FORMAT -O z -o $species_sift/${species}.snp.${snp_type}.simple.vcf.gz $species_sift/${species}.snp.${snp_type}.recode.vcf.gz && rm $species_sift/${species}.snp.${snp_type}.recode.vcf.gz
$bcftools index $species_sift/${species}.snp.${snp_type}.simple.vcf.gz 

#$bcftools isec -i- -i 'MAF[0]>=0' -i 'MAF[0]>=0' -c all -n=3 $species_sift/${species}.snp.${snp_type}.simple.vcf.gz $OutDir_outgroup/peu.cds_polymorphic.no_miss.vcf.gz $OutDir_outgroup/pti.cds_polymorphic.no_miss.vcf.gz -w 1 |grep -v "#" | cut -f 3 > $species_sift/${snp_type}.intersect.snp.txt
$bcftools isec -c all -n=3 $species_sift/${species}.snp.${snp_type}.simple.vcf.gz $OutDir_outgroup/peu.cds_polymorphic.no_miss.vcf.gz $OutDir_outgroup/pti.cds_polymorphic.no_miss.vcf.gz -w 1 |grep -v "#" | cut -f 3 > $species_sift/${snp_type}.intersect.snp.txt

$vcftools --gzvcf $species_sift/${species}.snp.${snp_type}.simple.vcf.gz --snps $species_sift/${snp_type}.intersect.snp.txt --recode --recode-INFO-all --out $species_sift/${species}.${snp_type}.intersect
$bcftools view $species_sift/${species}.${snp_type}.intersect.recode.vcf -O z -o $species_sift/${species}.${snp_type}.intersect.recode.vcf.gz && rm $species_sift/${species}.${snp_type}.intersect.recode.vcf
$bcftools index $species_sift/${species}.${snp_type}.intersect.recode.vcf.gz
#
$vcftools --gzvcf $OutDir_outgroup/peu.cds_polymorphic.no_miss.vcf.gz --snps $species_sift/${snp_type}.intersect.snp.txt --recode --recode-INFO-all --out $species_sift/peu.${snp_type}.intersect
$bcftools view $species_sift/peu.${snp_type}.intersect.recode.vcf -O z -o $species_sift/peu.${snp_type}.intersect.recode.vcf.gz && rm $species_sift/peu.${snp_type}.intersect.recode.vcf
$bcftools index $species_sift/peu.${snp_type}.intersect.recode.vcf.gz
#
$vcftools --gzvcf $OutDir_outgroup/pti.cds_polymorphic.no_miss.vcf.gz --snps $species_sift/${snp_type}.intersect.snp.txt --recode --recode-INFO-all --out $species_sift/pti.${snp_type}.intersect
$bcftools view $species_sift/pti.${snp_type}.intersect.recode.vcf -O z -o $species_sift/pti.${snp_type}.intersect.recode.vcf.gz && rm $species_sift/pti.${snp_type}.intersect.recode.vcf
$bcftools index $species_sift/pti.${snp_type}.intersect.recode.vcf.gz

$bcftools merge -m snps $species_sift/${species}.${snp_type}.intersect.recode.vcf.gz $species_sift/peu.${snp_type}.intersect.recode.vcf.gz $species_sift/pti.${snp_type}.intersect.recode.vcf.gz -O z -o $species_sift/${species}.${snp_type}.outgroup.vcf.gz


elif [ $step == "4" ]; then
##use script to extract the results to the input of est-sfs
species_est_sfs=$species_sift/est_sfs
if [ ! -d "$species_est_sfs" ]; then
mkdir -p $species_est_sfs
fi

#extract the specific species
$vcftools --gzvcf $species_sift/${species}.${snp_type}.outgroup.vcf.gz --keep /w/user260/multiple_aspen/populus162_vcf/samples/populus162.samples.txt --freq --out $species_est_sfs/${species}.${snp_type}.est_sfs
#extract outgroup1 
$vcftools --gzvcf $species_sift/${species}.${snp_type}.outgroup.vcf.gz --indv PeuXBY08 --freq --out $species_est_sfs/${species}.${snp_type}.outgroup1.est_sfs
#extract outgroup2
$vcftools --gzvcf $species_sift/${species}.${snp_type}.outgroup.vcf.gz --indv trichocarpa1 --freq --out $species_est_sfs/${species}.${snp_type}.outgroup2.est_sfs


elif [ $step == "5" ]; then
##use perl script to produce the input file of est-sfs
config_jc="/UserHome/user260/tools/est-sfs-release-2.03/config-JC.txt"
config_kimura="/UserHome/user260/tools/est-sfs-release-2.03/config-kimura.txt"
config_rate6="/UserHome/user260/tools/est-sfs-release-2.03/config-rate6.txt"
seed="/UserHome/user260/tools/est-sfs-release-2.03/seedfile.txt"

perl est_sfs_input.${species}.pl $species_est_sfs/${species}.${snp_type}.est_sfs.frq $species_est_sfs/${species}.${snp_type}.outgroup1.est_sfs.frq $species_est_sfs/${species}.${snp_type}.outgroup2.est_sfs.frq
#rate6
$est_sfs $config_rate6 $species_est_sfs/${species}.${snp_type}.est_sfs.input.txt $seed $species_est_sfs/${species}.${snp_type}.est_sfs.output.txt $species_est_sfs/${species}.${snp_type}.est_sfs.output.p_anc.txt
#jc
$est_sfs $config_jc $species_est_sfs/${species}.${snp_type}.est_sfs.input.txt $seed $species_est_sfs/${species}.${snp_type}.est_sfs.jc.output.txt $species_est_sfs/${species}.${snp_type}.est_sfs.jc.output.p_anc.txt
#kimura
$est_sfs $config_kimura $species_est_sfs/${species}.${snp_type}.est_sfs.input.txt $seed $species_est_sfs/${species}.${snp_type}.est_sfs.kimura.output.txt $species_est_sfs/${species}.${snp_type}.est_sfs.kimura.output.p_anc.txt


elif [ $step == "6" ]; then
###determine the major/minor allele and check the ancestral/derived state of snps
perl est_sfs_ancestral_derived.pl $species_est_sfs/${species}.${snp_type}.est_sfs.frq $species_est_sfs/${species}.${snp_type}.est_sfs.output.p_anc.txt 

elif [ $step == "7" ]; then
##for each indivdiual, extract the homozygotes and heterzygotes snps for various snp types
#the output format can be: individual, derived_heterozygotes, derived_homozygotes

awk '$3=="ALT"' $species_est_sfs/${species}.${snp_type}.est_sfs.ances_derive.txt |cut -f 1,2 |sed 's/\t/:/g' > $species_derived_hom_het/${species}.${snp_type}.derived.snps

$vcftools --gzvcf $species_sift/${species}.${snp_type}.outgroup.vcf.gz --snps $species_derived_hom_het/${species}.${snp_type}.derived.snps --recode --recode-INFO-all --out $species_derived_hom_het/${species}.${snp_type}.derived
$bcftools view $species_derived_hom_het/${species}.${snp_type}.derived.recode.vcf -O z -o $species_derived_hom_het/${species}.${snp_type}.derived.recode.vcf.gz && rm $species_derived_hom_het/${species}.${snp_type}.derived.recode.vcf

$plink2 --vcf $species_derived_hom_het/${species}.${snp_type}.derived.recode.vcf.gz --sample-counts --allow-extra-chr --out $species_derived_hom_het/${species}.${snp_type}.derived.hom_het

fi





