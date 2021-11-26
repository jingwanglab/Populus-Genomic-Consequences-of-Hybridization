#!/bin/sh
#SBATCH -c 1 --mem 10G
#SBATCH -o /UserHome/user260/pipeline/snp_effects/dfe_mutation_load/dfe.0_4fold.out
#SBATCH -e /UserHome/user260/pipeline/snp_effects/dfe_mutation_load/dfe.0_4fold.err


###The main aim of this script is:
#1.extract the intersection of 0-fold and 4-fold sites and the filtered regions 
#2.extract the  sites of the above filtered 0-fold and 4-fold sites for the outgroup species
#3.use est-sfs to predict the derived/ancestral allele for all 0-fold and 4-fold sites with at least one outgroup species information

vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
est_sfs="/UserHome/user260/tools/est-sfs-release-2.03/est-sfs"
est_dfe="/UserHome/user260/tools/dfe_alpha/dfe-alpha-release-2.16/est_dfe"
est_alpha_omega="/UserHome/user260/tools/dfe_alpha/dfe-alpha-release-2.16/est_alpha_omega"
prop_muts_in_s_ranges="/UserHome/user260/tools/dfe_alpha/dfe-alpha-release-2.16/prop_muts_in_s_ranges"


zero_fold="/w/user260/multiple_aspen/gff/0_4_fold/Pte_degenerate_0fold.no_scaff.bed"
four_fold="/w/user260/multiple_aspen/gff/0_4_fold/Pte_degenerate_4fold.no_scaff.bed"
filter_bed="/w/user260/multiple_aspen/snpable/Potra02_genome.mask_mappability.no_scaffold.bed"
peu_vcf="/w/user260/multiple_aspen/outgroups/PeuXBY08.snp.vcf.gz"
pti_vcf="/w/user260/multiple_aspen/outgroups/trichocarpa1.snp.vcf.gz"

Inputvcf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf"


step=$1
species=$2

OutDir="/w/user260/multiple_aspen/populus162_vcf/dfe_alpha"
if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

OutDir_outgroup=$OutDir/outgroup_vcf
if [ ! -d "$OutDir_outgroup" ]; then
mkdir -p $OutDir_outgroup
fi


OutDir_species=$OutDir/$species
if [ ! -d "$OutDir_species" ]; then
mkdir -p $OutDir_species
fi

OutDir_est_sfs=$OutDir/$species/est_sfs
if [ ! -d "$OutDir_est_sfs" ]; then
mkdir -p $OutDir_est_sfs
fi

OutDir_dfe=$OutDir/$species/est_sfs/dfe_alpha
if [ ! -d "$OutDir_dfe" ]; then
mkdir -p $OutDir_dfe
fi

OutDir_bootstrap=$OutDir_est_sfs/bootstrap
if [ ! -d "$OutDir_bootstrap" ]; then
mkdir -p $OutDir_bootstrap
fi

OutDir_summary=$OutDir_est_sfs/summary
if [ ! -d "$OutDir_summary" ]; then
mkdir -p $OutDir_summary
fi



filter_folder="/w/user260/multiple_aspen/gff/0_4_fold/filtered"

config_jc="/UserHome/user260/tools/est-sfs-release-2.03/config-JC.txt"
config_kimura="/UserHome/user260/tools/est-sfs-release-2.03/config-kimura.txt"
config_rate6="/UserHome/user260/tools/est-sfs-release-2.03/config-rate6.txt"
seed="/UserHome/user260/tools/est-sfs-release-2.03/seedfile.txt"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/UserHome/user260/tools/dfe_alpha/

if [ $step == "1" ]; then

$bedtools intersect -a $filter_bed -b $zero_fold > $filter_folder/Pte_degenerate_0fold.filter.bed
$bedtools intersect -a $filter_bed -b $four_fold > $filter_folder/Pte_degenerate_4fold.filter.bed

elif [ $step == "2" ]; then
#extracting the 0_fold sites 
#$vcftools --gzvcf $peu_vcf --bed $filter_folder/Pte_degenerate_0fold.filter.bed --recode --recode-INFO-all --out $OutDir_outgroup/peu.0fold.filter
#$bcftools index $peu_vcf
#$bcftools view $peu_vcf -R $filter_folder/Pte_degenerate_0fold.filter.bed -O z -o $OutDir_outgroup/peu.0fold.filter.recode.vcf.gz
#$bcftools view $OutDir_outgroup/peu.0fold.filter.recode.vcf -O z -o $OutDir_outgroup/peu.0fold.filter.recode.vcf.gz && rm $OutDir_outgroup/peu.0fold.filter.recode.vcf

#$vcftools --gzvcf $pti_vcf --bed $filter_folder/Pte_degenerate_0fold.filter.bed --recode --recode-INFO-all --out $OutDir_outgroup/pti.0fold.filter
#$bcftools view $OutDir_outgroup/pti.0fold.filter.recode.vcf -O z -o $OutDir_outgroup/pti.0fold.filter.recode.vcf.gz && rm $OutDir_outgroup/pti.0fold.filter.recode.vcf

$bcftools index $pti_vcf
$bcftools view $pti_vcf -R $filter_folder/Pte_degenerate_0fold.filter.bed -O z -o $OutDir_outgroup/pti.0fold.filter.2.recode.vcf.gz

elif [ $step == "2.1" ]; then
##extracting the 4_fold sites
#$vcftools --gzvcf $peu_vcf --bed $filter_folder/Pte_degenerate_4fold.filter.bed --recode --recode-INFO-all --out $OutDir_outgroup/peu.4fold.filter
#$bcftools view $OutDir_outgroup/peu.4fold.filter.recode.vcf -O z -o $OutDir_outgroup/peu.4fold.filter.recode.vcf.gz && rm $OutDir_outgroup/peu.4fold.filter.recode.vcf
#$bcftools index $pti_vcf
#$bcftools view $peu_vcf -R $filter_folder/Pte_degenerate_4fold.filter.bed -O z -o $OutDir_outgroup/peu.4fold.filter.recode.vcf.gz

#$vcftools --gzvcf $pti_vcf --bed $filter_folder/Pte_degenerate_4fold.filter.bed --recode --recode-INFO-all --out $OutDir_outgroup/pti.4fold.filter
#$bcftools view $OutDir_outgroup/pti.4fold.filter.recode.vcf -O z -o $OutDir_outgroup/pti.4fold.filter.recode.vcf.gz && rm $OutDir_outgroup/pti.4fold.filter.recode.vcf
$bcftools view $pti_vcf -R $filter_folder/Pte_degenerate_4fold.filter.bed -O z -o $OutDir_outgroup/pti.4fold.filter.2.recode.vcf.gz


elif [ $step == "3.0" ]; then
###index vcf
$bcftools view $Inputvcf/${species}.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf -O z -o $Inputvcf/${species}.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz && rm $Inputvcf/${species}.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf
$bcftools index $Inputvcf/${species}.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz

elif [ $step == "3.1" ]; then
###extract the 4-fold polymorphic sites for each species

$bcftools view $Inputvcf/${species}.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz -R $filter_folder/Pte_degenerate_4fold.filter.bed -O z -o $OutDir_species/${species}.snp.4fold.vcf.gz

elif [ $step == "3.2" ]; then
###extract the 0-fold polymorphic sites for each species

$bcftools view $Inputvcf/${species}.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz -R $filter_folder/Pte_degenerate_0fold.filter.bed -O z -o $OutDir_species/${species}.snp.0fold.vcf.gz

elif [ $step == "4" ]; then
###estimate frequency of alternative allele
$vcftools --gzvcf $OutDir_species/${species}.snp.4fold.vcf.gz --freq --out $OutDir_est_sfs/${species}.snp.4fold
$vcftools --gzvcf $OutDir_species/${species}.snp.0fold.vcf.gz --freq --out $OutDir_est_sfs/${species}.snp.0fold

elif [ $step == "4.1" ]; then
###the frequency of outgroup species

#peu
##remove missing and indel snps
$vcftools --gzvcf $OutDir_outgroup/peu.4fold.filter.recode.vcf.gz --remove-indels --recode --recode-INFO-all --out $OutDir_outgroup/peu.4fold.filter.no_indel
$bcftools view -g ^miss -O z -o $OutDir_outgroup/peu.4fold.filter.no_miss.vcf.gz $OutDir_outgroup/peu.4fold.filter.no_indel.recode.vcf && rm $OutDir_outgroup/peu.4fold.filter.no_indel.recode.vcf

$vcftools --gzvcf $OutDir_outgroup/peu.0fold.filter.recode.vcf.gz --remove-indels --recode --recode-INFO-all --out $OutDir_outgroup/peu.0fold.filter.no_indel
$bcftools view -g ^miss -O z -o $OutDir_outgroup/peu.0fold.filter.no_miss.vcf.gz $OutDir_outgroup/peu.0fold.filter.no_indel.recode.vcf && rm $OutDir_outgroup/peu.0fold.filter.no_indel.recode.vcf


elif [ $step == "4.2" ]; then
#pti
$vcftools --gzvcf $OutDir_outgroup/pti.4fold.filter.recode.vcf.gz --remove-indels --recode --recode-INFO-all --out $OutDir_outgroup/pti.4fold.filter.no_indel
$bcftools view -g ^miss -O z -o $OutDir_outgroup/pti.4fold.filter.no_miss.vcf.gz $OutDir_outgroup/pti.4fold.filter.no_indel.recode.vcf && rm $OutDir_outgroup/pti.4fold.filter.no_indel.recode.vcf

$vcftools --gzvcf $OutDir_outgroup/pti.0fold.filter.recode.vcf.gz --remove-indels --recode --recode-INFO-all --out $OutDir_outgroup/pti.0fold.filter.no_indel
$bcftools view -g ^miss -O z -o $OutDir_outgroup/pti.0fold.filter.no_miss.vcf.gz $OutDir_outgroup/pti.0fold.filter.no_indel.recode.vcf && rm $OutDir_outgroup/pti.0fold.filter.no_indel.recode.vcf


elif [ $step == "4.3" ]; then
#merge two outgroup vcfs
#4fold
#$bcftools index $OutDir_outgroup/peu.4fold.filter.no_miss.vcf.gz
#$bcftools index  $OutDir_outgroup/pti.4fold.filter.no_miss.vcf.gz

#$bcftools isec -c all -n=2 $OutDir_outgroup/peu.4fold.filter.no_miss.vcf.gz $OutDir_outgroup/pti.4fold.filter.no_miss.vcf.gz -w 1 |grep -v "#" |cut -f 3 > $OutDir_outgroup/peu.pti.4fold.intersect.snp.txt

$vcftools --gzvcf $OutDir_outgroup/peu.4fold.filter.no_miss.vcf.gz --snps $OutDir_outgroup/peu.pti.4fold.intersect.snp.txt --recode --recode-INFO-all --out $OutDir_outgroup/peu.4fold.filter.no_miss.intersect
$bcftools view $OutDir_outgroup/peu.4fold.filter.no_miss.intersect.recode.vcf -O z -o $OutDir_outgroup/peu.4fold.filter.no_miss.intersect.vcf.gz && rm $OutDir_outgroup/peu.4fold.filter.no_miss.intersect.recode.vcf
$bcftools index $OutDir_outgroup/peu.4fold.filter.no_miss.intersect.vcf.gz

$vcftools --gzvcf $OutDir_outgroup/pti.4fold.filter.no_miss.vcf.gz --snps $OutDir_outgroup/peu.pti.4fold.intersect.snp.txt --recode --recode-INFO-all --out $OutDir_outgroup/pti.4fold.filter.no_miss.intersect
$bcftools view $OutDir_outgroup/pti.4fold.filter.no_miss.intersect.recode.vcf -O z -o $OutDir_outgroup/pti.4fold.filter.no_miss.intersect.vcf.gz && rm $OutDir_outgroup/pti.4fold.filter.no_miss.intersect.recode.vcf
$bcftools index $OutDir_outgroup/pti.4fold.filter.no_miss.intersect.vcf.gz

#merge
$bcftools merge -m snps $OutDir_outgroup/peu.4fold.filter.no_miss.intersect.vcf.gz $OutDir_outgroup/pti.4fold.filter.no_miss.intersect.vcf.gz -O z -o $OutDir_outgroup/peu_pti.4fold.filter.no_miss.intersect.vcf.gz && rm $OutDir_outgroup/peu.4fold.filter.no_miss.intersect.vcf.gz $OutDir_outgroup/pti.4fold.filter.no_miss.intersect.vcf.gz 


elif [ $step == "4.4" ]; then
#0fold
#$bcftools index $OutDir_outgroup/peu.0fold.filter.no_miss.vcf.gz
#$bcftools index  $OutDir_outgroup/pti.0fold.filter.no_miss.vcf.gz

#$bcftools isec -c all -n=2 $OutDir_outgroup/peu.0fold.filter.no_miss.vcf.gz $OutDir_outgroup/pti.0fold.filter.no_miss.vcf.gz -w 1 |grep -v "#" |cut -f 3 > $OutDir_outgroup/peu.pti.0fold.intersect.snp.txt

$vcftools --gzvcf $OutDir_outgroup/peu.0fold.filter.no_miss.vcf.gz --snps $OutDir_outgroup/peu.pti.0fold.intersect.snp.txt --recode --recode-INFO-all --out $OutDir_outgroup/peu.0fold.filter.no_miss.intersect
$bcftools view $OutDir_outgroup/peu.0fold.filter.no_miss.intersect.recode.vcf -O z -o $OutDir_outgroup/peu.0fold.filter.no_miss.intersect.vcf.gz && rm $OutDir_outgroup/peu.0fold.filter.no_miss.intersect.recode.vcf
$bcftools index $OutDir_outgroup/peu.0fold.filter.no_miss.intersect.vcf.gz

$vcftools --gzvcf $OutDir_outgroup/pti.0fold.filter.no_miss.vcf.gz --snps $OutDir_outgroup/peu.pti.0fold.intersect.snp.txt --recode --recode-INFO-all --out $OutDir_outgroup/pti.0fold.filter.no_miss.intersect
$bcftools view $OutDir_outgroup/pti.0fold.filter.no_miss.intersect.recode.vcf -O z -o $OutDir_outgroup/pti.0fold.filter.no_miss.intersect.vcf.gz && rm $OutDir_outgroup/pti.0fold.filter.no_miss.intersect.recode.vcf
$bcftools index $OutDir_outgroup/pti.0fold.filter.no_miss.intersect.vcf.gz

#merge
$bcftools merge -m snps $OutDir_outgroup/peu.0fold.filter.no_miss.intersect.vcf.gz $OutDir_outgroup/pti.0fold.filter.no_miss.intersect.vcf.gz -O z -o $OutDir_outgroup/peu_pti.0fold.filter.no_miss.intersect.vcf.gz && rm $OutDir_outgroup/peu.0fold.filter.no_miss.intersect.vcf.gz $OutDir_outgroup/pti.0fold.filter.no_miss.intersect.vcf.gz 


elif [ $step == "4.5" ]; then
##estimate the frequency of the ref and alt alleles for the two outgroups

#extract outgroup1
$vcftools --gzvcf $OutDir_outgroup/peu_pti.4fold.filter.no_miss.intersect.vcf.gz --indv PeuXBY08 --freq --out $OutDir_outgroup/est_sfs/peu.4fold.est_sfs
$vcftools --gzvcf $OutDir_outgroup/peu_pti.0fold.filter.no_miss.intersect.vcf.gz --indv PeuXBY08 --freq --out $OutDir_outgroup/est_sfs/peu.0fold.est_sfs

#extract outgroup2
$vcftools --gzvcf $OutDir_outgroup/peu_pti.4fold.filter.no_miss.intersect.vcf.gz --indv trichocarpa1 --freq --out $OutDir_outgroup/est_sfs/pti.4fold.est_sfs
$vcftools --gzvcf $OutDir_outgroup/peu_pti.0fold.filter.no_miss.intersect.vcf.gz --indv trichocarpa1 --freq --out $OutDir_outgroup/est_sfs/pti.0fold.est_sfs


elif [ $step == "5" ]; then
##use perl script to create the input file of est_sfs
export PERL5LIB=/data/apps/cpan/lib/perl5:/data/apps/cpan/lib64/perl5:/data/apps/cpan/share/perl5
perl est_sfs_input.dfe.${species}.pl $OutDir_est_sfs/${species}.snp.4fold.frq $OutDir_outgroup/est_sfs/peu.4fold.est_sfs.frq $OutDir_outgroup/est_sfs/pti.4fold.est_sfs.frq
perl est_sfs_input.dfe.${species}.pl $OutDir_est_sfs/${species}.snp.0fold.frq $OutDir_outgroup/est_sfs/peu.0fold.est_sfs.frq $OutDir_outgroup/est_sfs/pti.0fold.est_sfs.frq


elif [ $step == "6" ]; then
##running est_sfs to infer the ancestral state and unfold sfs file

#rate6
$est_sfs $config_rate6 $OutDir_est_sfs/${species}.snp.4fold.input.txt $seed $OutDir_est_sfs/${species}.snp.4fold.est_sfs.output.txt $OutDir_est_sfs/${species}.snp.4fold.est_sfs.output.p_anc.txt

###because there is core dump problem for 0fold sites, we randomly sampled 5322249 sites (the same number as 4fold sites) to run the analysis
shuf -n 5322249 $OutDir_est_sfs/${species}.snp.0fold.input.txt > $OutDir_est_sfs/${species}.snp.0fold.random.input.txt
$est_sfs $config_rate6 $OutDir_est_sfs/${species}.snp.0fold.random.input.txt $seed $OutDir_est_sfs/${species}.snp.0fold.est_sfs.random.output.txt $OutDir_est_sfs/${species}.snp.0fold.est_sfs.random.output.p_anc.txt

elif [ $step == "7" ]; then
##create the input file for dfe-alpha
n=$(cat /w/user260/multiple_aspen/populus190_vcf/pixy/samples/${species}.txt |wc -l)
alleles=$(echo ''$n'*2'| bc)
echo $alleles

echo -e "1" > $OutDir_dfe/${species}.dfe.input.txt
echo -e "$alleles" >> $OutDir_dfe/${species}.dfe.input.txt
sed 's/,/ /g' $OutDir_est_sfs/${species}.snp.0fold.est_sfs.random.output.txt >> $OutDir_dfe/${species}.dfe.input.txt
sed 's/,/ /g' $OutDir_est_sfs/${species}.snp.4fold.est_sfs.output.txt >> $OutDir_dfe/${species}.dfe.input.txt
#sed 's/,/\n/g' $OutDir_est_sfs/${species}.snp.0fold.est_sfs.random.output.txt | sed 's/\..*//g' |sed ':a;N;$!ba;s/\n/ /g'  >> $OutDir_dfe/${species}.dfe.input.txt
#sed 's/,/\n/g' $OutDir_est_sfs/${species}.snp.4fold.est_sfs.output.txt | sed 's/\..*//g' |sed ':a;N;$!ba;s/\n/ /g' >> $OutDir_dfe/${species}.dfe.input.txt


elif [ $step == "7.1" ]; then
##run est_dfe on site_class 0

if [ ! -d "$OutDir_dfe/results_dir_neut" ]; then
mkdir -p $OutDir_dfe/results_dir_neut
fi

echo "data_path_1    /UserHome/user260/tools/dfe_alpha/data" > $OutDir_dfe/${species}.dfe.site_class-0.txt
echo "data_path_2    /UserHome/user260/tools/dfe_alpha/data/data-three-epoch" >> $OutDir_dfe/${species}.dfe.site_class-0.txt
echo "sfs_input_file    $OutDir_dfe/${species}.dfe.input.txt" >> $OutDir_dfe/${species}.dfe.site_class-0.txt
echo "est_dfe_results_dir    $OutDir_dfe/results_dir_neut" >> $OutDir_dfe/${species}.dfe.site_class-0.txt
echo "site_class    0" >> $OutDir_dfe/${species}.dfe.site_class-0.txt
echo "fold    1" >> $OutDir_dfe/${species}.dfe.site_class-0.txt
echo "epochs    2"  >> $OutDir_dfe/${species}.dfe.site_class-0.txt
echo "search_n2    1" >> $OutDir_dfe/${species}.dfe.site_class-0.txt
echo "t2_variable    1" >> $OutDir_dfe/${species}.dfe.site_class-0.txt
echo "t2    50" >> $OutDir_dfe/${species}.dfe.site_class-0.txt

$est_dfe -c $OutDir_dfe/${species}.dfe.site_class-0.txt

elif [ $step == "7.2" ]; then
##run est_dfe on site_class 1
if [ ! -d "$OutDir_dfe/results_dir_sel" ]; then
mkdir -p $OutDir_dfe/results_dir_sel
fi

echo "data_path_1    /UserHome/user260/tools/dfe_alpha/data" > $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "data_path_2    /UserHome/user260/tools/dfe_alpha/data/data-three-epoch" >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "sfs_input_file    $OutDir_dfe/${species}.dfe.input.txt" >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "est_dfe_results_dir    $OutDir_dfe/results_dir_sel" >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "est_dfe_demography_results_file    $OutDir_dfe/results_dir_neut/est_dfe.out" >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "site_class    1" >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "fold    1" >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "epochs    2"  >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "mean_s_variable    1" >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "mean_s    -0.1" >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "beta_variable    1" >> $OutDir_dfe/${species}.dfe.site_class-1.txt
echo "beta    0.5" >> $OutDir_dfe/${species}.dfe.site_class-1.txt

$est_dfe -c $OutDir_dfe/${species}.dfe.site_class-1.txt
$prop_muts_in_s_ranges -c $OutDir_dfe/results_dir_sel/est_dfe.out -o $OutDir_dfe/results_dir_sel/prop_muls_in_s_ranges.output

elif [ $step == "7.3" ]; then
##run est_alpha_omega

echo "data_path_1    /UserHome/user260/tools/dfe_alpha/data" >  $OutDir_dfe/${species}.dfe.alpha_omega.txt
echo "divergence_file    $OutDir_dfe/${species}.omega.txt" >> $OutDir_dfe/${species}.dfe.alpha_omega.txt
echo "est_alpha_omega_results_file    $OutDir_dfe/est_alpha_omega.${species}.out" >> $OutDir_dfe/${species}.dfe.alpha_omega.txt
echo "est_dfe_results_file    $OutDir_dfe/results_dir_sel/est_dfe.out"  >> $OutDir_dfe/${species}.dfe.alpha_omega.txt
echo "neut_egt_file    $OutDir_dfe/results_dir_neut/neut_egf.out" >> $OutDir_dfe/${species}.dfe.alpha_omega.txt
echo "sel_egf_file    $OutDir_dfe/results_dir_sel/sel_egf.out"  >> $OutDir_dfe/${species}.dfe.alpha_omega.txt
echo "do_jukes_cantor    1" >> $OutDir_dfe/${species}.dfe.alpha_omega.txt
echo "remove_poly    0" >> $OutDir_dfe/${species}.dfe.alpha_omega.txt

$est_alpha_omega -c $OutDir_dfe/${species}.dfe.alpha_omega.txt


elif [ $step == "8" ]; then
###running bootstrap 200 times

#rate6
#4fold

boot_n=200

i=$3
echo $i

OutDir_boot_n=$OutDir_est_sfs/bootstrap/bootstrap$i
if [ ! -d "$OutDir_boot_n" ]; then
mkdir -p $OutDir_boot_n
fi

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/UserHome/user260/tools/dfe_alpha/

#4fold
shuf -r -n 5322249 $OutDir_est_sfs/${species}.snp.4fold.input.txt > $OutDir_boot_n/${species}.snp.4fold.bootstrap$i.input.txt
$est_sfs $config_rate6 $OutDir_boot_n/${species}.snp.4fold.bootstrap$i.input.txt $seed $OutDir_boot_n/${species}.snp.4fold.bootstrap$i.est_sfs.output.txt $OutDir_boot_n/${species}.snp.4fold.bootstrap$i.est_sfs.p_anc.output.txt
rm $OutDir_boot_n/${species}.snp.4fold.bootstrap$i.input.txt $OutDir_boot_n/${species}.snp.4fold.bootstrap$i.est_sfs.p_anc.output.txt

#0fold
shuf -r -n 5322249 $OutDir_est_sfs/${species}.snp.0fold.input.txt > $OutDir_boot_n/${species}.snp.0fold.bootstrap$i.input.txt
$est_sfs $config_rate6 $OutDir_boot_n/${species}.snp.0fold.bootstrap$i.input.txt $seed $OutDir_boot_n/${species}.snp.0fold.bootstrap$i.est_sfs.output.txt $OutDir_boot_n/${species}.snp.0fold.bootstrap$i.est_sfs.p_anc.output.txt
rm $OutDir_boot_n/${species}.snp.0fold.bootstrap$i.input.txt $OutDir_boot_n/${species}.snp.0fold.bootstrap$i.est_sfs.p_anc.output.txt

#input for dfe_alpha
n=$(cat /w/user260/multiple_aspen/populus190_vcf/pixy/samples/${species}.txt |wc -l)
alleles=$(echo ''$n'*2'| bc)
echo $alleles

echo -e "1" > $OutDir_boot_n/${species}.dfe.input.txt
echo -e "$alleles" >> $OutDir_boot_n/${species}.dfe.input.txt
sed 's/,/ /g' $OutDir_boot_n/${species}.snp.0fold.bootstrap$i.est_sfs.output.txt >> $OutDir_boot_n/${species}.dfe.input.txt
sed 's/,/ /g' $OutDir_boot_n/${species}.snp.4fold.bootstrap$i.est_sfs.output.txt >> $OutDir_boot_n/${species}.dfe.input.txt

###run est_dfe on site_class 0


if [ ! -d "$OutDir_boot_n/results_dir_neut" ]; then
mkdir -p $OutDir_boot_n/results_dir_neut
fi

echo "data_path_1    /UserHome/user260/tools/dfe_alpha/data" > $OutDir_boot_n/${species}.dfe.site_class-0.txt
echo "data_path_2    /UserHome/user260/tools/dfe_alpha/data/data-three-epoch" >> $OutDir_boot_n/${species}.dfe.site_class-0.txt
echo "sfs_input_file    $OutDir_boot_n/${species}.dfe.input.txt" >> $OutDir_boot_n/${species}.dfe.site_class-0.txt
echo "est_dfe_results_dir    $OutDir_boot_n/results_dir_neut" >> $OutDir_boot_n/${species}.dfe.site_class-0.txt
echo "site_class    0" >> $OutDir_boot_n/${species}.dfe.site_class-0.txt
echo "fold    1" >> $OutDir_boot_n/${species}.dfe.site_class-0.txt
echo "epochs    2"  >> $OutDir_boot_n/${species}.dfe.site_class-0.txt
echo "search_n2    1" >> $OutDir_boot_n/${species}.dfe.site_class-0.txt
echo "t2_variable    1" >> $OutDir_boot_n/${species}.dfe.site_class-0.txt
echo "t2    50" >> $OutDir_boot_n/${species}.dfe.site_class-0.txt

$est_dfe -c $OutDir_boot_n/${species}.dfe.site_class-0.txt

##run est_dfe on site_class 1
if [ ! -d "$OutDir_boot_n/results_dir_sel" ]; then
mkdir -p $OutDir_boot_n/results_dir_sel
fi

echo "data_path_1    /UserHome/user260/tools/dfe_alpha/data" > $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "data_path_2    /UserHome/user260/tools/dfe_alpha/data/data-three-epoch" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "sfs_input_file    $OutDir_boot_n/${species}.dfe.input.txt" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "est_dfe_results_dir    $OutDir_boot_n/results_dir_sel" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "est_dfe_demography_results_file    $OutDir_boot_n/results_dir_neut/est_dfe.out" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "site_class    1" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "fold    1" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "epochs    2"  >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "mean_s_variable    1" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "mean_s    -0.1" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "beta_variable    1" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt
echo "beta    0.5" >> $OutDir_boot_n/${species}.dfe.site_class-1.txt

$est_dfe -c $OutDir_boot_n/${species}.dfe.site_class-1.txt
$prop_muts_in_s_ranges -c $OutDir_boot_n/results_dir_sel/est_dfe.out -o $OutDir_boot_n/results_dir_sel/prop_muls_in_s_ranges.output


elif [ $step == "9" ]; then
##summarize the results
echo -e "name\tN1\tN2\tt2\tNw\tb\tEs\tf0\tL\tNes_1\tNes_1_10\tNes_10_100\tNes_100" > $OutDir_summary/${species}.dfe.summary.txt

##the real data
N1=$(cut -f 2 -d " " $OutDir_dfe/results_dir_sel/est_dfe.out)
N2=$(cut -f 4 -d " " $OutDir_dfe/results_dir_sel/est_dfe.out)
t2=$(cut -f 6 -d " " $OutDir_dfe/results_dir_sel/est_dfe.out)
Nw=$(cut -f 8 -d " " $OutDir_dfe/results_dir_sel/est_dfe.out)
b=$(cut -f 10 -d " " $OutDir_dfe/results_dir_sel/est_dfe.out)
Es=$(cut -f 12 -d " " $OutDir_dfe/results_dir_sel/est_dfe.out)
f0=$(cut -f 14 -d " " $OutDir_dfe/results_dir_sel/est_dfe.out)
L=$(cut -f 16 -d " " $OutDir_dfe/results_dir_sel/est_dfe.out)
Nes_1=$(cut -f 3 -d " "  $OutDir_dfe/results_dir_sel/prop_muls_in_s_ranges.output)
Nes_1_10=$(cut -f 6 -d " "  $OutDir_dfe/results_dir_sel/prop_muls_in_s_ranges.output)
Nes_10_100=$(cut -f 9 -d " "  $OutDir_dfe/results_dir_sel/prop_muls_in_s_ranges.output)
Nes_100=$(cut -f 12 -d " "  $OutDir_dfe/results_dir_sel/prop_muls_in_s_ranges.output)

echo -e "real\t$N1\t$N2\t$t2\t$Nw\t$b\t$Es\t$f0\t$L\t$Nes_1\t$Nes_1_10\t$Nes_10_100\t$Nes_100" >> $OutDir_summary/${species}.dfe.summary.txt

##bootstrap
for i in {1..200}
do
OutDir_boot_n=$OutDir_est_sfs/bootstrap/bootstrap$i

N1=$(cut -f 2 -d " " $OutDir_boot_n/results_dir_sel/est_dfe.out)
N2=$(cut -f 4 -d " " $OutDir_boot_n/results_dir_sel/est_dfe.out)
t2=$(cut -f 6 -d " " $OutDir_boot_n/results_dir_sel/est_dfe.out)
Nw=$(cut -f 8 -d " " $OutDir_boot_n/results_dir_sel/est_dfe.out)
b=$(cut -f 10 -d " " $OutDir_boot_n/results_dir_sel/est_dfe.out)
Es=$(cut -f 12 -d " " $OutDir_boot_n/results_dir_sel/est_dfe.out)
f0=$(cut -f 14 -d " " $OutDir_boot_n/results_dir_sel/est_dfe.out)
L=$(cut -f 16 -d " " $OutDir_boot_n/results_dir_sel/est_dfe.out)
Nes_1=$(cut -f 3 -d " "  $OutDir_boot_n/results_dir_sel/prop_muls_in_s_ranges.output)
Nes_1_10=$(cut -f 6 -d " "  $OutDir_boot_n/results_dir_sel/prop_muls_in_s_ranges.output)
Nes_10_100=$(cut -f 9 -d " "  $OutDir_boot_n/results_dir_sel/prop_muls_in_s_ranges.output)
Nes_100=$(cut -f 12 -d " "  $OutDir_boot_n/results_dir_sel/prop_muls_in_s_ranges.output)


echo -e "bootstrap\t$N1\t$N2\t$t2\t$Nw\t$b\t$Es\t$f0\t$L\t$Nes_1\t$Nes_1_10\t$Nes_10_100\t$Nes_100" >> $OutDir_summary/${species}.dfe.summary.txt
done

fi


