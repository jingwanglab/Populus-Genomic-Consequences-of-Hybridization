#!/bin/sh
#SBATCH -c 1 --mem 10G
#SBATCH -o /UserHome/user260/pipeline/hybridization/volcanofinder/volcanofinder.populus.out
#SBATCH -e /UserHome/user260/pipeline/hybridization/volcanofinder/volcanofinder.populus.err

###Main aim:
#1.outgroup vcf, remove missing and intersect with mask.bed, to extract the sites without missing in both outgroups
#2.create the input of est_sfs to infer the ancestral and derived allele and esimate the unnormalized site frequency spectrum
#3.create the input of allele frequency file


vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
est_sfs="/UserHome/user260/tools/est-sfs-release-2.03/est-sfs"
plink2="/UserHome/user260/tools/plink/plink2"
volcanofinder="/UserHome/user260/tools/volcanofinder_v1.0/VolcanoFinder"

config_jc="/UserHome/user260/tools/est-sfs-release-2.03/config-JC.txt"
config_kimura="/UserHome/user260/tools/est-sfs-release-2.03/config-kimura.txt"
config_rate6="/UserHome/user260/tools/est-sfs-release-2.03/config-rate6.txt"
seed="/UserHome/user260/tools/est-sfs-release-2.03/seedfile.txt"


peu_vcf="/w/user260/multiple_aspen/outgroups/PeuXBY08.snp.vcf.gz"
pti_vcf="/w/user260/multiple_aspen/outgroups/trichocarpa1.snp.vcf.gz"

outgroup_filter="/w/user260/multiple_aspen/outgroups/filter"
if [ ! -d "$outgroup_filter" ]; then
mkdir -p $outgroup_filter
fi

mask_bed="/w/user260/multiple_aspen/snpable/Potra02_genome.mask_mappability.no_scaffold.bed"

populus162="/w/user260/multiple_aspen/populus162_vcf/vcf/populus162.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz"

volcan_output="/w/user260/multiple_aspen/populus162_vcf/vcf/volcanofinder"


step=$1
species=$2
chr=$3

species_sample="/w/user260/multiple_aspen/populus162_vcf/vcf/volcanofinder/samples/${species}.txt"

volcan_species=$volcan_output/$species
if [ ! -d "$volcan_species" ]; then
mkdir -p $volcan_species
fi

volcan_outgroup=$volcan_output/outgroup
if [ ! -d "$volcan_outgroup" ]; then
mkdir -p $volcan_outgroup
fi

volcan_species_run=$volcan_output/$species/run
if [ ! -d "$volcan_species_run" ]; then
mkdir -p $volcan_species_run
fi

volcan_species_chr=$volcan_output/$species/run/$chr
if [ ! -d "$volcan_species_chr" ]; then
mkdir -p $volcan_species_chr
fi


if [ $step == "1" ]; then
##create the no miss file
$bcftools view -g ^miss -R $mask_bed -M2 -V indels -O z -o $outgroup_filter/PeuXBY08.snp.no_miss.mask.vcf.gz $peu_vcf

elif [ $step == "2" ]; then

$bcftools view -g ^miss -R $mask_bed  -M2 -V indels -O z -o $outgroup_filter/trichocarpa1.snp.no_miss.mask.vcf.gz $pti_vcf

elif [ $step == "3" ]; then

$vcftools --gzvcf $populus162 --keep $species_sample --freq2 --out $volcan_species/$species 
$vcftools --gzvcf $populus162 --keep $species_sample --counts2 --out $volcan_species/$species 

elif [ $step == "4" ]; then
#merge samples
#$bcftools view /w/user260/multiple_aspen/populus162_vcf/vcf/populus162.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf -O z -o $populus162 
# rm /w/user260/multiple_aspen/populus162_vcf/vcf/populus162.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf
#$bcftools index $populus162
#$bcftools index $outgroup_filter/PeuXBY08.snp.no_miss.mask.vcf.gz
#$bcftools index $outgroup_filter/trichocarpa1.snp.no_miss.mask.vcf.gz

#$bcftools annotate -x FORMAT -O z -o $volcan_output/populus162.outgroup.clean.vcf.gz $populus162 
$bcftools index $volcan_output/populus162.outgroup.clean.vcf.gz

$bcftools merge -O z -o $volcan_output/populus162.outgroup.merge.vcf.gz $volcan_output/populus162.outgroup.clean.vcf.gz $outgroup_filter/PeuXBY08.snp.no_miss.mask.vcf.gz $outgroup_filter/trichocarpa1.snp.no_miss.mask.vcf.gz


elif [ $step == "5" ]; then
#extract the freq of each species 
$vcftools --gzvcf $volcan_output/populus162.outgroup.merge.vcf.gz --keep $species_sample --freq --out $volcan_species/${species}.merge

elif [ $step == "5.1" ]; then
#outgroup
$vcftools --gzvcf $volcan_output/populus162.outgroup.merge.vcf.gz --indv PeuXBY08 --freq --out $volcan_outgroup/PeuXBY08.merge
$vcftools --gzvcf $volcan_output/populus162.outgroup.merge.vcf.gz --indv trichocarpa1 --freq --out $volcan_outgroup/trichocarpa1.merge


elif [ $step == "6" ]; then
#extract the SNPs with 2 alleles
awk '$3<3' $volcan_species/${species}.merge.frq > $volcan_species/temp && mv $volcan_species/temp $volcan_species/${species}.merge.frq


elif [ $step == "6.1" ]; then
awk '$3<3' $volcan_outgroup/PeuXBY08.merge.frq > $volcan_outgroup/PeuXBY08.temp && mv $volcan_outgroup/PeuXBY08.temp $volcan_outgroup/PeuXBY08.merge.frq
awk '$3<3' $volcan_outgroup/trichocarpa1.merge.frq > $volcan_outgroup/trichocarpa1.temp && mv $volcan_outgroup/trichocarpa1.temp $volcan_outgroup/trichocarpa1.merge.frq

elif [ $step == "7" ]; then

perl est_sfs_input.volcanofinder.${species}.pl $volcan_species/${species}.merge.frq 

elif [ $step == "7.1" ]; then

perl est_sfs_input.volcanofinder.outgroup.pl $volcan_outgroup/trichocarpa1.merge.frq 
perl est_sfs_input.volcanofinder.outgroup.pl $volcan_outgroup/PeuXBY08.merge.frq

elif [ $step == "8" ]; then
#use est_sfs to estimate the ancestral/derived state of 
#paste $volcan_species/${species}.merge.input.txt $volcan_outgroup/trichocarpa1.merge.input.txt $volcan_outgroup/PeuXBY08.merge.input.txt > $volcan_species/${species}.merge.est_sfs.input.txt
shuf -n 1000000 $volcan_species/${species}.merge.est_sfs.input.txt > $volcan_species/${species}.shuf.est_sfs.input.txt

$est_sfs $config_rate6 $volcan_species/${species}.shuf.est_sfs.input.txt $seed $volcan_species/${species}.shuf.est_sfs.output.txt $volcan_species/${species}.shuf.est_sfs.output.p_anc.txt && rm $volcan_species/${species}.shuf.est_sfs.output.p_anc.txt

elif [ $step == "9.1" ]; then

awk '$3=="2"' $volcan_outgroup/trichocarpa1.merge.frq |sed 's/:/\t/g' |cut -f 6,8 >  $volcan_outgroup/trichocarpa1.merge.temp
awk '$3=="2"' $volcan_outgroup/PeuXBY08.merge.frq |sed 's/:/\t/g'|cut -f 6,8 >  $volcan_outgroup/PeuXBY08.merge.temp


elif [ $step == "9" ]; then
#extract the allele frequency file
#1.fixed snp-id of outgroup trichocarpa
#2.fixed snp-id of outgroup peu
#3.polymorphic snp of each species
#awk '$3=="2"' $volcan_species/${species}.merge.frq | sed 's/:/\t/g' > $volcan_species/${species}.merge.outgroup.temp

#paste $volcan_species/${species}.merge.outgroup.temp $volcan_outgroup/trichocarpa1.merge.temp $volcan_outgroup/PeuXBY08.merge.temp >  $volcan_species/${species}.merge.outgroup.frq && rm $volcan_species/${species}.merge.outgroup.temp 

#sed 's/-nan/1/g' $volcan_species/${species}.merge.outgroup.frq > $volcan_species/${species}.merge.outgroup.frq2

#perl volcanofinder.alle_freq.input1.pl $volcan_species/${species}.merge.outgroup.frq $volcan_species/${species}.merge.outgroup.frq2

#only extract the polymorhpic sites for the analyses
awk '$6!="-nan"' $volcan_species/${species}.merge.outgroup.input.txt |awk '$6!="0"' |awk '$6!="1"' > $volcan_species/${species}.merge.outgroup.poly.input.txt

perl volcanofinder.alle_freq.est_sfs.${species}.pl $volcan_species/${species}.merge.outgroup.poly.input.txt

elif [ $step == "10" ]; then

split -l 50000 $volcan_species/${species}.merge.outgroup.poly.est_sfs.input.txt $volcan_species/${species}.merge.outgroup.poly.est_sfs.split.

for x in $volcan_species/${species}.merge.outgroup.poly.est_sfs.split.*
do
$est_sfs $config_rate6 $x $seed $x.output.txt $x.output.p_anc.txt
done

#$est_sfs $config_rate6 $volcan_species/${species}.merge.outgroup.poly.est_sfs.input.txt $seed $volcan_species/${species}.merge.outgroup.poly.est_sfs.input.output.txt $volcan_species/${species}.merge.outgroup.poly.est_sfs.input.output.p_anc.txt

elif [ $step == "11" ]; then

input=
for x in $volcan_species/${species}.merge.outgroup.poly.est_sfs.split.*.output.p_anc.txt
do
grep -v "^0" $x > $x.temp
input="$input $x.temp"
done
cat $input > $volcan_species/${species}.merge.outgroup.poly.est_sfs.output.p_anc.txt && rm $input

perl est_sfs_ancestral_derived.alle_freq.${species}.pl $volcan_species/${species}.merge.outgroup.poly.input.txt $volcan_species/${species}.merge.outgroup.poly.est_sfs.output.p_anc.txt

mv $volcan_species/${species}.merge.outgroup.poly.volcanofinder.input.txt $volcan_species_run/

#elif [ $step == "12" ]; then
##input to run the volcan

sed 's/,/\n/g' $volcan_species/${species}.shuf.est_sfs.output.txt | awk '{printf("%f\n",$1/1000000)}' |sed '1d' > $volcan_species/${species}.shuf.est_sfs.output.temp

n=$(cat $volcan_species/${species}.shuf.est_sfs.output.temp |wc -l)
n2=$(echo "$n-1" | bc)


for ((i=1;i<=$n;i++))
do
echo $i >> $volcan_species/${species}.shuf.est_sfs.output.temp2
done

paste $volcan_species/${species}.shuf.est_sfs.output.temp2 $volcan_species/${species}.shuf.est_sfs.output.temp | head -n $n2 > $volcan_species_run/${species}.shuf.volcan.sfs.input.txt&& rm $volcan_species/${species}.shuf.est_sfs.output.temp2 $volcan_species/${species}.shuf.est_sfs.output.temp

echo "5" > $volcan_species_run/${species}.d.grid.txt
echo "0.00065" >> $volcan_species_run/${species}.d.grid.txt
echo "0.00075" >> $volcan_species_run/${species}.d.grid.txt
echo "0.00085" >> $volcan_species_run/${species}.d.grid.txt
echo "0.00095" >> $volcan_species_run/${species}.d.grid.txt
echo "0.00102" >> $volcan_species_run/${species}.d.grid.txt


elif [ $step == "13" ]; then
#run for each chromosome

awk '$1=="'$chr'"' $volcan_species_run/${species}.merge.outgroup.poly.volcanofinder.input.txt |cut -f 2- > $volcan_species_chr/${species}.$chr.alle_freq.temp

echo -e "position\tx\tn\tfolded" > $volcan_species_chr/${species}.$chr.alle_freq.txt
cat $volcan_species_chr/${species}.$chr.alle_freq.temp >> $volcan_species_chr/${species}.$chr.alle_freq.txt && rm $volcan_species_chr/${species}.$chr.alle_freq.temp

cut -f 1 $volcan_species_chr/${species}.$chr.alle_freq.txt > $volcan_species_chr/${species}.$chr.grid.txt

$volcanofinder -ig 10000 $volcan_species_chr/${species}.$chr.alle_freq.txt $volcan_species_run/${species}.shuf.volcan.sfs.input.txt -1 1 1 $volcan_species_chr/${species}.$chr.volcan.out.txt

elif [ $step == "14" ]; then
#summarize the results across the chromosomes

input=
for chr in chr{1..19}
do
volcan_species_chr=$volcan_output/$species/run/$chr
awk '{$5="'$chr'"}1' $volcan_species_chr/${species}.$chr.volcan.out.txt |awk '{print $5,$1,$2,$3,$4}' |sed '1d' > $volcan_species_chr/${species}.$chr.volcan.out.new.txt
input="$input $volcan_species_chr/${species}.$chr.volcan.out.new.txt"
done

cat $input > $volcan_species_run/${species}.volcan.out.txt

elif [ $step == "15" ]; then
#change the formate to bed

perl volcanofinder.summarize.pl $volcan_species_run/${species}.volcan.out.txt

cp $volcan_species_run/${species}.volcan.out.bed /w/user260/multiple_aspen/populus162_vcf/vcf/volcanofinder/summary

fi



















