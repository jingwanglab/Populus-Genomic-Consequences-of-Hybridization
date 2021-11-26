#!/bin/sh
#SBATCH -c 1 --mem 8G
#SBATCH -o /UserHome/user260/pipeline/selection/selscan/selscan.populus.out
#SBATCH -e /UserHome/user260/pipeline/selection/selscan/selscan.populus.err

vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
selscan="/UserHome/user260/tools/selscan-master/bin/linux/selscan"
norm="/UserHome/user260/tools/selscan-master/bin/linux/norm"



vcf="/w/user260/multiple_aspen/populus190_vcf/beagle/populus162.phased.recode.vcf.gz"

species_samples="/w/user260/multiple_aspen/populus190_vcf/ldhat/species_samples"


OutDir="/w/user260/multiple_aspen/populus190_vcf/beagle"

OutDir_selscan=$OutDir/selscan
if [ ! -d "$OutDir_selscan" ]; then
mkdir -p $OutDir_selscan
fi

step=$1
species=$2
OutDir_species=$OutDir_selscan/$species
if [ ! -d "$OutDir_species" ]; then
mkdir -p $OutDir_species
fi

OutDir_out=$OutDir_selscan/$species/output
if [ ! -d "$OutDir_out" ]; then
mkdir -p $OutDir_out
fi

OutDir_bed=$OutDir_out/bed
if [ ! -d "$OutDir_bed" ]; then
mkdir -p $OutDir_bed
fi


Out=${vcf##*/}
OutSuffix=${Out%.vcf.gz}.$species
echo $OutSuffix


non_beagle_vcf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/${species}.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz"

if [ $step == "1" ]; then

#$vcftools --gzvcf $vcf --keep $species_samples/${species}.txt --maf 0.000001 --recode --recode-INFO-all --out $OutDir_species/$OutSuffix && $pigz $OutDir_species/$OutSuffix.recode.vcf

#zcat $non_beagle_vcf|grep -v "#" |cut -f 3 > $OutDir_species/$species.snp.txt
#$vcftools --gzvcf $OutDir_species/$OutSuffix.recode.vcf.gz --snps $OutDir_species/$species.snp.txt --recode --recode-INFO-all --out $OutDir_species/$OutSuffix.snp

chr=$3
$vcftools --vcf $OutDir_species/$OutSuffix.snp.recode.vcf --chr $chr --recode --recode-INFO-all --out $OutDir_species/$OutSuffix.snp.$chr
awk '{sub("chr", "")}1' $OutDir_species/$OutSuffix.snp.$chr.recode.vcf > $OutDir_species/$OutSuffix.snp.$chr.temp && mv $OutDir_species/$OutSuffix.snp.$chr.temp $OutDir_species/$OutSuffix.snp.$chr.recode.vcf
$vcftools --vcf $OutDir_species/$OutSuffix.snp.$chr.recode.vcf --plink  --out $OutDir_species/$OutSuffix.snp.$chr
awk '$3=$4/10000' $OutDir_species/$OutSuffix.snp.$chr.map |sed 's/ /\t/g' > $OutDir_species/$OutSuffix.snp.$chr.map.temp && mv $OutDir_species/$OutSuffix.snp.$chr.map.temp $OutDir_species/$OutSuffix.snp.$chr.map

##running the selscan
$selscan --ihh12 --vcf $OutDir_species/$OutSuffix.snp.$chr.recode.vcf --map $OutDir_species/$OutSuffix.snp.$chr.map --out $OutDir_out/$OutSuffix.snp.$chr

elif [ $step == "2" ]; then
#normalize the values

input=
for chr in chr{1..19}
do
input="$input $OutDir_out/$OutSuffix.snp.$chr.ihh12.out"
done

$norm --ihh12 --files $input


elif [ $step == "3" ]; then
##change the normalize output into bed file
#ihh12
for chr in chr{1..19}
do
sed '1d' $OutDir_out/$OutSuffix.snp.$chr.ihh12.out.norm |sed 's/:/\t/g' |awk '$2=$2-1' | sed 's/ /\t/g' > $OutDir_out/$OutSuffix.snp.$chr.ihh12.out.norm.temp
done

input=
for chr in chr{1..19}
do
input="$input $OutDir_out/$OutSuffix.snp.$chr.ihh12.out.norm.temp"
done

#cat $input | awk 'function abs(x){return (x < 0) ? -x : x;} {print $1,$2,$3,abs($6)}' |sed 's/ /\t/g'  > $OutDir_bed/$OutSuffix.snp.all.ihh12.out.100bins.norm.bed && rm $OutDir_out/*temp
cat $input |cut -f 1,2,3,6 > $OutDir_bed/$OutSuffix.snp.all.ihh12.out.100bins.norm.bed && rm $OutDir_out/*temp

fi
