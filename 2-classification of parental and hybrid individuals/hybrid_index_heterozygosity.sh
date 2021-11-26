#!/bin/sh
#SBATCH -c 1 --mem 8G
#SBATCH -o /UserHome/user260/pipeline/hybridization/hybrid_index/hybrid_index_heterozygosity.out
#SBATCH -e /UserHome/user260/pipeline/hybridization/hybrid_index/hybrid_index_heterozygosity.err


#####this script is used to calcualte the hybrid index and heterzygosity for hybrid individuals
###The steps includes:
#1. Extract individuals for species1 (parent1), species2(parent2) and their hybrids
#2. only retain the SNPs that are fixed between species1 and species2 (Fst=1)
#3. use script to calculate the hybrid index and heterozygosity for the hybrid individuals

vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
bgzip="/UserHome/user260/tools/bcftools-1.2/htslib-1.2.1/bgzip"
beagle="/data/apps/beagle/beagle.08Jun17.d8b.jar"
java="/data/apps/jdk/jdk1.8.0_131/bin/java"


InputDir="/w/user260/multiple_aspen/populus227_vcf"
populus227_vcf="$InputDir/populus227.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz"

OutDir=$InputDir/hybrid_index

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

samples="$OutDir/samples/populus227.samples.remove_clones.txt"

step=$1 ##the step to perform
species1=$2 #the potential parent1
species2=$3  #the potential parent2
hybrids=$4  ##the hybrids

OutDir2=$OutDir/${species1}_${species2}

if [ ! -d "$OutDir2" ]; then
mkdir -p $OutDir2
fi

OutDir3=$OutDir/${species1}_${species2}/output

if [ ! -d "$OutDir3" ]; then
mkdir -p $OutDir3
fi

##########################################################
##########################################################
###the step1 is to calculate the fst between species 1 and species 2 
if [ $step == "1" ]; then

##species1: tremula_china
##species2: alba
##hybrids: alba_tremula
awk -F "\t" '{if($9~/'$species1'/ && $7~/Pure/)print $4}' $samples > $OutDir2/${species1}.txt
awk -F "\t" '{if($9~/'$species2'/ && $7~/Pure/)print $4}' $samples > $OutDir2/${species2}.txt
awk -F "\t" '{if($9~/'$hybrids'/ && $7~/hybrids/)print $4}' $samples > $OutDir2/${hybrids}.txt

$vcftools --gzvcf $populus227_vcf --weir-fst-pop $OutDir2/${species1}.txt --weir-fst-pop $OutDir2/${species2}.txt --out $OutDir2/${species1}.${species2}

##################extract the SNPs with Fst=1 
awk '$3=="1"' $OutDir2/${species1}.${species2}.weir.fst|cut -f 1,2 |sed 's/\t/:/g' > $OutDir2/${species1}.${species2}.fst1.snp.txt

####extracting these SNPs among the two parental species and their hybrids and only consider those with missing proportion <20%
$vcftools --gzvcf $populus227_vcf --snps $OutDir2/${species1}.${species2}.fst1.snp.txt --keep $OutDir2/${species1}.txt --keep $OutDir2/${species2}.txt --keep $OutDir2/${hybrids}.txt --max-missing 0.8 --maf 0.0000001 --recode --recode-INFO-all --out $OutDir2/${species1}.${species2}.${hybrids}


#elif [ $step == "2" ]; then
#################################################
grep -v "#" $OutDir2/${species1}.${species2}.${hybrids}.recode.vcf | cut -f 3 > $OutDir2/${species1}.${species2}.fst1.miss20.snp.txt

$vcftools --gzvcf $populus227_vcf --snps $OutDir2/${species1}.${species2}.fst1.miss20.snp.txt --keep $OutDir2/${species1}.txt --recode --recode-INFO-all --out $OutDir2/${species1}
$vcftools --gzvcf $populus227_vcf --snps $OutDir2/${species1}.${species2}.fst1.miss20.snp.txt --keep $OutDir2/${species2}.txt --recode --recode-INFO-all --out $OutDir2/${species2}
$vcftools --gzvcf $populus227_vcf --snps $OutDir2/${species1}.${species2}.fst1.miss20.snp.txt --keep $OutDir2/${hybrids}.txt --recode --recode-INFO-all --out $OutDir2/${hybrids}


$bcftools annotate -x INFO,^FORMAT/GT -o $OutDir2/${species1}.gt.vcf $OutDir2/${species1}.recode.vcf
$bcftools annotate -x INFO,^FORMAT/GT -o $OutDir2/${species2}.gt.vcf $OutDir2/${species2}.recode.vcf
$bcftools annotate -x INFO,^FORMAT/GT -o $OutDir2/${hybrids}.gt.vcf $OutDir2/${hybrids}.recode.vcf

####output the allele frequency for the two parent species to determine the status of the allele
$vcftools --vcf $OutDir2/${species1}.gt.vcf --freq2 --out $OutDir2/${species1} 
$vcftools --vcf $OutDir2/${species2}.gt.vcf --freq2 --out $OutDir2/${species2} 


#elif [ $step == "3" ]; then
#################################################

##species1
cut -f 6 $OutDir2/${species1}.frq |sed '1d' > $OutDir3/${species1}.non_ref.freq
paste $OutDir3/${species1}.non_ref.freq $OutDir3/${species1}.non_ref.freq |sed 's/\t/\//g' > $OutDir3/${species1}.temp1
echo -e "$species1" > $OutDir3/${species1}.temp2
cat $OutDir3/${species1}.temp2 $OutDir3/${species1}.temp1 > $OutDir3/${species1}.geno.txt && rm $OutDir3/${species1}.non_ref.freq $OutDir3/${species1}.temp1 $OutDir3/${species1}.temp2

##species2
cut -f 6 $OutDir2/${species2}.frq |sed '1d' > $OutDir3/${species2}.non_ref.freq
paste $OutDir3/${species2}.non_ref.freq $OutDir3/${species2}.non_ref.freq |sed 's/\t/\//g' > $OutDir3/${species2}.temp1
echo -e "$species2" > $OutDir3/${species2}.temp2
cat $OutDir3/${species2}.temp2 $OutDir3/${species2}.temp1 > $OutDir3/${species2}.geno.txt && rm $OutDir3/${species2}.non_ref.freq $OutDir3/${species2}.temp1 $OutDir3/${species2}.temp2

grep -v "##" $OutDir2/${hybrids}.gt.vcf |cut -f 10- > $OutDir3/${hybrids}.temp
grep -v "##" $OutDir2/${hybrids}.gt.vcf |cut -f 1,2 > $OutDir3/${hybrids}.temp2

paste $OutDir3/${hybrids}.temp2 $OutDir3/${species1}.geno.txt $OutDir3/${species2}.geno.txt $OutDir3/${hybrids}.temp > $OutDir3/${species1}.${species2}.${hybrids}.geno.txt && rm $OutDir3/${hybrids}.temp2 $OutDir3/${species1}.geno.txt $OutDir3/${species2}.geno.txt $OutDir3/${hybrids}.temp

sed 's/#//g' $OutDir3/${species1}.${species2}.${hybrids}.geno.txt > $OutDir3/temp && mv $OutDir3/temp $OutDir3/${species1}.${species2}.${hybrids}.geno.txt 

elif [ $step == "4" ]; then
#################################################
###running the python script to calculate the hybrid index and heterozygosity for the potential hybrids
python3.6 hybrid_index.ancestral_proportion.py $OutDir3 ${species1}.${species2}.${hybrids}.geno.txt

elif [ $step == "5" ]; then
#################################################
###remove other unused files
rm $OutDir2/*



fi


