#!/bin/sh
#SBATCH -c 4 --mem 60G
#SBATCH -o /UserHome/user260/pipeline/population_genetics/pi_fst_dxy/populus227.pixy.out
#SBATCH -e /UserHome/user260/pipeline/population_genetics/pi_fst_dxy/populus227.pixy.err

vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
pixy="/UserHome/user260/miniconda3/envs/pixy/bin/pixy"

step=$1
chr=$2
species=$3

InputDir="/w/user261/Jointvcf/filtering"
vcf=$InputDir/$chr.g.vcf.gz
#bed="/w/user261/data/snpable/bed/Potra02_genome.$chr.mask_mappability.bed"
bed="/w/user261/data/snpable/bed/Potra02_genome.$chr.mask_mappability.bed"

OutDir="/w/user260/multiple_aspen/populus190_vcf"
populus190_samples="/w/user260/multiple_aspen/populus190_vcf/samples/populus190.samples.remove_clones.csv"


OutDir_pixy=$OutDir/pixy
if [ ! -d "$OutDir_pixy" ]; then
mkdir -p $OutDir_pixy
fi

OutDir_pixy_samples=$OutDir/pixy/samples
if [ ! -d "$OutDir_pixy_samples" ]; then
mkdir -p $OutDir_pixy_samples
fi

OutDir_pixy_filtering=$OutDir/pixy/filtering
if [ ! -d "$OutDir_pixy_filtering" ]; then
mkdir -p $OutDir_pixy_filtering
fi

OutDir_pixy_temp=$OutDir/pixy/temp
if [ ! -d "$OutDir_pixy_temp" ]; then
mkdir -p $OutDir_pixy_temp
fi

OutDir_pixy_out=$OutDir/pixy/out
if [ ! -d "$OutDir_pixy_out" ]; then
mkdir -p $OutDir_pixy_out
fi

##pi directory
OutDir_pixy_out_pi=$OutDir/pixy/out/pi/100kb
if [ ! -d "$OutDir_pixy_out_pi" ]; then
mkdir -p $OutDir_pixy_out_pi
fi

##fst directory
OutDir_pixy_out_fst=$OutDir/pixy/out/fst/100kb
if [ ! -d "$OutDir_pixy_out_fst" ]; then
mkdir -p $OutDir_pixy_out_fst
fi

##dxy directory
OutDir_pixy_out_dxy=$OutDir/pixy/out/dxy/100kb
if [ ! -d "$OutDir_pixy_out_dxy" ]; then
mkdir -p $OutDir_pixy_out_dxy
fi

Out=${vcf##*/}
OutSNPSuffix=${Out%.vcf.gz}.rm_indel
echo $OutSNPSuffix


if [ $step == "0" ]; then

for species in {Pade,Palb,Pdav,Pqio,Prot,Ptra,Ptrs}
do
awk -F "," '{if($2~/'$species'/ && $7~/Pure/)print $3,$2}' $populus190_samples |sed 's/ /\t/g'  > $OutDir_pixy_samples/$species.txt
done


elif [ $step == "1" ]; then
####The first script is to remove indels and multi-allele SNPs
$vcftools --gzvcf $vcf --remove-indels --max-alleles 2 --recode --recode-INFO-all --out $OutDir_pixy/$OutSNPSuffix && $pigz $OutDir_pixy/$OutSNPSuffix.recode.vcf

elif [ $step == "2" ]; then
####set the sites (both invariant and variant) that does not have high quality
$bcftools filter -i '(FMT/DP>=5 & FMT/RGQ>=20) | (FMT/DP>=5 & FMT/GQ>=20)' -S . -O z -o $OutDir_pixy/$OutSNPSuffix.filtering.vcf.gz $OutDir_pixy/$OutSNPSuffix.recode.vcf.gz 

##extract the 190 samples
cut -f 3 -d "," $populus190_samples |sed '1d' > $OutDir_pixy_samples/populus190.txt
$vcftools --gzvcf $OutDir_pixy/$OutSNPSuffix.filtering.vcf.gz --keep $OutDir_pixy_samples/populus190.txt --bed $bed --max-missing 0.8 --recode --recode-INFO-all --out $OutDir_pixy_filtering/populus190.$chr.filtering
$bcftools view $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf -O z -o $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf.gz && rm $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf

elif [ $step == "3" ]; then
##species: Pade, Palb, Pdav, Pqio, Prot, Ptra, Ptrs
##estimating pi
OutDir_pixy_temp_pi=$OutDir/pixy/temp/pi/100kb/$species
if [ ! -d "$OutDir_pixy_temp_pi" ]; then
mkdir -p $OutDir_pixy_temp_pi
fi

$pixy --stats pi \
	--vcf $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf.gz \
	--zarr_path $OutDir_pixy_temp_pi \
	--window_size 100000 \
	--populations $OutDir_pixy_samples/$species.txt \
	--bypass_filtration yes \
	--outfile_prefix $OutDir_pixy_out_pi/$species.$chr.pixy.pi.out

elif [ $step == "3.1" ]; then
##merging
species=$2

##pi summary_directory
OutDir_pixy_out_pi_summary=$OutDir/pixy/out/pi/100kb/summary
if [ ! -d "$OutDir_pixy_out_pi_summary" ]; then
mkdir -p $OutDir_pixy_out_pi_summary
fi

input=
for chr in chr{2..19}
do
sed '1d' $OutDir_pixy_out_pi/$species.$chr.pixy.pi.out_pi.txt > $OutDir_pixy_out_pi/$species.$chr.pixy.pi.out.temp
input="$input $OutDir_pixy_out_pi/$species.$chr.pixy.pi.out.temp"
done

cat $OutDir_pixy_out_pi/$species.chr1.pixy.pi.out_pi.txt $input > $OutDir_pixy_out_pi_summary/$species.pixy.pi.out_pi.txt
rm $input



elif [ $step == "4" ]; then
###estimating fst and dxy between pairs of species, in total of 21 pairs
species=(Pade Palb Pdav Pqio Prot Ptra Ptrs)

for ((i=0; i<=5; i++))
do
b=`expr $i + 1`
    for ((j=$b; j<=6; j++))
    do
    cat $OutDir_pixy_samples/${species[i]}.txt $OutDir_pixy_samples/${species[j]}.txt > $OutDir_pixy_samples/${species[i]}_${species[j]}.txt

    done
done

elif [ $step == "5" ]; then

##estimating fst between pairs of species
species1=$3
species2=$4

OutDir_pixy_temp_fst=$OutDir/pixy/temp/fst/100kb/${species1}_${species2}
if [ ! -d "$OutDir_pixy_temp_fst" ]; then
mkdir -p $OutDir_pixy_temp_fst
fi

$pixy --stats fst \
        --vcf $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf.gz \
        --zarr_path $OutDir_pixy_temp_fst \
        --window_size 100000 \
        --populations $OutDir_pixy_samples/${species1}_${species2}.txt \
        --bypass_filtration yes \
        --outfile_prefix $OutDir_pixy_out_fst/${species1}.${species2}.$chr.pixy.fst.out
###remove teh temperary folders
rm -r $OutDir_pixy_temp_fst/$chr

elif [ $step == "5.1" ]; then
##merging
species1=$2
species2=$3

##pi summary_directory
OutDir_pixy_out_fst_summary=$OutDir/pixy/out/fst/100kb/summary
if [ ! -d "$OutDir_pixy_out_fst_summary" ]; then
mkdir -p $OutDir_pixy_out_fst_summary
fi

input=
for chr in chr{2..19}
do
sed '1d' $OutDir_pixy_out_fst/$species1.$species2.$chr.pixy.fst.out_fst.txt > $OutDir_pixy_out_fst/$species1.$species2.$chr.pixy.fst.out.temp
input="$input $OutDir_pixy_out_fst/$species1.$species2.$chr.pixy.fst.out.temp"
done

cat $OutDir_pixy_out_fst/$species1.$species2.chr1.pixy.fst.out_fst.txt $input > $OutDir_pixy_out_fst_summary/$species1.$species2.pixy.fst.out_fst.txt
rm $input


elif [ $step == "6" ]; then

##estimating dxy between pairs of species
species1=$3
species2=$4

OutDir_pixy_temp_dxy=$OutDir/pixy/temp/dxy/100kb/${species1}_${species2}
if [ ! -d "$OutDir_pixy_temp_dxy" ]; then
mkdir -p $OutDir_pixy_temp_dxy
fi

$pixy --stats dxy \
        --vcf $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf.gz \
        --zarr_path $OutDir_pixy_temp_dxy \
        --window_size 100000 \
        --populations $OutDir_pixy_samples/${species1}_${species2}.txt \
        --bypass_filtration yes \
        --outfile_prefix $OutDir_pixy_out_dxy/${species1}.${species2}.$chr.pixy.dxy.out

rm -r $OutDir_pixy_temp_dxy/$chr

elif [ $step == "6.1" ]; then
##merging
species1=$2
species2=$3

##pi summary_directory
OutDir_pixy_out_dxy_summary=$OutDir/pixy/out/dxy/100kb/summary
if [ ! -d "$OutDir_pixy_out_dxy_summary" ]; then
mkdir -p $OutDir_pixy_out_dxy_summary
fi

input=
for chr in chr{2..19}
do
sed '1d' $OutDir_pixy_out_dxy/$species1.$species2.$chr.pixy.dxy.out_dxy.txt > $OutDir_pixy_out_dxy/$species1.$species2.$chr.pixy.dxy.out.temp
input="$input $OutDir_pixy_out_dxy/$species1.$species2.$chr.pixy.dxy.out.temp"
done

cat $OutDir_pixy_out_dxy/$species1.$species2.chr1.pixy.dxy.out_dxy.txt $input > $OutDir_pixy_out_dxy_summary/$species1.$species2.pixy.dxy.out_dxy.txt
rm $input


fi




