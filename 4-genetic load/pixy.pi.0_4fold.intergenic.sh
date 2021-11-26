#!/bin/sh
#SBATCH -c 1 --mem 15G
#SBATCH -o /UserHome/user260/pipeline/population_genetics/pi_fst_dxy/populus227.pixy.out
#SBATCH -e /UserHome/user260/pipeline/population_genetics/pi_fst_dxy/populus227.pixy.err

vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
pixy="/UserHome/user260/miniconda3/envs/pixy/bin/pixy"
mean="/UserHome/user260/tools/small_tools/mean"

step=$1
species=$2
chr=$3

InputDir="/w/user261/Jointvcf/filtering"
vcf=$InputDir/$chr.g.vcf.gz
#bed="/w/user261/data/snpable/bed/Potra02_genome.$chr.mask_mappability.bed"
bed="/w/user261/data/snpable/bed/Potra02_genome.$chr.mask_mappability.bed"
bed_0fold="/w/user260/multiple_aspen/gff/0_4_fold/Pte_degenerate_0fold.bed"
bed_4fold="/w/user260/multiple_aspen/gff/0_4_fold/Pte_degenerate_4fold.bed"


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

##individual output directory
OutDir_individual=$OutDir/pixy/individual/$species/pi
if [ ! -d "$OutDir_individual" ]; then
mkdir -p $OutDir_individual
fi


##pi directory
OutDir_pixy_out_pi=$OutDir/pixy/out/pi
if [ ! -d "$OutDir_pixy_out_pi" ]; then
mkdir -p $OutDir_pixy_out_pi
fi

##pi 0-fold diversity
OutDir_pixy_out_pi_0fold=$OutDir_individual/0fold
if [ ! -d "$OutDir_pixy_out_p_0fold" ]; then
mkdir -p $OutDir_pixy_out_pi_0fold
fi

##pi 4-fold diversity
OutDir_pixy_out_pi_4fold=$OutDir_individual/4fold
if [ ! -d "$OutDir_pixy_out_p_4fold" ]; then
mkdir -p $OutDir_pixy_out_pi_4fold
fi

##pi intergenic diversity
OutDir_pixy_out_pi_intergenic=$OutDir_individual/intergenic
if [ ! -d "$OutDir_pixy_out_p_intergenic" ]; then
mkdir -p $OutDir_pixy_out_pi_intergenic
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

#elif [ $step == "3" ]; then
##extract the 190 samples
cut -f 3 -d "," $populus190_samples |sed '1d' > $OutDir_pixy_samples/populus190.txt
$vcftools --gzvcf $OutDir_pixy/$OutSNPSuffix.filtering.vcf.gz --keep $OutDir_pixy_samples/populus190.txt --bed $bed --max-missing 0.8 --recode --recode-INFO-all --out $OutDir_pixy_filtering/populus190.$chr.filtering
$bcftools view $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf -O z -o $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf.gz && rm $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf


elif [ $step == "3" ]; then
#0fold
$vcftools --gzvcf $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf.gz --bed $bed_0fold --recode --recode-INFO-all --out $OutDir_pixy_filtering/0_fold/populus190.$chr.0_fold.filtering
$bcftools view $OutDir_pixy_filtering/0_fold/populus190.$chr.0_fold.filtering.recode.vcf -O z -o $OutDir_pixy_filtering/0_fold/populus190.$chr.0_fold.filtering.recode.vcf.gz && rm $OutDir_pixy_filtering/0_fold/populus190.$chr.0_fold.filtering.recode.vcf

elif [ $step == "3.3" ]; then
#4fold
$vcftools --gzvcf $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf.gz --bed $bed_4fold --recode --recode-INFO-all --out $OutDir_pixy_filtering/4_fold/populus190.$chr.4_fold.2.filtering
$bcftools view $OutDir_pixy_filtering/4_fold/populus190.$chr.4_fold.2.filtering.recode.vcf -O z -o $OutDir_pixy_filtering/4_fold/populus190.$chr.4_fold.2.filtering.recode.vcf.gz && rm $OutDir_pixy_filtering/4_fold/populus190.$chr.4_fold.2.filtering.recode.vcf


elif [ $step == "3.1" ]; then
#intergenic
mask_bed="/w/user260/multiple_aspen/snpable/Potra02_genome.mask_mappability.bed"
gene_gff="/w/user260/multiple_aspen/gff/Potra02_genes.gff"
intergenic_OutDir="/w/user260/multiple_aspen/gff/intergenic"

##removing scaffolds first
#grep -v "scaffold" $mask_bed > $intergenic_OutDir/mask.no_scaffold.bed

#extract genes
#mask_bed="$intergenic_OutDir/Potra02_genome.mask_mappability.bed"
#grep -v "#" $gene_gff | grep "gene" ||sort -k1,1V -k4,4n | $bedtools merge -i - > $intergenic_OutDir/Potra02.gene.bed
#extract 1kb regions flaning genes
#awk -F "\t" 'BEGIN{OFS="\t"}{print $1,$2-100,$3+100}' $intergenic_OutDir/Potra02.gene.bed  > $intergenic_OutDir/Potra02.gene_flank1kb.bed
#$bedtools slop -i $intergenic_OutDir/Potra02.gene.bed -g $mask_bed -b 100  > $intergenic_OutDir/Potra02.gene_flank1kb.bed
##extract intergenic region
#$bedtools subtract -a $intergenic_OutDir/mask.no_scaffold.bed -b $intergenic_OutDir/Potra02.gene.bed > $intergenic_OutDir/Potra02.intergenic.bed
##extract only intergenic with 1kb from genes
#$bedtools sort -i $intergenic_OutDir/Potra02.gene.flank1kb.bed
grep -v "#" $gene_gff | grep "gene" ||sort -k1,1V -k4,4n | $bedtools merge -i - > $intergenic_OutDir/Potra02.gene.bed
$bedtools subtract -a $mask_bed -b $intergenic_OutDir/Potra02.gene.bed > $intergenic_OutDir/Potra02.intergenic.bed


elif [ $step == "3.2" ]; then
#intergenic
intergenic_OutDir="/w/user260/multiple_aspen/gff/intergenic"
awk '$1=="'$chr'"' $intergenic_OutDir/Potra02.intergenic.bed > $intergenic_OutDir/Potra02.$chr.intergenic.bed

$vcftools --gzvcf $OutDir_pixy_filtering/populus190.$chr.filtering.recode.vcf.gz --bed $intergenic_OutDir/Potra02.$chr.intergenic.bed --recode --recode-INFO-all --out $OutDir_pixy_filtering/intergenic/populus190.$chr.intergenic.filtering
$bcftools view $OutDir_pixy_filtering/intergenic/populus190.$chr.intergenic.filtering.recode.vcf -O z -o $OutDir_pixy_filtering/intergenic/populus190.$chr.intergenic.filtering.recode.vcf.gz && rm $OutDir_pixy_filtering/intergenic/populus190.$chr.intergenic.filtering.recode.vcf


elif [ $step == "4.0" ]; then
##extracting individuals first
species_samples=$OutDir_pixy_samples/$species
if [ ! -d "$species_samples" ]; then
mkdir -p $species_samples
fi

n=$(cat $OutDir_pixy_samples/$species.txt |wc -l)
echo $n

for i in $(seq 1 $n)
do
echo $i
ind=$(head -n $i $OutDir_pixy_samples/$species.txt |tail -n 1 |cut -f 1)

head -n $i $OutDir_pixy_samples/$species.txt |tail -n 1 > $species_samples/${ind}.txt

done

elif [ $step == "4.1" ]; then
##species: Pade, Palb, Pdav, Pqio, Prot, Ptra, Ptrs
##estimating 4-fold sites pi

species_samples=$OutDir_pixy_samples/$species
if [ ! -d "$species_samples" ]; then
mkdir -p $species_samples
fi

n=$(cat $OutDir_pixy_samples/$species.txt |wc -l)
echo $n

for i in $(seq 1 $n)
do
echo $i
ind=$(head -n $i $OutDir_pixy_samples/$species.txt |tail -n 1 |cut -f 1)


OutDir_pixy_temp_pi=$OutDir/pixy/temp/pi/$species
if [ ! -d "$OutDir_pixy_temp_pi" ]; then
mkdir -p $OutDir_pixy_temp_pi
fi

##4_fold synonymous pi estimation
$pixy --stats pi \
	--vcf $OutDir_pixy_filtering/4_fold/populus190.$chr.4_fold.filtering.recode.vcf.gz \
	--zarr_path $OutDir_pixy_temp_pi \
	--window_size 100000 \
	--populations $species_samples/${ind}.txt \
	--bypass_filtration yes \
	--outfile_prefix $OutDir_pixy_out_pi_4fold/${ind}.$chr.4_fold.pixy.pi.out

done


elif [ $step == "4.2" ]; then
##species: Pade, Palb, Pdav, Pqio, Prot, Ptra, Ptrs
##estimating 0-fold pi

species_samples=$OutDir_pixy_samples/$species
if [ ! -d "$species_samples" ]; then
mkdir -p $species_samples
fi

n=$(cat $OutDir_pixy_samples/$species.txt |wc -l)
echo $n

for i in $(seq 1 $n)
do
echo $i
ind=$(head -n $i $OutDir_pixy_samples/$species.txt |tail -n 1 |cut -f 1)


OutDir_pixy_temp_pi=$OutDir/pixy/temp/pi/$species
if [ ! -d "$OutDir_pixy_temp_pi" ]; then
mkdir -p $OutDir_pixy_temp_pi
fi

##0_fold synonymous pi estimation
$pixy --stats pi \
        --vcf $OutDir_pixy_filtering/0_fold/populus190.$chr.0_fold.filtering.recode.vcf.gz \
        --zarr_path $OutDir_pixy_temp_pi \
        --window_size 100000 \
        --populations $species_samples/${ind}.txt \
        --bypass_filtration yes \
        --outfile_prefix $OutDir_pixy_out_pi_4fold/${ind}.$chr.0_fold.pixy.pi.out

done

elif [ $step == "4.3" ]; then
##species: Pade, Palb, Pdav, Pqio, Prot, Ptra, Ptrs
##estimating intergenic pi

species_samples=$OutDir_pixy_samples/$species
if [ ! -d "$species_samples" ]; then
mkdir -p $species_samples
fi

n=$(cat $OutDir_pixy_samples/$species.txt |wc -l)
echo $n

for i in $(seq 1 $n)
do
echo $i
ind=$(head -n $i $OutDir_pixy_samples/$species.txt |tail -n 1 |cut -f 1)


OutDir_pixy_temp_pi=$OutDir/pixy/temp/pi/$species
if [ ! -d "$OutDir_pixy_temp_pi" ]; then
mkdir -p $OutDir_pixy_temp_pi
fi

##intergenic synonymous pi estimation
$pixy --stats pi \
        --vcf $OutDir_pixy_filtering/intergenic/populus190.$chr.intergenic.filtering.recode.vcf.gz \
        --zarr_path $OutDir_pixy_temp_pi \
        --window_size 100000 \
        --populations $species_samples/${ind}.txt \
        --bypass_filtration yes \
        --outfile_prefix $OutDir_pixy_out_pi_4fold/${ind}.$chr.intergenic.pixy.pi.out

done


elif [ $step == "5" ]; then
species=$2

OutDir_pixy_out_pi_summary=$OutDir_individual/summary
if [ ! -d "$OutDir_pixy_out_pi_summary" ]; then
mkdir -p $OutDir_pixy_out_pi_summary
fi

species_samples=$OutDir_pixy_samples/$species
if [ ! -d "$species_samples" ]; then
mkdir -p $species_samples
fi

n=$(cat $OutDir_pixy_samples/$species.txt |wc -l)
echo $n

for i in $(seq 1 $n)
do
echo $i
ind=$(head -n $i $OutDir_pixy_samples/$species.txt |tail -n 1 |cut -f 1)


input=
for chr in chr{2..19}
do
sed '1d' $OutDir_pixy_out_pi_4fold/${ind}.$chr.intergenic.pixy.pi.out_pi.txt > $OutDir_pixy_out_pi_summary/${ind}.$chr.intergenic.pixy.pi.out.temp
input="$input $OutDir_pixy_out_pi_summary/${ind}.$chr.intergenic.pixy.pi.out.temp"
done

cat $OutDir_pixy_out_pi_4fold/${ind}.chr1.intergenic.pixy.pi.out_pi.txt $input > $OutDir_pixy_out_pi_summary/${ind}.intergenic.pixy.pi.txt
rm $input

done

elif [ $step == "5.1" ]; then
species=$2

OutDir_pixy_out_pi_summary=$OutDir_individual/summary
if [ ! -d "$OutDir_pixy_out_pi_summary" ]; then
mkdir -p $OutDir_pixy_out_pi_summary
fi

species_samples=$OutDir_pixy_samples/$species
if [ ! -d "$species_samples" ]; then
mkdir -p $species_samples
fi

n=$(cat $OutDir_pixy_samples/$species.txt |wc -l)
echo $n

for i in $(seq 1 $n)
do
echo $i
ind=$(head -n $i $OutDir_pixy_samples/$species.txt |tail -n 1 |cut -f 1)


input=
for chr in chr{2..19}
do
sed '1d' $OutDir_pixy_out_pi_4fold/${ind}.$chr.4_fold.pixy.pi.out_pi.txt > $OutDir_pixy_out_pi_summary/${ind}.$chr.4_fold.pixy.pi.out.temp
input="$input $OutDir_pixy_out_pi_summary/${ind}.$chr.4_fold.pixy.pi.out.temp"
done

cat $OutDir_pixy_out_pi_4fold/${ind}.chr1.4_fold.pixy.pi.out_pi.txt $input > $OutDir_pixy_out_pi_summary/${ind}.4_fold.pixy.pi.txt
rm $input

done

elif [ $step == "5.2" ]; then
species=$2

OutDir_pixy_out_pi_summary=$OutDir_individual/summary
if [ ! -d "$OutDir_pixy_out_pi_summary" ]; then
mkdir -p $OutDir_pixy_out_pi_summary
fi

species_samples=$OutDir_pixy_samples/$species
if [ ! -d "$species_samples" ]; then
mkdir -p $species_samples
fi

n=$(cat $OutDir_pixy_samples/$species.txt |wc -l)
echo $n

for i in $(seq 1 $n)
do
echo $i
ind=$(head -n $i $OutDir_pixy_samples/$species.txt |tail -n 1 |cut -f 1)


input=
for chr in chr{2..19}
do
sed '1d' $OutDir_pixy_out_pi_4fold/${ind}.$chr.0_fold.pixy.pi.out_pi.txt > $OutDir_pixy_out_pi_summary/${ind}.$chr.0_fold.pixy.pi.out.temp
input="$input $OutDir_pixy_out_pi_summary/${ind}.$chr.0_fold.pixy.pi.out.temp"
done

cat $OutDir_pixy_out_pi_4fold/${ind}.chr1.0_fold.pixy.pi.out_pi.txt $input > $OutDir_pixy_out_pi_summary/${ind}.0_fold.pixy.pi.txt
rm $input

done


elif [ $step == "6" ]; then
##summarize the results into tables with four columns: individual, theta_0_fold, theta_4_fold, theta_intergenic
species=$2

OutDir_pixy_out_pi_summary=$OutDir_individual/summary
if [ ! -d "$OutDir_pixy_out_pi_summary" ]; then
mkdir -p $OutDir_pixy_out_pi_summary
fi

species_samples=$OutDir_pixy_samples/$species
if [ ! -d "$species_samples" ]; then
mkdir -p $species_samples
fi

n=$(cat $OutDir_pixy_samples/$species.txt |wc -l)
echo $n

echo -e "ind\t0_fold_theta\t4_fold_theta\tintergenic_theta" > $OutDir_pixy_out_pi_summary/$species.theta.summary.txt

for i in $(seq 1 $n)
do
echo $i
ind=$(head -n $i $OutDir_pixy_samples/$species.txt |tail -n 1 |cut -f 1)
zero_mean=$($mean col=5 $OutDir_pixy_out_pi_summary/${ind}.0_fold.pixy.pi.txt)
four_mean=$($mean col=5 $OutDir_pixy_out_pi_summary/${ind}.4_fold.pixy.pi.txt)
intergenic_mean=$($mean col=5 $OutDir_pixy_out_pi_summary/${ind}.intergenic.pixy.pi.txt)

echo -e "$ind\t$zero_mean\t$four_mean\t$intergenic_mean" >> $OutDir_pixy_out_pi_summary/$species.theta.summary.txt

done


fi



