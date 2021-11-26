#!/bin/sh
#SBATCH -c 2 --mem 20G
#SBATCH -o /UserHome/user260/pipeline/hybridization/Dsuite/introgressed_regions/introgressed.dsuite.out
#SBATCH -e /UserHome/user260/pipeline/hybridization/Dsuite/introgressed_regions/introgressed.dsuite.err

#The main aim of this study:
#1. identify the introgressed regions for pairs of species from Dsuite results
#2. compare the genomic characteristic between the introgressed regions and the shuffled regions with the same size
#including (1) deleterious mutation (deleterious+loss function)/synonymous); also defer by heterozygous and homozygous
#          (2) diversity
#          (3) recombination rate
#          (4) divergence (fst and dxy)
#	   (5) coding density
#	



vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
small_tools="/UserHome/user260/tools/small_tools"
plink2="/UserHome/user260/tools/plink/plink2"

step=$1
pairs=$2  ##Ptra_Pdav_Pqio
top_window=$3 ##206

#species_pairs window_size
#1.Ptra_Pdav_Pqio 206
#2.Pade_Pqio_Ptra 133
#3.Pade_Pqio_Ptrs 165
#4.Ptra_Palb_Pqio 598
#5.Ptra_Prot_Pqio 273
#6.Ptrs_Palb_Pade 582
#7.Ptrs_Ptra_Pade 93
#8.Ptrs_Pdav_Pade 459
#9.Ptrs_Prot_Pade 503
#10.Ptrs_Ptra_Palb 1141
#11.Ptrs_Pdav_Palb 1314
#12.Ptrs_Prot_Palb 1288
#13.Prot_Ptra_Ptrs 3108
#14.Prot_Pdav_Ptrs 1133
#15.Prot_Pdav_Ptra 3237


P1="$(echo $pairs | cut -d '_' -f 1)"
P2="$(echo $pairs | cut -d '_' -f 2)"
P3="$(echo $pairs | cut -d '_' -f 3)"

echo $P1
echo $P2
echo $P3

###sift deleterious derived allele frequency
sift_input="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/derived_allele_frequency"

###dsuite input
dsuite_input="/w/user260/multiple_aspen/populus190_vcf/dsuite/dsuite/dinvestigate/out/50_20"

#diversity input
pi_input="/w/user260/multiple_aspen/populus190_vcf/pixy/out/pi/summary/bed"

#recombination input
ldhat_input="/w/user260/multiple_aspen/populus190_vcf/ldhat/summary/bed"

#dxy
dxy_input="/w/user260/multiple_aspen/populus190_vcf/pixy/out/dxy/summary/bed"

#fst
fst_input="/w/user260/multiple_aspen/populus190_vcf/pixy/out/fst/summary/bed"

#volcan
volcan_input="/w/user260/multiple_aspen/populus162_vcf/vcf/volcanofinder/summary"

#selscan
selscan_input="/w/user260/multiple_aspen/populus190_vcf/beagle/selscan/summary"


OutDir_introgressed=$dsuite_input/introgressed_islands/$pairs
if [ ! -d "$OutDir_introgressed" ]; then
mkdir -p $OutDir_introgressed
fi

OutDir_permutation=$dsuite_input/introgressed_islands/$pairs/permutation
if [ ! -d "$OutDir_permutation" ]; then
mkdir -p $OutDir_permutation
fi

OutDir_introgressed_summary=$dsuite_input/introgressed_islands/$pairs/summary
if [ ! -d "$OutDir_introgressed_summary" ]; then
mkdir -p $OutDir_introgressed_summary
fi


OutDir_introgressed_adaptive=$dsuite_input/introgressed_islands/$pairs/adaptive_comparison
if [ ! -d "$OutDir_introgressed_adaptive" ]; then
mkdir -p $OutDir_introgressed_adaptive
fi



potra_genome="/w/user260/multiple_aspen/snpable/Potra02.genome"
cds_bed="/w/user260/multiple_aspen/gff/Potra02_genes.CDS.bed"

if [ $step == "0" ]; then
#step 0: create the merged introgression islands
#total genome size: 361795894
head -n1 $dsuite_input/${pairs}_localFstats__50_20.txt > $OutDir_introgressed/$pairs.introgressed.50_20.windows.txt
sort -k6,6nr $dsuite_input/${pairs}_localFstats__50_20.txt |head -n $top_window |sed '1d'| awk '{$2=$2-1;print}' |sed 's/ /\t/g' >> $OutDir_introgressed/$pairs.introgressed.50_20.windows.txt

sort -k1,1V -k2,2n $OutDir_introgressed/$pairs.introgressed.50_20.windows.txt |cut -f 1-3 |sed '1d' > $OutDir_introgressed/$pairs.introgressed.50_20.windows.bed
$bedtools merge -i $OutDir_introgressed/$pairs.introgressed.50_20.windows.bed -d 1000 > $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed 


elif [ $step == "1" ]; then
#step1: create bed file for CDS regions and compare the CDS density within introgressed regions and other regions
#grep -v "#" Potra02_genes.gff |grep "CDS" |cut -f 1,4,5 |awk '{$2=$2-1;print}' |sed 's/ /\t/g' |grep -v "scaffold" > Potra02_genes.CDS.bed

$bedtools coverage -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $cds_bed > $OutDir_introgressed/$pairs.introgressed_islands.cds_density.txt

##summarize
echo -e "data\tcds_density" > $OutDir_introgressed/$pairs.cds_density.summary.txt
cds_den=$(cat $OutDir_introgressed/$pairs.introgressed_islands.cds_density.txt |$small_tools/mean col=7)

echo -e "real\t$cds_den" >> $OutDir_introgressed/$pairs.cds_density.summary.txt

#permutation
for i in {1..1000}
do

$bedtools coverage -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $cds_bed > $OutDir_permutation/$pairs.introgressed_islands.cds_density.permu.$i.txt
permu_cds_den=$(cat $OutDir_permutation/$pairs.introgressed_islands.cds_density.permu.$i.txt |$small_tools/mean col=7)

echo -e "permutation_$i\t$permu_cds_den" >> $OutDir_introgressed/$pairs.cds_density.summary.txt
done

elif [ $step == "2" ]; then
#step2: create bed file for the sift output file awk -F '\t' '{$2=$2-1"\t"$2}1' OFS='\t' $file |sed '1d' |cut -f 1-3
#P2
#1.deleteirous
$bedtools intersect -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $sift_input/deleterious/bed/$P2.deleterious.derived_frq.txt.bed -c > $OutDir_introgressed/$pairs.introgressed_islands.$P2.deleterious.txt
#2.loss_of_funtion
$bedtools intersect -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $sift_input/loss_function/bed/$P2.loss_of_function.derived_frq.txt.bed -c > $OutDir_introgressed/$pairs.introgressed_islands.$P2.loss_of_function.txt
#3.synonymous
$bedtools intersect -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $sift_input/synonymous/bed/$P2.synonymous.derived_frq.txt.bed -c > $OutDir_introgressed/$pairs.introgressed_islands.$P2.synonymous.txt
#4.tolerated
$bedtools intersect -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $sift_input/tolerated/bed/$P2.tolerated.derived_frq.txt.bed -c > $OutDir_introgressed/$pairs.introgressed_islands.$P2.tolerated.txt

##summarize
echo -e "data\tdeleterious\tloss_of_function\tsynonymous\ttolerated" > $OutDir_introgressed/$pairs.$P2.sift.summary.txt
p2_deleterious=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.deleterious.txt |$small_tools/sum col=4)
p2_loss_of_function=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.loss_of_function.txt |$small_tools/sum col=4)
p2_synonymous=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.synonymous.txt |$small_tools/sum col=4)
p2_tolerated=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.tolerated.txt |$small_tools/sum col=4)
 
echo -e "real\t$p2_deleterious\t$p2_loss_of_function\t$p2_synonymous\t$p2_tolerated" >> $OutDir_introgressed/$pairs.$P2.sift.summary.txt

##create 1000 permutation of the introgressed windows and check the enrichment of the deleterious mutations
for i in {1..1000}
do
$bedtools shuffle -i $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -g $potra_genome  > $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed

#P2
#1.deleteirous
$bedtools intersect -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $sift_input/deleterious/bed/$P2.deleterious.derived_frq.txt.bed -c > $OutDir_permutation/$pairs.introgressed_islands.$P2.deleterious.permu.$i.txt
#2.loss_of_funtion
$bedtools intersect -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $sift_input/loss_function/bed/$P2.loss_of_function.derived_frq.txt.bed -c > $OutDir_permutation/$pairs.introgressed_islands.$P2.loss_of_function.permu.$i.txt
#3.synonymous
$bedtools intersect -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $sift_input/synonymous/bed/$P2.synonymous.derived_frq.txt.bed -c > $OutDir_permutation/$pairs.introgressed_islands.$P2.synonymous.permu.$i.txt
#4.tolerated
$bedtools intersect -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $sift_input/tolerated/bed/$P2.tolerated.derived_frq.txt.bed -c > $OutDir_permutation/$pairs.introgressed_islands.$P2.tolerated.permu.$i.txt

#summarize
permu_p2_deleterious=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.deleterious.permu.$i.txt |$small_tools/sum col=4)
permu_p2_loss_of_function=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.loss_of_function.permu.$i.txt |$small_tools/sum col=4)
permu_p2_synonymous=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.synonymous.permu.$i.txt |$small_tools/sum col=4)
permu_p2_tolerated=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.tolerated.permu.$i.txt |$small_tools/sum col=4)

echo -e "permutation_$i\t$permu_p2_deleterious\t$permu_p2_loss_of_function\t$permu_p2_synonymous\t$permu_p2_tolerated" >> $OutDir_introgressed/$pairs.$P2.sift.summary.txt

done


#P3
#1.deleteirous
$bedtools intersect -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $sift_input/deleterious/bed/$P3.deleterious.derived_frq.txt.bed -c > $OutDir_introgressed/$pairs.introgressed_islands.$P3.deleterious.txt
#2.loss_of_funtion
$bedtools intersect -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $sift_input/loss_function/bed/$P3.loss_of_function.derived_frq.txt.bed -c > $OutDir_introgressed/$pairs.introgressed_islands.$P3.loss_of_function.txt
#3.synonymous
$bedtools intersect -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $sift_input/synonymous/bed/$P3.synonymous.derived_frq.txt.bed -c > $OutDir_introgressed/$pairs.introgressed_islands.$P3.synonymous.txt
#4.tolerated
$bedtools intersect -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $sift_input/tolerated/bed/$P3.tolerated.derived_frq.txt.bed -c > $OutDir_introgressed/$pairs.introgressed_islands.$P3.tolerated.txt
#summarize
echo -e "data\tdeleterious\tloss_of_function\tsynonymous\ttolerated" > $OutDir_introgressed/$pairs.$P3.sift.summary.txt
p3_deleterious=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.deleterious.txt |$small_tools/sum col=4)
p3_loss_of_function=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.loss_of_function.txt |$small_tools/sum col=4)
p3_synonymous=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.synonymous.txt |$small_tools/sum col=4)
p3_tolerated=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.tolerated.txt |$small_tools/sum col=4)

echo -e "real\t$p3_deleterious\t$p3_loss_of_function\t$p3_synonymous\t$p3_tolerated" >> $OutDir_introgressed/$pairs.$P3.sift.summary.txt

##create 1000 permutation of the introgressed windows and check the enrichment of the deleterious mutations
for i in {1..1000}
do

#P3
#1.deleteirous
$bedtools intersect -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $sift_input/deleterious/bed/$P3.deleterious.derived_frq.txt.bed -c > $OutDir_permutation/$pairs.introgressed_islands.$P3.deleterious.permu.$i.txt
#2.loss_of_funtion
$bedtools intersect -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $sift_input/loss_function/bed/$P3.loss_of_function.derived_frq.txt.bed -c > $OutDir_permutation/$pairs.introgressed_islands.$P3.loss_of_function.permu.$i.txt
#3.synonymous
$bedtools intersect -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $sift_input/synonymous/bed/$P3.synonymous.derived_frq.txt.bed -c > $OutDir_permutation/$pairs.introgressed_islands.$P3.synonymous.permu.$i.txt
#4.tolerated
$bedtools intersect -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $sift_input/tolerated/bed/$P3.tolerated.derived_frq.txt.bed -c > $OutDir_permutation/$pairs.introgressed_islands.$P3.tolerated.permu.$i.txt

#summarize
permu_p3_deleterious=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.deleterious.permu.$i.txt |$small_tools/sum col=4)
permu_p3_loss_of_function=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.loss_of_function.permu.$i.txt |$small_tools/sum col=4)
permu_p3_synonymous=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.synonymous.permu.$i.txt |$small_tools/sum col=4)
permu_p3_tolerated=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.tolerated.permu.$i.txt |$small_tools/sum col=4)

echo -e "permutation_$i\t$permu_p3_deleterious\t$permu_p3_loss_of_function\t$permu_p3_synonymous\t$permu_p3_tolerated" >> $OutDir_introgressed/$pairs.$P3.sift.summary.txt

done


elif [ $step == "3" ]; then
#compare diversity
#create the pi bed file for each species
#cut -f 2-5 $file |sed '1d' |awk '{$2=$2-1}1' |sed 's/ /\t/g'|grep -v "NA" > bed/$file.bed
#P2
sed 's/chr//g' $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed | sort -k1,1n -k 2,2n |awk '{$1="chr"$1}1' |sed 's/ /\t/g' > $OutDir_introgressed/temp && mv $OutDir_introgressed/temp $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $pi_input/$P2.pixy.pi.out_pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.pi.txt

##summarize
echo -e "data\taverage_pi" > $OutDir_introgressed/$pairs.$P2.pi.summary.txt
p2_pi=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.pi.txt |$small_tools/mean col=4)

echo -e "real\t$p2_pi" >> $OutDir_introgressed/$pairs.$P2.pi.summary.txt

#permutation
for i in {1..1000}
do

sed 's/chr//g' $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed |sort -k1,1n -k 2,2n |awk '{$1="chr"$1}1' |sed 's/ /\t/g' > $OutDir_permutation/temp && mv $OutDir_permutation/temp $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed


$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $pi_input/$P2.pixy.pi.out_pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.pi.permu.$i.txt

permu_p2_pi=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.pi.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_pi" >> $OutDir_introgressed/$pairs.$P2.pi.summary.txt 
done

#P3
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $pi_input/$P3.pixy.pi.out_pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.pi.txt

##summarize
echo -e "data\taverage_pi" > $OutDir_introgressed/$pairs.$P3.pi.summary.txt
p3_pi=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.pi.txt |$small_tools/mean col=4)

echo -e "real\t$p3_pi" >> $OutDir_introgressed/$pairs.$P3.pi.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $pi_input/$P3.pixy.pi.out_pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.pi.permu.$i.txt

permu_p3_pi=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.pi.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p3_pi" >> $OutDir_introgressed/$pairs.$P3.pi.summary.txt
done


elif [ $step == "3.1" ]; then
#compare 4-fold diversity
#create the pi bed file for each species
#cut -f 2-5 $file |sed '1d' |awk '{$2=$2-1}1' |sed 's/ /\t/g'|grep -v "NA" > bed/$file.bed

#P2
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $pi_input/$P2.4_fold.pixy.pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.pi.4fold.txt

##summarize
echo -e "data\taverage_pi_4fold" > $OutDir_introgressed/$pairs.$P2.pi.4fold.summary.txt
p2_pi_4fold=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.pi.4fold.txt |$small_tools/mean col=4)

echo -e "real\t$p2_pi_4fold" >> $OutDir_introgressed/$pairs.$P2.pi.4fold.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $pi_input/$P2.4_fold.pixy.pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.pi.4fold.permu.$i.txt

permu_p2_pi_4fold=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.pi.4fold.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_pi_4fold" >> $OutDir_introgressed/$pairs.$P2.pi.4fold.summary.txt
done


#P3
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $pi_input/$P3.4_fold.pixy.pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.pi.4fold.txt

##summarize
echo -e "data\taverage_pi_4fold" > $OutDir_introgressed/$pairs.$P3.pi.4fold.summary.txt
p3_pi_4fold=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.pi.4fold.txt |$small_tools/mean col=4)

echo -e "real\t$p3_pi_4fold" >> $OutDir_introgressed/$pairs.$P3.pi.4fold.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $pi_input/$P3.4_fold.pixy.pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.pi.4fold.permu.$i.txt

permu_p3_pi_4fold=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.pi.4fold.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p3_pi_4fold" >> $OutDir_introgressed/$pairs.$P3.pi.4fold.summary.txt
done


elif [ $step == "3.2" ]; then

#compare 0-fold diversity
#create the pi bed file for each species
#cut -f 2-5 $file |sed '1d' |awk '{$2=$2-1}1' |sed 's/ /\t/g'|grep -v "NA" > bed/$file.bed

#P2
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $pi_input/$P2.0_fold.pixy.pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.pi.0fold.txt

##summarize
echo -e "data\taverage_pi_0fold" > $OutDir_introgressed/$pairs.$P2.pi.0fold.summary.txt
p2_pi_0fold=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.pi.0fold.txt |$small_tools/mean col=4)

echo -e "real\t$p2_pi_0fold" >> $OutDir_introgressed/$pairs.$P2.pi.0fold.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $pi_input/$P2.0_fold.pixy.pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.pi.0fold.permu.$i.txt

permu_p2_pi_0fold=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.pi.0fold.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_pi_0fold" >> $OutDir_introgressed/$pairs.$P2.pi.0fold.summary.txt
done


#P3
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $pi_input/$P3.0_fold.pixy.pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.pi.0fold.txt

##summarize
echo -e "data\taverage_pi_0fold" > $OutDir_introgressed/$pairs.$P3.pi.0fold.summary.txt
p3_pi_0fold=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.pi.0fold.txt |$small_tools/mean col=4)

echo -e "real\t$p3_pi_0fold" >> $OutDir_introgressed/$pairs.$P3.pi.0fold.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $pi_input/$P3.0_fold.pixy.pi.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.pi.0fold.permu.$i.txt

permu_p3_pi_0fold=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.pi.0fold.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p3_pi_0fold" >> $OutDir_introgressed/$pairs.$P3.pi.0fold.summary.txt
done


elif [ $step == "4" ]; then
##recombination
#awk '{$2=$2-1" "$2}1' Pade.ldhat.summary.txt |sed '1d' |sed 's/ /\t/g' |cut -f 1-4

#P2

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $ldhat_input/$P2.ldhat.summary.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.ldhat.txt

##summarize
echo -e "data\taverage_ldhat" > $OutDir_introgressed/$pairs.$P2.ldhat.summary.txt
p2_ldhat=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.ldhat.txt |$small_tools/mean col=4)

echo -e "real\t$p2_ldhat" >> $OutDir_introgressed/$pairs.$P2.ldhat.summary.txt

#permutation
for i in {1..1000}
do


$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $ldhat_input/$P2.ldhat.summary.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.ldhat.permu.$i.txt

permu_p2_ldhat=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.ldhat.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_ldhat" >> $OutDir_introgressed/$pairs.$P2.ldhat.summary.txt
done

#P3
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $ldhat_input/$P3.ldhat.summary.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.ldhat.txt

##summarize
echo -e "data\taverage_ldhat" > $OutDir_introgressed/$pairs.$P3.ldhat.summary.txt
p3_ldhat=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.ldhat.txt |$small_tools/mean col=4)

echo -e "real\t$p3_ldhat" >> $OutDir_introgressed/$pairs.$P3.ldhat.summary.txt

#permutation
for i in {1..1000}
do
$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $ldhat_input/$P3.ldhat.summary.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.ldhat.permu.$i.txt

permu_p3_ldhat=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.ldhat.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p3_ldhat" >> $OutDir_introgressed/$pairs.$P3.ldhat.summary.txt
done



elif [ $step == "5" ]; then
#divergence fst

#determine the order P2-P3 or P3-P2 

if [ -f $fst_input/$P2.$P3.pixy.fst.out_fst.txt.bed ];then

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $fst_input/$P2.$P3.pixy.fst.out_fst.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.$P3.fst.txt

##summarize
echo -e "data\taverage_fst" > $OutDir_introgressed/$pairs.$P2.$P3.fst.summary.txt
p2_fst=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.$P3.fst.txt |$small_tools/mean col=4)

echo -e "real\t$p2_fst" >> $OutDir_introgressed/$pairs.$P2.$P3.fst.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $fst_input/$P2.$P3.pixy.fst.out_fst.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.$P3.fst.permu.$i.txt

permu_p2_fst=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.$P3.fst.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_fst" >> $OutDir_introgressed/$pairs.$P2.$P3.fst.summary.txt
done

else

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $fst_input/$P3.$P2.pixy.fst.out_fst.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.$P2.fst.txt

##summarize
echo -e "data\taverage_fst" > $OutDir_introgressed/$pairs.$P3.$P2.fst.summary.txt
p2_fst=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.$P2.fst.txt |$small_tools/mean col=4)

echo -e "real\t$p2_fst" >> $OutDir_introgressed/$pairs.$P3.$P2.fst.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $fst_input/$P3.$P2.pixy.fst.out_fst.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.$P2.fst.permu.$i.txt

permu_p2_fst=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.$P2.fst.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_fst" >> $OutDir_introgressed/$pairs.$P3.$P2.fst.summary.txt
done

fi


elif [ $step == "6" ]; then
#divergence dxy

#determine the order P2-P3 or P3-P2

if [ -f $dxy_input/$P2.$P3.pixy.dxy.out_dxy.txt.bed ];then

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $dxy_input/$P2.$P3.pixy.dxy.out_dxy.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.$P3.dxy.txt

##summarize
echo -e "data\taverage_dxy" > $OutDir_introgressed/$pairs.$P2.$P3.dxy.summary.txt
p2_dxy=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.$P3.dxy.txt |$small_tools/mean col=4)

echo -e "real\t$p2_dxy" >> $OutDir_introgressed/$pairs.$P2.$P3.dxy.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $dxy_input/$P2.$P3.pixy.dxy.out_dxy.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.$P3.dxy.permu.$i.txt

permu_p2_dxy=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.$P3.dxy.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_dxy" >> $OutDir_introgressed/$pairs.$P2.$P3.dxy.summary.txt
done


else

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $dxy_input/$P3.$P2.pixy.dxy.out_dxy.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.$P2.dxy.txt

##summarize
echo -e "data\taverage_dxy" > $OutDir_introgressed/$pairs.$P3.$P2.dxy.summary.txt
p2_dxy=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.$P2.dxy.txt |$small_tools/mean col=4)

echo -e "real\t$p2_dxy" >> $OutDir_introgressed/$pairs.$P3.$P2.dxy.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $dxy_input/$P3.$P2.pixy.dxy.out_dxy.txt.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.$P2.dxy.permu.$i.txt

permu_p2_dxy=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.$P2.dxy.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_dxy" >> $OutDir_introgressed/$pairs.$P3.$P2.dxy.summary.txt
done

fi


elif [ $step == "7" ]; then
#introgression sweep
#P2

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $volcan_input/$P2.volcan.out.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.volcan.txt

##summarize
echo -e "data\taverage_volcan" > $OutDir_introgressed/$pairs.$P2.volcan.summary.txt
p2_volcan=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.volcan.txt |$small_tools/mean col=4)

echo -e "real\t$p2_volcan" >> $OutDir_introgressed/$pairs.$P2.volcan.summary.txt

#permutation
for i in {1..1000}
do


$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $volcan_input/$P2.volcan.out.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.volcan.permu.$i.txt

permu_p2_volcan=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.volcan.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_volcan" >> $OutDir_introgressed/$pairs.$P2.volcan.summary.txt
done


#P3
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $volcan_input/$P3.volcan.out.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.volcan.txt

echo -e "data\taverage_volcan" > $OutDir_introgressed/$pairs.$P3.volcan.summary.txt
p3_volcan=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.volcan.txt |$small_tools/mean col=4)

echo -e "real\t$p3_volcan" >> $OutDir_introgressed/$pairs.$P3.volcan.summary.txt

#permutation
for i in {1..1000}
do
$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $volcan_input/$P3.volcan.out.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.volcan.permu.$i.txt

permu_p3_volcan=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.volcan.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p3_volcan" >> $OutDir_introgressed/$pairs.$P3.volcan.summary.txt
done


elif [ $step == "8" ]; then
#examine the enrichment of adaptive introgression in the introgressed regions

species=$2
##choose the top 5% (1806) adaptive introgressed regions for each speceis

sort -k4,4nr $volcan_input/$species.volcan.out.bed |head -n 1806 |sort -k1,1V -k2,2n > $volcan_input/$species.volcan.out.top5percentile.bed


elif [ $step == "8.1" ]; then

#introgression sweep
#P2

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $volcan_input/$P2.volcan.out.top5percentile.bed -c 4 -o count -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.volcan.count.txt

##summarize
echo -e "data\tcount_volcan" > $OutDir_introgressed/$pairs.$P2.volcan.count.summary.txt
p2_volcan=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.volcan.count.txt |$small_tools/sum col=4)

echo -e "real\t$p2_volcan" >> $OutDir_introgressed/$pairs.$P2.volcan.count.summary.txt

#permutation
for i in {1..1000}
do


$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $volcan_input/$P2.volcan.out.top5percentile.bed -c 4 -o count -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.volcan.count.permu.$i.txt

permu_p2_volcan=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.volcan.count.permu.$i.txt |$small_tools/sum col=4)

echo -e "permutation_$i\t$permu_p2_volcan" >> $OutDir_introgressed/$pairs.$P2.volcan.count.summary.txt
done

#P3
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $volcan_input/$P3.volcan.out.top5percentile.bed -c 4 -o count -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.volcan.count.txt

echo -e "data\tcount_volcan" > $OutDir_introgressed/$pairs.$P3.volcan.count.summary.txt
p3_volcan=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.volcan.count.txt |$small_tools/sum col=4)

echo -e "real\t$p3_volcan" >> $OutDir_introgressed/$pairs.$P3.volcan.count.summary.txt

#permutation
for i in {1..1000}
do
$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $volcan_input/$P3.volcan.out.top5percentile.bed -c 4 -o count -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.volcan.count.permu.$i.txt

permu_p3_volcan=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.volcan.count.permu.$i.txt |$small_tools/sum col=4)

echo -e "permutation_$i\t$permu_p3_volcan" >> $OutDir_introgressed/$pairs.$P3.volcan.count.summary.txt
done


elif [ $step == "9" ]; then
#P2
deleterious_p2_vcf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift_${P2}/est_sfs/derived_hom_het/$P2.deleterious.derived.recode.vcf.gz"
loss_of_function_p2_vcf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift_${P2}/est_sfs/derived_hom_het/$P2.loss_of_function.derived.recode.vcf.gz"
synonymous_p2_vcf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift_${P2}/est_sfs/derived_hom_het/$P2.synonymous.derived.recode.vcf.gz"

#deleterious
$bcftools index $deleterious_p2_vcf
$bcftools view -R $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed $deleterious_p2_vcf -O z -o $OutDir_introgressed/$P2.introgressed.deleterious.derived.recode.vcf.gz
$plink2 --vcf $OutDir_introgressed/$P2.introgressed.deleterious.derived.recode.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_introgressed/$P2.introgressed.deleterious.derived.hom_het

#loss-of-function
$bcftools index $loss_of_function_p2_vcf
$bcftools view -R $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed $loss_of_function_p2_vcf -O z -o $OutDir_introgressed/$P2.introgressed.loss_of_function.derived.recode.vcf.gz
$plink2 --vcf $OutDir_introgressed/$P2.introgressed.loss_of_function.derived.recode.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_introgressed/$P2.introgressed.loss_of_function.derived.hom_het

#synonymous
$bcftools index $synonymous_p2_vcf
$bcftools view -R $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed $synonymous_p2_vcf -O z -o $OutDir_introgressed/$P2.introgressed.synonymous.derived.recode.vcf.gz
$plink2 --vcf $OutDir_introgressed/$P2.introgressed.synonymous.derived.recode.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_introgressed/$P2.introgressed.synonymous.derived.hom_het

##summarize
echo -e "data\tdele_hom\tdele_het\tlof_hom\tlof_het\tsyn_hom\tsyn_het" > $OutDir_introgressed/$pairs.$P2.hom_het.summary.txt
p2_dele_hom=$(sed '1d' $OutDir_introgressed/$P2.introgressed.deleterious.derived.hom_het.scount |head -n -2 |$small_tools/sum col=3)
p2_dele_het=$(sed '1d' $OutDir_introgressed/$P2.introgressed.deleterious.derived.hom_het.scount |head -n -2 |$small_tools/sum col=4)
p2_lof_hom=$(sed '1d' $OutDir_introgressed/$P2.introgressed.loss_of_function.derived.hom_het.scount |head -n -2 |$small_tools/sum col=3)
p2_lof_het=$(sed '1d' $OutDir_introgressed/$P2.introgressed.loss_of_function.derived.hom_het.scount |head -n -2 |$small_tools/sum col=4)
p2_syno_hom=$(sed '1d' $OutDir_introgressed/$P2.introgressed.synonymous.derived.hom_het.scount |head -n -2 |$small_tools/sum col=3)
p2_syno_het=$(sed '1d' $OutDir_introgressed/$P2.introgressed.synonymous.derived.hom_het.scount |head -n -2 |$small_tools/sum col=4)

echo -e "real\t$p2_dele_hom\t$p2_dele_het\t$p2_lof_hom\t$p2_lof_het\t$p2_syno_hom\t$p2_syno_het" >> $OutDir_introgressed/$pairs.$P2.hom_het.summary.txt

#permutation
for i in {1..1000}
do

$bcftools view -R $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed $deleterious_p2_vcf -O z -o $OutDir_permutation/$P2.deleterious.hom_het.permu.$i.vcf.gz
$plink2 --vcf $OutDir_permutation/$P2.deleterious.hom_het.permu.$i.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_permutation/$P2.deleterious.hom_het.permu.$i

$bcftools view -R $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed $loss_of_function_p2_vcf -O z -o $OutDir_permutation/$P2.loss_of_function.hom_het.permu.$i.vcf.gz
$plink2 --vcf $OutDir_permutation/$P2.loss_of_function.hom_het.permu.$i.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_permutation/$P2.loss_of_function.hom_het.permu.$i

$bcftools view -R $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed $synonymous_p2_vcf -O z -o $OutDir_permutation/$P2.synonymous.hom_het.permu.$i.vcf.gz
$plink2 --vcf $OutDir_permutation/$P2.synonymous.hom_het.permu.$i.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_permutation/$P2.synonymous.hom_het.permu.$i

permu_p2_dele_hom=$(sed '1d' $OutDir_permutation/$P2.deleterious.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=3)
permu_p2_dele_het=$(sed '1d' $OutDir_permutation/$P2.deleterious.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=4)

permu_p2_lof_hom=$(sed '1d' $OutDir_permutation/$P2.loss_of_function.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=3)
permu_p2_lof_het=$(sed '1d' $OutDir_permutation/$P2.loss_of_function.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=4)

permu_p2_syno_hom=$(sed '1d' $OutDir_permutation/$P2.synonymous.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=3)
permu_p2_syno_het=$(sed '1d' $OutDir_permutation/$P2.synonymous.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=4)


echo -e "permutation_$i\t$permu_p2_dele_hom\t$permu_p2_dele_het\t$permu_p2_lof_hom\t$permu_p2_lof_het\t$permu_p2_syno_hom\t$permu_p2_syno_het" >> $OutDir_introgressed/$pairs.$P2.hom_het.summary.txt

done


#P3
deleterious_p3_vcf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift_${P3}/est_sfs/derived_hom_het/$P3.deleterious.derived.recode.vcf.gz"
loss_of_function_p3_vcf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift_${P3}/est_sfs/derived_hom_het/$P3.loss_of_function.derived.recode.vcf.gz"
synonymous_p3_vcf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift_${P3}/est_sfs/derived_hom_het/$P3.synonymous.derived.recode.vcf.gz"

$bcftools index $deleterious_p3_vcf
#deleteroius
$bcftools view -R $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed $deleterious_p3_vcf -O z -o $OutDir_introgressed/$P3.introgressed.deleterious.derived.recode.vcf.gz
$plink2 --vcf $OutDir_introgressed/$P3.introgressed.deleterious.derived.recode.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_introgressed/$P3.introgressed.deleterious.derived.hom_het

#loss-of-function
$bcftools index $loss_of_function_p3_vcf
$bcftools view -R $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed $loss_of_function_p3_vcf -O z -o $OutDir_introgressed/$P3.introgressed.loss_of_function.derived.recode.vcf.gz
$plink2 --vcf $OutDir_introgressed/$P3.introgressed.loss_of_function.derived.recode.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_introgressed/$P3.introgressed.loss_of_function.derived.hom_het

#synonymous
$bcftools index $synonymous_p3_vcf
$bcftools view -R $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed $synonymous_p3_vcf -O z -o $OutDir_introgressed/$P3.introgressed.synonymous.derived.recode.vcf.gz
$plink2 --vcf $OutDir_introgressed/$P3.introgressed.synonymous.derived.recode.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_introgressed/$P3.introgressed.synonymous.derived.hom_het

##summarize
echo -e "data\tdele_hom\tdele_het\tlof_hom\tlof_het\tsyn_hom\tsyn_het" > $OutDir_introgressed/$pairs.$P3.hom_het.summary.txt
p3_dele_hom=$(sed '1d' $OutDir_introgressed/$P3.introgressed.deleterious.derived.hom_het.scount |head -n -2 |$small_tools/sum col=3)
p3_dele_het=$(sed '1d' $OutDir_introgressed/$P3.introgressed.deleterious.derived.hom_het.scount |head -n -2 |$small_tools/sum col=4)
p3_lof_hom=$(sed '1d' $OutDir_introgressed/$P3.introgressed.loss_of_function.derived.hom_het.scount |head -n -2 |$small_tools/sum col=3)
p3_lof_het=$(sed '1d' $OutDir_introgressed/$P3.introgressed.loss_of_function.derived.hom_het.scount |head -n -2 |$small_tools/sum col=4)
p3_syno_hom=$(sed '1d' $OutDir_introgressed/$P3.introgressed.synonymous.derived.hom_het.scount |head -n -2 |$small_tools/sum col=3)
p3_syno_het=$(sed '1d' $OutDir_introgressed/$P3.introgressed.synonymous.derived.hom_het.scount |head -n -2 |$small_tools/sum col=4)

echo -e "real\t$p3_dele_hom\t$p3_dele_het\t$p3_lof_hom\t$p3_lof_het\t$p3_syno_hom\t$p3_syno_het" >> $OutDir_introgressed/$pairs.$P3.hom_het.summary.txt


#permutation
for i in {1..1000}
do

$bcftools view -R $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed $deleterious_p3_vcf -O z -o $OutDir_permutation/$P3.deleterious.hom_het.permu.$i.vcf.gz
$plink2 --vcf $OutDir_permutation/$P3.deleterious.hom_het.permu.$i.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_permutation/$P3.deleterious.hom_het.permu.$i

$bcftools view -R $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed $loss_of_function_p3_vcf -O z -o $OutDir_permutation/$P3.loss_of_function.hom_het.permu.$i.vcf.gz
$plink2 --vcf $OutDir_permutation/$P3.loss_of_function.hom_het.permu.$i.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_permutation/$P3.loss_of_function.hom_het.permu.$i

$bcftools view -R $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed $synonymous_p3_vcf -O z -o $OutDir_permutation/$P3.synonymous.hom_het.permu.$i.vcf.gz
$plink2 --vcf $OutDir_permutation/$P3.synonymous.hom_het.permu.$i.vcf.gz --sample-counts --allow-extra-chr --out $OutDir_permutation/$P3.synonymous.hom_het.permu.$i

permu_p3_dele_hom=$(sed '1d' $OutDir_permutation/$P3.deleterious.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=3)
permu_p3_dele_het=$(sed '1d' $OutDir_permutation/$P3.deleterious.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=4)

permu_p3_lof_hom=$(sed '1d' $OutDir_permutation/$P3.loss_of_function.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=3)
permu_p3_lof_het=$(sed '1d' $OutDir_permutation/$P3.loss_of_function.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=4)

permu_p3_syno_hom=$(sed '1d' $OutDir_permutation/$P3.synonymous.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=3)
permu_p3_syno_het=$(sed '1d' $OutDir_permutation/$P3.synonymous.hom_het.permu.$i.scount |head -n -2 |$small_tools/sum col=4)


echo -e "permutation_$i\t$permu_p3_dele_hom\t$permu_p3_dele_het\t$permu_p3_lof_hom\t$permu_p3_lof_het\t$permu_p3_syno_hom\t$permu_p3_syno_het" >> $OutDir_introgressed/$pairs.$P3.hom_het.summary.txt

done


elif [ $step == "10" ]; then
##selection sweeps
#ihs
#P2

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $selscan_input/populus162.phased.recode.$P2.snp.all.ihs.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.ihs.txt

##summarize
echo -e "data\taverage_ihs" > $OutDir_introgressed/$pairs.$P2.ihs.summary.txt
p2_ihs=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.ihs.txt |$small_tools/mean col=4)

echo -e "real\t$p2_ihs" >> $OutDir_introgressed/$pairs.$P2.ihs.summary.txt

#permutation
for i in {1..1000}
do


$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $selscan_input/populus162.phased.recode.$P2.snp.all.ihs.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.ihs.permu.$i.txt

permu_p2_ihs=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.ihs.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_ihs" >> $OutDir_introgressed/$pairs.$P2.ihs.summary.txt
done


#P3
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $selscan_input/populus162.phased.recode.$P3.snp.all.ihs.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.ihs.txt

echo -e "data\taverage_ihs" > $OutDir_introgressed/$pairs.$P3.ihs.summary.txt
p3_ihs=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.ihs.txt |$small_tools/mean col=4)

echo -e "real\t$p3_ihs" >> $OutDir_introgressed/$pairs.$P3.ihs.summary.txt

#permutation
for i in {1..1000}
do
$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $selscan_input/populus162.phased.recode.$P3.snp.all.ihs.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.ihs.permu.$i.txt

permu_p3_ihs=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.ihs.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p3_ihs" >> $OutDir_introgressed/$pairs.$P3.ihs.summary.txt
done



elif [ $step == "11" ]; then

##selection sweeps
#nsl
#P2

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $selscan_input/populus162.phased.recode.$P2.snp.all.nsl.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.nsl.txt

##summarize
echo -e "data\taverage_nsl" > $OutDir_introgressed/$pairs.$P2.nsl.summary.txt
p2_nsl=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.nsl.txt |$small_tools/mean col=4)

echo -e "real\t$p2_nsl" >> $OutDir_introgressed/$pairs.$P2.nsl.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $selscan_input/populus162.phased.recode.$P2.snp.all.nsl.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.nsl.permu.$i.txt

permu_p2_nsl=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.nsl.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_nsl" >> $OutDir_introgressed/$pairs.$P2.nsl.summary.txt
done

#P3
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $selscan_input/populus162.phased.recode.$P3.snp.all.nsl.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.nsl.txt

echo -e "data\taverage_nsl" > $OutDir_introgressed/$pairs.$P3.nsl.summary.txt
p3_nsl=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.nsl.txt |$small_tools/mean col=4)

echo -e "real\t$p3_nsl" >> $OutDir_introgressed/$pairs.$P3.nsl.summary.txt

#permutation
for i in {1..1000}
do
$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $selscan_input/populus162.phased.recode.$P3.snp.all.nsl.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.nsl.permu.$i.txt

permu_p3_nsl=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.nsl.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p3_nsl" >> $OutDir_introgressed/$pairs.$P3.nsl.summary.txt
done


elif [ $step == "12" ]; then

##selection sweeps
#ihh12
#P2

$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $selscan_input/populus162.phased.recode.$P2.snp.all.ihh12.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P2.ihh12.txt

##summarize
echo -e "data\taverage_ihh12" > $OutDir_introgressed/$pairs.$P2.ihh12.summary.txt
p2_ihh12=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P2.ihh12.txt |$small_tools/mean col=4)

echo -e "real\t$p2_ihh12" >> $OutDir_introgressed/$pairs.$P2.ihh12.summary.txt

#permutation
for i in {1..1000}
do

$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $selscan_input/populus162.phased.recode.$P2.snp.all.ihh12.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P2.ihh12.permu.$i.txt

permu_p2_ihh12=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P2.ihh12.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p2_ihh12" >> $OutDir_introgressed/$pairs.$P2.ihh12.summary.txt
done

#P3
$bedtools map -a $OutDir_introgressed/$pairs.introgressed.50_20.merge.bed -b $selscan_input/populus162.phased.recode.$P3.snp.all.ihh12.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_introgressed/$pairs.introgressed_islands.$P3.ihh12.txt

echo -e "data\taverage_ihh12" > $OutDir_introgressed/$pairs.$P3.ihh12.summary.txt
p3_ihh12=$(cat $OutDir_introgressed/$pairs.introgressed_islands.$P3.ihh12.txt |$small_tools/mean col=4)

echo -e "real\t$p3_ihh12" >> $OutDir_introgressed/$pairs.$P3.ihh12.summary.txt

#permutation
for i in {1..1000}
do
$bedtools map -a $OutDir_permutation/$pairs.introgressed.50_20.merge.permu.$i.bed -b $selscan_input/populus162.phased.recode.$P3.snp.all.ihh12.out.100bins.norm.bed -c 4 -o mean -g $potra_genome > $OutDir_permutation/$pairs.introgressed_islands.$P3.ihh12.permu.$i.txt

permu_p3_ihh12=$(cat $OutDir_permutation/$pairs.introgressed_islands.$P3.ihh12.permu.$i.txt |$small_tools/mean col=4)

echo -e "permutation_$i\t$permu_p3_ihh12" >> $OutDir_introgressed/$pairs.$P3.ihh12.summary.txt
done


elif [ $step == "13" ]; then
##compare the adaptive and non-adaptive introgressed regions

mv $OutDir_introgressed/*summary.txt $OutDir_introgressed_summary/

for file in $OutDir_introgressed/*.txt
do
Out=${file##*/}
Out_adaptive=${Out%.txt}.adaptive
Out_nonadaptive=${Out%.txt}.nonadaptive

#P2
$bedtools intersect -a $file -b $volcan_input/$P2.volcan.out.top5percentile.bed -wa | awk '! a[$1" "$2]++' > $OutDir_introgressed_adaptive/${Out_adaptive}.$P2.txt

$bedtools intersect -a $file -b $volcan_input/$P2.volcan.out.top5percentile.bed -v | awk '! a[$1" "$2]++' > $OutDir_introgressed_adaptive/${Out_nonadaptive}.$P2.txt

#P3
$bedtools intersect -a $file -b $volcan_input/$P3.volcan.out.top5percentile.bed -v | awk '! a[$1" "$2]++' > $OutDir_introgressed_adaptive/${Out_nonadaptive}.$P3.txt
$bedtools intersect -a $file -b $volcan_input/$P3.volcan.out.top5percentile.bed -wa | awk '! a[$1" "$2]++' > $OutDir_introgressed_adaptive/${Out_adaptive}.$P3.txt

done


fi


