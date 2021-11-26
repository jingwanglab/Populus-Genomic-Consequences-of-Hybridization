#!/bin/sh
#SBATCH -c 1 --mem 16G
#SBATCH -o /UserHome/user260/pipeline/hybridization/Dsuite/dsuite.populus190.out
#SBATCH -e /UserHome/user260/pipeline/hybridization/Dsuite/dsuite.populus190.err

###running Dsuite among trios of species pairs

vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
dsuite="/UserHome/user260/tools/dsuite/Dsuite/Build/Dsuite"
intersect_R="/UserHome/user260/pipeline/hybridization/Dsuite/intersect.snp.R"

step=$1

######The first step is to merge the two outgroup vcf with the populus 190 vcf files
tri_vcf="/w/user260/multiple_aspen/outgroups/trichocarpa1.snp.vcf.gz"
peu_vcf="/w/user260/multiple_aspen/outgroups/PeuXBY08.snp.vcf.gz"

populus190_vcf="/w/user260/multiple_aspen/populus190_vcf/populus190.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz"

populus190_samples="/w/user260/multiple_aspen/populus190_vcf/samples/populus190.samples.remove_clones.csv"


OutDir="/w/user260/multiple_aspen/populus190_vcf/dsuite"

OutDir_vcf=$OutDir/vcf
if [ ! -d "$OutDir_vcf" ]; then
mkdir -p $OutDir_vcf
fi

OutDir_dsuite=$OutDir/dsuite
if [ ! -d "$OutDir_dsuite" ]; then
mkdir -p $OutDir_dsuite
fi

OutDir_tree=$OutDir_dsuite/tree_based
if [ ! -d "$OutDir_tree" ]; then
mkdir -p $OutDir_tree
fi

OutDir_tree2=$OutDir_dsuite/tree_based2
if [ ! -d "$OutDir_tree2" ]; then
mkdir -p $OutDir_tree2
fi

OutDir_dtrios=$OutDir_dsuite/dtrios
if [ ! -d "$OutDir_dtrios" ]; then
mkdir -p $OutDir_dtrios
fi

OutDir_dquartets=$OutDir_dsuite/dquartets
if [ ! -d "$OutDir_dquartets" ]; then
mkdir -p $OutDir_dquartets
fi

OutDir_dinvestigate=$OutDir_dsuite/dinvestigate
if [ ! -d "$OutDir_dinvestigate" ]; then
mkdir -p $OutDir_dinvestigate
fi


if [ $step == "1" ]; then
#extract the snp from the populus190 samples
#zcat $populus190_vcf | grep -v "#" | cut -f 3  > $OutDir_vcf/populus190.vcf.snp

##trichocarpa1
$vcftools --gzvcf $tri_vcf --snps $OutDir_vcf/populus190.vcf.snp --recode --recode-INFO-all --out $OutDir_vcf/tri.snp.populus190
$bcftools view $OutDir_vcf/tri.snp.populus190.recode.vcf -O z -o $OutDir_vcf/tri.snp.populus190.recode.vcf.gz && rm $OutDir_vcf/tri.snp.populus190.recode.vcf

##PeuXBY08
$vcftools --gzvcf $peu_vcf --snps $OutDir_vcf/populus190.vcf.snp --recode --recode-INFO-all --out $OutDir_vcf/PeuXBY08.snp.populus190
$bcftools view $OutDir_vcf/PeuXBY08.snp.populus190.recode.vcf -O z -o $OutDir_vcf/PeuXBY08.snp.populus190.recode.vcf.gz && rm $OutDir_vcf/PeuXBY08.snp.populus190.recode.vcf

elif [ $step == "2" ]; then
####only keep the non-missing sites of the 

###for the outgroup, only keep the non-missing and biallelic SNPs
#$bcftools view -g ^miss -O z -o $OutDir_vcf/tri.no_miss.populus190.vcf.gz $OutDir_vcf/tri.snp.populus190.recode.vcf.gz
#$bcftools view -m0 -M2 -V indels -O z -o $OutDir_vcf/tri.no_miss.biallelic.populus190.vcf.gz $OutDir_vcf/tri.no_miss.populus190.vcf.gz
#$bcftools index $OutDir_vcf/tri.no_miss.biallelic.populus190.vcf.gz
#zcat $OutDir_vcf/tri.no_miss.biallelic.populus190.vcf.gz |grep -v "#" |cut -f 3 > $OutDir_vcf/tri.final.snp

#$bcftools view -g ^miss -O z -o $OutDir_vcf/PeuXBY08.no_miss.populus190.vcf.gz $OutDir_vcf/PeuXBY08.snp.populus190.recode.vcf.gz
#$bcftools view -m0 -M2 -V indels -O z -o $OutDir_vcf/PeuXBY08.no_miss.biallelic.populus190.vcf.gz $OutDir_vcf/PeuXBY08.no_miss.populus190.vcf.gz
#$bcftools index $OutDir_vcf/PeuXBY08.no_miss.biallelic.populus190.vcf.gz
#zcat $OutDir_vcf/PeuXBY08.no_miss.biallelic.populus190.vcf.gz |grep -v "#" |cut -f 3 > $OutDir_vcf/PeuXBY08.final.snp

#Rscript $intersect_R $OutDir_vcf
#$bcftools index $populus190_vcf

$bcftools annotate -x FORMAT -O z -o $OutDir_vcf/populus190.gt.vcf.gz $populus190_vcf
$vcftools --gzvcf $OutDir_vcf/populus190.gt.vcf.gz --snps $OutDir_vcf/intersect.snp --recode --recode-INFO-all --out $OutDir_vcf/populus190.isec
$bcftools view $OutDir_vcf/populus190.isec.recode.vcf -O z -o $OutDir_vcf/populus190.isec.recode.vcf.gz && rm $OutDir_vcf/populus190.isec.recode.vcf
$bcftools index $OutDir_vcf/populus190.isec.recode.vcf.gz
#
$bcftools annotate -x FORMAT -O z -o $OutDir_vcf/tri.no_miss.gt.populus190.vcf.gz $OutDir_vcf/tri.no_miss.biallelic.populus190.vcf.gz
$vcftools --gzvcf $OutDir_vcf/tri.no_miss.gt.populus190.vcf.gz --snps $OutDir_vcf/intersect.snp --recode --recode-INFO-all --out $OutDir_vcf/tri.no_miss.populus190.isec
$bcftools view $OutDir_vcf/tri.no_miss.populus190.isec.recode.vcf -O z -o $OutDir_vcf/tri.no_miss.populus190.isec.recode.vcf.gz && rm $OutDir_vcf/tri.no_miss.populus190.isec.recode.vcf
$bcftools index $OutDir_vcf/tri.no_miss.populus190.isec.recode.vcf.gz

$bcftools annotate -x FORMAT -O z -o $OutDir_vcf/PeuXBY08.no_miss.gt.populus190.vcf.gz $OutDir_vcf/PeuXBY08.no_miss.biallelic.populus190.vcf.gz
$vcftools --gzvcf $OutDir_vcf/PeuXBY08.no_miss.gt.populus190.vcf.gz --snps $OutDir_vcf/intersect.snp --recode --recode-INFO-all --out $OutDir_vcf/PeuXBY08.no_miss.populus190.isec
$bcftools view $OutDir_vcf/PeuXBY08.no_miss.populus190.isec.recode.vcf -O z -o $OutDir_vcf/PeuXBY08.no_miss.populus190.isec.recode.vcf.gz && rm $OutDir_vcf/PeuXBY08.no_miss.populus190.isec.recode.vcf
$bcftools index $OutDir_vcf/PeuXBY08.no_miss.populus190.isec.recode.vcf.gz

$bcftools merge $OutDir_vcf/populus190.isec.recode.vcf.gz $OutDir_vcf/tri.no_miss.populus190.isec.recode.vcf.gz $OutDir_vcf/PeuXBY08.no_miss.populus190.isec.recode.vcf.gz -O z -o $OutDir_vcf/populus190.outgroup.vcf.gz


elif [ $step == "2.1" ]; then

$bcftools view -m2 -M2 -v snps -O z -o $OutDir_vcf/populus190.outgroup.biallelic.vcf.gz $OutDir_vcf/populus190.outgroup.vcf.gz

elif [ $step == "3" ]; then

### writing the input SETS.TXT file
awk -F "," '{if($7~/Pure/)print $4,$2}' $populus190_samples | sed 's/ /\t/g'  > $OutDir_dsuite/populus190.dsuite.txt
echo -e "trichocarpa1\tOutgroup" >> $OutDir_dsuite/populus190.dsuite.txt
echo -e "PeuXBY08\tOutgroup" >> $OutDir_dsuite/populus190.dsuite.txt

### writing the optional Newick file
echo -e "(Pqio,(Pade,(Palb,(Ptrs,(Ptra,(Pdav,Prot))))))" >> $OutDir_dsuite/tree.nwk

elif [ $step == "4" ]; then

###running Dsuite with tree based files
export LD_LIBRARY_PATH=/data/apps/gmp/gmp-6.1.2/lib:/data/apps/mpfr/mpfr-3.1.5/lib:/data/apps/mpc/mpc-1.0.3/lib:/data/apps/isl/isl-0.19/lib:/data/apps/gcc/8.3.0/lib64
export PATH=/data/apps/gcc/8.3.0/bin:$PATH

$dsuite Dtrios -t $OutDir_dsuite/tree.nwk -o $OutDir_tree/populus190_outgroup $OutDir_vcf/populus190.outgroup.vcf.gz $OutDir_dsuite/populus190.dsuite.txt

elif [ $step == "4.2" ]; then

###using the tree nwk file with (Pade,Pqio) being together 
export LD_LIBRARY_PATH=/data/apps/gmp/gmp-6.1.2/lib:/data/apps/mpfr/mpfr-3.1.5/lib:/data/apps/mpc/mpc-1.0.3/lib:/data/apps/isl/isl-0.19/lib:/data/apps/gcc/8.3.0/lib64
export PATH=/data/apps/gcc/8.3.0/bin:$PATH

$dsuite Dtrios -t $OutDir_dsuite/tree2.nwk -o $OutDir_tree2/populus190_outgroup $OutDir_vcf/populus190.outgroup.vcf.gz $OutDir_dsuite/populus190.dsuite.txt


elif [ $step == "5" ]; then

###running Dsuite withnot tree based files
export LD_LIBRARY_PATH=/data/apps/gmp/gmp-6.1.2/lib:/data/apps/mpfr/mpfr-3.1.5/lib:/data/apps/mpc/mpc-1.0.3/lib:/data/apps/isl/isl-0.19/lib:/data/apps/gcc/8.3.0/lib64
export PATH=/data/apps/gcc/8.3.0/bin:$PATH

$dsuite Dtrios -o $OutDir_dtrios/populus190_outgroup $OutDir_vcf/populus190.outgroup.vcf.gz $OutDir_dsuite/populus190.dsuite.txt

elif [ $step == "6" ]; then

export LD_LIBRARY_PATH=/data/apps/gmp/gmp-6.1.2/lib:/data/apps/mpfr/mpfr-3.1.5/lib:/data/apps/mpc/mpc-1.0.3/lib:/data/apps/isl/isl-0.19/lib:/data/apps/gcc/8.3.0/lib64
export PATH=/data/apps/gcc/8.3.0/bin:$PATH

###running Dsuite Dquartets
$dsuite Dquartets -o $OutDir_dquartets/populus190_dquartets $OutDir_vcf/populus190.outgroup.vcf.gz $OutDir_dsuite/populus190.dsuite.no_outgroup.txt

elif [ $step == "7" ]; then

###running Dinvestigate for trios of populations with significantly elevated D; with p-value<0.001
export LD_LIBRARY_PATH=/data/apps/gmp/gmp-6.1.2/lib:/data/apps/mpfr/mpfr-3.1.5/lib:/data/apps/mpc/mpc-1.0.3/lib:/data/apps/isl/isl-0.19/lib:/data/apps/gcc/8.3.0/lib64
export PATH=/data/apps/gcc/8.3.0/bin:$PATH

####extracting the trios information of the output file 
#aawk '$6<0.001' $OutDir_tree2/populus190_outgroup_BBAA.txt |cut -f 1-3 > $OutDir_dinvestigate/dinvestigate.trios.txt
n=$2
trios=$(head -n $n $OutDir_dinvestigate/dinvestigate.trios.txt |tail -n 1 |sed 's/\t/_/g')
echo $trios
head -n $n $OutDir_dinvestigate/dinvestigate.trios.txt |tail -n 1 > $OutDir_dinvestigate/$trios.txt
#$dsuite Dinvestigate -w 100 50 

#$dsuite Dinvestigate -w 100,50 $OutDir_vcf/populus190.outgroup.vcf.gz $OutDir_dsuite/populus190.dsuite.txt  $OutDir_dinvestigate/$trios.txt
#$dsuite Dinvestigate -w 50,20 $OutDir_vcf/populus190.outgroup.vcf.gz $OutDir_dsuite/populus190.dsuite.txt  $OutDir_dinvestigate/$trios.txt
$dsuite Dinvestigate -w 20,10 $OutDir_vcf/populus190.outgroup.vcf.gz $OutDir_dsuite/populus190.dsuite.txt  $OutDir_dinvestigate/$trios.txt

elif [ $step == "8" ]; then

##running Fbranch
export LD_LIBRARY_PATH=/data/apps/gmp/gmp-6.1.2/lib:/data/apps/mpfr/mpfr-3.1.5/lib:/data/apps/mpc/mpc-1.0.3/lib:/data/apps/isl/isl-0.19/lib:/data/apps/gcc/8.3.0/lib64
export PATH=/data/apps/gcc/8.3.0/bin:$PATH

$dsuite Fbranch $OutDir_dsuite/tree.nwk $OutDir_tree/populus190_outgroup__tree.txt > $OutDir_tree/populus190_outgroup.fbranch.txt 
python3.6 /UserHome/user260/tools/dsuite/Dsuite/utils/dtools.py $OutDir_tree/populus190_outgroup.fbranch.txt $OutDir_dsuite/tree.nwk

fi




