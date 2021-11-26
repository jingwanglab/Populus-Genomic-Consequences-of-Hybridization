#!/bin/sh
#SBATCH -c 2 --mem 30G
#SBATCH -o /UserHome/user260/pipeline/hybridization/twisst/twisst.populus.out
#SBATCH -e /UserHome/user260/pipeline/hybridization/twisst/twisst.populus.err


vcftools="/data/apps/vcftools/vcftools-0.1.15/bin/vcftools"
bcftools="/data/apps/bcftools/1.9/bin/bcftools"
bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
pigz="/UserHome/user260/tools/pigz-2.4/pigz"
raxml="/data/apps/raxml/RAxML-8.2.11/raxmlHPC-PTHREADS"
raxml_tool="/UserHome/user260/tools/genomics_general-master/phylo/raxml_sliding_windows.py"
phyml="/data/apps/phyml/3.3.20190909/bin/phyml"
phyml_tool="/UserHome/user260/tools/genomics_general-master/phylo/phyml_sliding_windows.py"
twisst_tool="/UserHome/user260/tools/twisst-master/twisst.py"
parse_vcf="/UserHome/user260/tools/genomics_general-master/VCF_processing/parseVCF.py"
python="/data/apps/python/3.8.0/bin/python3.8"


Inputvcf="/w/user260/multiple_aspen/populus190_vcf/beagle/populus190.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.gt.vcf.gz"
Out=${Inputvcf##*/}
InputDir=`dirname $Inputvcf`

populus162="/w/user260/multiple_aspen/populus162_vcf/samples/populus162.samples.txt"


twisst_dir=$InputDir/twisst

if [ ! -d "$twisst_dir" ]; then
mkdir -p $twisst_dir
fi

weighting_dir=$twisst_dir/weightings

if [ ! -d "$weighting_dir" ]; then
mkdir -p $weighting_dir
fi



step=$1
chr=$2
window=$3

bed="/w/user260/multiple_aspen/snpable/Potra02_genome.$chr.bed"


if [ "$step" == "0" ]; then
###extract the 162 samples 

$vcftools --gzvcf $Inputvcf --keep $populus162 --maf 0.000001 --recode --recode-INFO-all --out $InputDir/populus162.phased && $pigz $InputDir/populus162.phased.recode.vcf

elif [ "$step" == "1" ]; then
###extract the chr separately
$vcftools --gzvcf $InputDir/populus162.phased.recode.vcf.gz --chr $chr --recode --recode-INFO-all --out $InputDir/populus162.phased.bed.$chr && $pigz $InputDir/populus162.phased.bed.$chr.recode.vcf

elif [ "$step" == "2" ]; then
###processing vcf files and transfer the .vcf to .geno file

$python $parse_vcf -i $InputDir/populus162.phased.bed.$chr.recode.vcf.gz |$pigz > $twisst_dir/populus162.$chr.geno.gz

elif [ "$step" == "3" ]; then
###using raxml to make phylogenetic trees

if [ "$window" == "100000" ]; then
##100kb
#$python $raxml_tool -T 10 --genoFile $twisst_dir/populus162.$chr.geno.gz --windType coordinate -w 100000 --raxml $raxml --model GTRCATI -Ms 200 --prefix $twisst_dir/populus162.$chr.raxml.w100kb
#$python $raxml_tool -T 10 --genoFile $twisst_dir/populus162.$chr.geno.gz --windType coordinate -w 100000 --raxml $raxml --model GTRCATI -Ms 200 --prefix $twisst_dir/populus162.$chr.raxml.w100kb
$python $phyml_tool -g $twisst_dir/populus162.$chr.geno.gz --prefix $twisst_dir/populus162.$chr.phyml.w100kb -w 100000 --windType coordinate --model GTR --phyml $phyml -T 16 -M 200

elif [ "$window" == "10000" ]; then
##10kb
$python $phyml_tool -g $twisst_dir/populus162.$chr.geno.gz --prefix $twisst_dir/populus162.$chr.phyml.w10kb -w 10000 --windType coordinate --model GTR --phyml $phyml -T 16 -M 50 

fi

elif [ "$step" == "4" ]; then
##sites
$python $phyml_tool -g $twisst_dir/populus162.$chr.geno.gz --prefix $twisst_dir/populus162.$chr.phyml.sites100 -w 200 --windType sites --model GTR --phyml $phyml -T 16

elif [ "$step" == "5" ]; then
##twisst
##use awk '{sub("$", "_A", $1)}; 1' populus162.species.txt > populus162.species.a.txt
##use awk '{sub("$", "_B", $1)}; 1' populus162.species.txt > populus162.species.b.txt to add _A, _B labels 
group="/w/user260/multiple_aspen/populus190_vcf/beagle/twisst/samples/populus162.species.a_b.txt"


gunzip $twisst_dir/populus162.$chr.phyml.w10kb.trees.gz
sed '1d' $twisst_dir/populus162.$chr.phyml.w10kb.data.tsv > $twisst_dir/populus162.$chr.phyml.w10kb.data
paste $twisst_dir/populus162.$chr.phyml.w10kb.data $twisst_dir/populus162.$chr.phyml.w10kb.trees |awk '$5>50' |cut -f 7- |gzip - > $twisst_dir/populus162.$chr.phyml.w10kb.shrink.trees.gz
gzip $twisst_dir/populus162.$chr.phyml.w10kb.trees
rm $twisst_dir/populus162.$chr.phyml.w10kb.data
#
cp $twisst_dir/populus162.$chr.phyml.w10kb.shrink.trees.gz $weighting_dir

#$python $twisst_tool -t $twisst_dir/populus162.$chr.phyml.w10kb.shrink.trees.gz -w $weighting_dir/populus162.$chr.phyml.w10kb.twisst.weights.csv.gz -o $weighting_dir/populus162.$chr.phyml.w10kb.twisst.topologies.trees -g Pade -g Pqio -g Pdav -g Prot -g Ptra -g Ptrs -g Palb --method fixed --iterations 50000 --groupsFile $group -D $weighting_dir/populus162.$chr.phyml.w10kb.twisst.distance.csv.gz

$python $twisst_tool -t $twisst_dir/populus162.$chr.phyml.w10kb.shrink.trees.gz -w $weighting_dir/populus162.$chr.phyml.w10kb.twisst.weights.csv.gz -g Pade -g Pqio -g Pdav -g Prot -g Ptra -g Ptrs -g Palb --method fixed --iterations 100000 --groupsFile $group -D $weighting_dir/populus162.$chr.phyml.w10kb.twisst.distance.csv.gz

elif [ "$step" == "6" ]; then
##summarize the twisst results

awk '$5>50' $twisst_dir/populus162.chr1.phyml.w10kb.data.tsv > $weighting_dir/chr1.data.tsv
zcat $weighting_dir/populus162.chr1.phyml.w10kb.twisst.weights.csv.gz |grep -v "#" > $weighting_dir/chr1.twisst.weights.csv

paste $weighting_dir/chr1.data.tsv $weighting_dir/chr1.twisst.weights.csv > $weighting_dir/populus162.chr_all.phyml.w10kb.twisst.weights.csv && rm $weighting_dir/chr1.data.tsv $weighting_dir/chr1.twisst.weights.csv


for chr in chr{2..19}
do
awk '$5>50' $twisst_dir/populus162.$chr.phyml.w10kb.data.tsv |sed '1d' > $weighting_dir/$chr.data.tsv
zcat $weighting_dir/populus162.$chr.phyml.w10kb.twisst.weights.csv.gz |grep -v "#" |sed '1d' > $weighting_dir/$chr.twisst.weights.csv

paste $weighting_dir/$chr.data.tsv $weighting_dir/$chr.twisst.weights.csv >> $weighting_dir/populus162.chr_all.phyml.w10kb.twisst.weights.csv && rm $weighting_dir/$chr.data.tsv $weighting_dir/$chr.twisst.weights.csv
done

fi






