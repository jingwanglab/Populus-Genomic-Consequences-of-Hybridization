#!bin/bash
PopLDdecay=/w/00/g/g00/lsy/software/PopLDdecay/bin/PopLDdecay
Inputdir=/w/00/g/g00/lsy/snp/
vcftools=/data/apps/vcftools/vcftools-0.1.15/bin/vcftools

for species in {Pqio,Pade,Palb,Prot,Ptrs,Ptra,Pdav}
do
if [ ! -d "$Inputdir/LD/$species" ]; then
mkdir -p $Inputdir/LD/$species
fi


##distance 100 maf 0.05
$PopLDdecay -InVCF $Inputdir/$species.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf -MAF 0.05 -MaxDist 100 -OutType 5 -OutStat $Inputdir/LD/$species/$species.100kb.maf5.poplddecay.5

if [ ! -d "$Inputdir/tajimad/100kb" ]; then
mkdir -p $Inputdir/tajimad/100kb
fi

##tajimad
$vcftools --gzvcf $Inputdir/$species.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf.gz --TajimaD 100000 --out $Inputdir/tajimad/100kb/$species.100kb

if [ ! -d "$Inputdir/singleton" ]; then
mkdir -p $Inputdir/singleton
fi

$vcftools --vcf  $Inputdir/$species.snp.rm_indel.para_filter.biallelic.GQ30.max_miss20.bed.recode.vcf  --singletons --out $Inputdir/singleton/$species.sample
sed 1d $species.sample.singleton|cut -f 3,5|grep -v '^D'|sort|uniq -c > $Inputdir/singleton/$species.sample.singleton.count

done
