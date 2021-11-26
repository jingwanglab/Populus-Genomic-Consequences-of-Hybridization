#!/bin/sh
#SBATCH -c 1 --mem 10G
#SBATCH -o /UserHome/user260/pipeline/snpable/snpable.out
#SBATCH -e /UserHome/user260/pipeline/snpable/snpable.err


apply_mask_l="/UserHome/user260/tools/seqbility-20091110/apply_mask_l"
apply_mask_s="/UserHome/user260/tools/seqbility-20091110/apply_mask_s"
gen_mask="/UserHome/user260/tools/seqbility-20091110/gen_mask"
splitfa="/UserHome/user260/tools/seqbility-20091110/splitfa"
gen_raw_mask="/UserHome/user260/tools/seqbility-20091110/gen_raw_mask.pl"

bedtools="/data/apps/bedtools/bedtools-2.26.0/bin/bedtools"
bwa="/data/apps/bwa/bwa-0.7.17/bwa"
samtools="/data/apps/samtools/1.9/bin/samtools"
picard="/data/apps/picard/2.18.11/picard.jar"

ref="/w/user261/data/reference.pte/fasta/Potra02_genome.fasta"
bed="/w/user261/data/reference.pte/fasta/Potra02_genome.bed"
OutDir="/w/user261/data/snpable"


#######################################
###The first step is to create bed file that based on the information of uniquely mapepd regions by reads created by SNPable
#######################################
###Step1: split the genome into reads, with the length of each read is 100 bases
#$splitfa $ref 100 > $OutDir/Potra02_genome.reads

#######################################
###Step2: merge the reads together and align the short reads back to the genome with bwa aln

#$bwa aln -R 1000000 -O 3 -E 3 $ref $OutDir/Potra02_genome.reads > $OutDir/Potra02_genome.sai
#$bwa samse $ref $OutDir/Potra02_genome.sai $OutDir/Potra02_genome.reads > $OutDir/Potra02_genome.sam

#######################################
###Step3: generate the raw mask file and then set the stringency to make the final mask files
#perl $gen_raw_mask $OutDir/Potra02_genome.sam > $OutDir/Potra02_genome.rawMask_100.fa
#$gen_mask -l 100 -r 0.5 $OutDir/Potra02_genome.rawMask_100.fa > $OutDir/Potra02_genome.Mask_100_0.5.fa

#rm $OutDir_chr/*sam $OutDir_chr/*sai $OutDir_chr/*reads

#######################################
###Step4: generate the chr.txt which contains two columns (chr pos)
chr=$1 ##chr1,chr2..chr19,scaffold
bed="/w/user261/data/bed/Potra02_genome.$chr.bed"
#$bedtools getfasta -fi $OutDir/Potra02_genome.Mask_100_0.5.fa -bed $bed > $OutDir/Potra02_genome.$chr.Mask_100_0.5.fa

#if [ $chr == "scaffold" ]; then
##########################################################
##for scaffold
#scaff_n=$(cat $bed |wc -l)
#
#for ((i=1;i<=$scaff_n;i++))
#do
#scaff=$(head -n $i $bed |tail -n 1 |cut -f 1 )
#length=$(head -n $i $bed |tail -n 1 |cut -f 3 )
##echo $length

#for ((j=1;j<=$length;j++))
#do
#echo -e "$scaff\t$j" >> $OutDir/Potra02_genome.$chr.sites.txt
#done
#
#done
#
#else
##########################################################
#for chromosome 
##for chr1-chr19
#length=$(cut -f 3 $bed)
#echo $length
#for ((i=1;i<=$length;i++))
#do
#echo -e "chr$chr\t$i" >> $OutDir/Potra02_genome.chr$chr.sites.txt
#done

#fi

###modify the header of the $OutDir_chr/salmon.ssa$chr.Mask_100_0.5.fa file
#sed 's/:0-[0-9]*//' $OutDir/Potra02_genome.$chr.Mask_100_0.5.fa > $OutDir/temp && mv $OutDir/temp $OutDir/Potra02_genome.$chr.Mask_100_0.5.fa


###Step5: use several commands to transper the output from the above step to the bed file, which contains the information of the uniquely mapped positions on the specific chromosome

if [ $chr == "scaffold" ]; then

#$apply_mask_l $OutDir/Potra02_genome.$chr.Mask_100_0.5.fa $OutDir/Potra02_genome.$chr.sites.txt > $OutDir/Potra02_genome.$chr.mask.txt

scaff_n=$(cat $bed |wc -l)
echo $scaff_n

for ((i=1;i<=$scaff_n;i++))
do
scaff=$(head -n $i $bed |tail -n 1 |cut -f 1 )
echo $scaff
awk '$1=="'$scaff'"' $OutDir/Potra02_genome.$chr.mask.txt | awk 'NR==1{chr=$1;start=$2;end=$2;next} $2 == end+1 {end=$2;next} {print chr,start,end;start=$2;end=start} END{print chr,start,end}' >>  $OutDir/Potra02_genome.$chr.mask_mappability.txt
done

awk 'BEGIN{OFS = "\t"}{$2=$2-1;print}' $OutDir/Potra02_genome.$chr.mask_mappability.txt > $OutDir/Potra02_genome.$chr.mask_mappability.bed

else
#$apply_mask_l $OutDir/Potra02_genome.$chr.Mask_100_0.5.fa $OutDir/Potra02_genome.$chr.sites.txt > $OutDir/Potra02_genome.$chr.mask.txt

awk 'NR==1{chr=$1;start=$2;end=$2;next} $2 == end+1 {end=$2;next} {print chr,start,end;start=$2;end=start} END{print chr,start,end}' $OutDir/Potra02_genome.$chr.mask.txt > $OutDir/Potra02_genome.$chr.mask_mappability.txt
awk 'BEGIN{OFS = "\t"}{$2=$2-1;print}' $OutDir/Potra02_genome.$chr.mask_mappability.txt > $OutDir/Potra02_genome.$chr.mask_mappability.bed

fi

######copy the bed to a specific folder
#/data/apps/bedtools/bedtools-2.26.0/bin/bedtools coverage -a Potra02_genome.bed -b Potra02_genome.mask_mappability.bed > Potra02_genome.mask_mappability.map.txt
#the average unmask proportion of the genome is 328895043/408834716=0.8044695


