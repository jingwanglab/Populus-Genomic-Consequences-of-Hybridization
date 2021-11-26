#!/bin/bash
#SBATCH -c 1 --mem 30G
#SBATCH -o /UserHome/user260/pipeline/snp_effects/sift4g/sift4g.pte.out
#SBATCH -e /UserHome/user260/pipeline/snp_effects/sift4g/sift4g.pte.err

#######################################################################
###The main aim of this script is to first build the database of the Populus tremula v2.2 reference genome

export PERL5LIB=/data/apps/cpan/lib/perl5:/data/apps/cpan/lib64/perl5:/data/apps/cpan/share/perl5
step=$1
species=$2


sift4g_annotator="/UserHome/user260/tools/sift4g_2/sift4g/SIFT4G_Annotator.jar"

OutDir="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift4g"

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

OutDir_gtf="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift4g/gene-annotation-src"

if [ ! -d "$OutDir_gtf" ]; then
mkdir -p $OutDir_gtf
fi

OutDir_fasta="/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift4g/chr-src"

if [ ! -d "$OutDir_fasta" ]; then
mkdir -p $OutDir_fasta
fi


if [ $step == "1" ]; then
##################################
###Step1: becuase the gff file of Populus tremula v2.2 do not have start_codon and stop_codon features, I need first to use GAG.py to add them
gag="/UserHome/user260/tools/GAG-master/gag.py"
gff="/w/user260/multiple_aspen/gff/Potra02_genes.gff"
fa="/w/user260/multiple_aspen/gff/Potra02_genome.fasta"

python2.7 $gag --fasta $fa --gff $gff --fix_start_stop --out /w/user260/multiple_aspen/gff/start_stop_fixed

####move the fasta and gff file in ./start_stop_fix to the corresponding folder

cat /w/user260/multiple_aspen/gff/start_stop_fixed/genome.comments.gff /w/user260/multiple_aspen/gff/start_stop_fixed/genome.gff > $OutDir_gtf/Potra02.genes.gff3
cp /w/user260/multiple_aspen/gff/Potra02_proteins.fasta $OutDir_gtf/Potra02.pep.fa
cp /w/user260/multiple_aspen/gff/Potra02_genome.fasta $OutDir_gtf/Potra02.genome.fa

elif [ $step == "2" ]; then
##################################
##Step2: transfer gff to gtf file to used by SIFT
gffread="/UserHome/user260/tools/gffread/gffread"

$gffread $OutDir_gtf/Potra02.genes.gff3 -T -o $OutDir_gtf/Potra02.genes.gtf

elif [ $step == "3" ]; then
##################################
##Step3:write config files to each chromsome folder

echo -e "" > $OutDir/pte.sift.config.txt
echo -e "GENETIC_CODE_TABLE=1" >> $OutDir/pte.sift.config.txt
echo -e "GENETIC_CODE_TABLENAME=Standard" >> $OutDir/pte.sift.config.txt
echo -e "MITO_GENETIC_CODE_TABLE=11" >> $OutDir/pte.sift.config.txt
echo -e "MITO_GENETIC_CODE_TABLENAME=Plant Plastid Code" >> $OutDir/pte.sift.config.txt
echo -e "" >> $OutDir/pte.sift.config.txt
echo -e "PARENT_DIR=$OutDir" >> $OutDir/pte.sift.config.txt
echo -e "ORG=Populus_tremula" >> $OutDir/pte.sift.config.txt
echo -e "ORG_VERSION=v2" >> $OutDir/pte.sift.config.txt
echo -e "DBSNP_VCF_FILE=" >> $OutDir/pte.sift.config.txt
echo -e "" >> $OutDir/pte.sift.config.txt

echo -e "#Running SIFT 4G" >> $OutDir/pte.sift.config.txt
echo -e "SIFT4G_PATH=/UserHome/user260/tools/sift4g_2/sift4g/bin/sift4g" >> $OutDir/pte.sift.config.txt
echo -e "PROTEIN_DB=/w/user260/multiple_aspen/populus162_vcf/vcf/species_vcf/sift4g/gene-annotation-src/uniref90.fasta" >> $OutDir/pte.sift.config.txt
echo -e "COMPUTER=GIS-KATNISS" >> $OutDir/pte.sift.config.txt
echo -e "" >> $OutDir/pte.sift.config.txt
echo -e "# Sub-directories, don't need to change" >> $OutDir/pte.sift.config.txt

echo -e "GENE_DOWNLOAD_DEST=gene-annotation-src" >> $OutDir/pte.sift.config.txt
echo -e "CHR_DOWNLOAD_DEST=chr-src" >> $OutDir/pte.sift.config.txt
echo -e "LOGFILE=Log.txt" >> $OutDir/pte.sift.config.txt
echo -e "ZLOGFILE=Log2.txt" >> $OutDir/pte.sift.config.txt
echo -e "FASTA_DIR=fasta" >> $OutDir/pte.sift.config.txt
echo -e "SUBST_DIR=subst" >> $OutDir/pte.sift.config.txt
echo -e "ALIGN_DIR=SIFT_alignments" >> $OutDir/pte.sift.config.txt
echo -e "SIFT_SCORE_DIR=SIFT_predictions" >> $OutDir/pte.sift.config.txt
echo -e "SINGLE_REC_BY_CHR_DIR=singleRecords" >> $OutDir/pte.sift.config.txt
echo -e "SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores" >> $OutDir/pte.sift.config.txt
echo -e "DBSNP_DIR=dbSNP" >> $OutDir/pte.sift.config.txt
echo -e "" >> $OutDir/pte.sift.config.txt
echo -e "# Doesn't need to change" >> $OutDir/pte.sift.config.txt
echo -e "FASTA_LOG=fasta.log" >> $OutDir/pte.sift.config.txt
echo -e "INVALID_LOG=invalid.log" >> $OutDir/pte.sift.config.txt
echo -e "PEPTIDE_LOG=peptide.log" >> $OutDir/pte.sift.config.txt
echo -e "ENS_PATTERN=ENS" >> $OutDir/pte.sift.config.txt
echo -e "SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord" >> $OutDir/pte.sift.config.txt

elif [ $step == "4" ]; then
#######################################

#run config files for each chromosome
perl_build="/UserHome/user260/tools/sift4g_2/SIFT4G_Create_Genomic_DB-master/make-SIFT-db-all.pl"
config="$OutDir/pte.sift.config.txt"

cd /UserHome/user260/tools/sift4g_2/SIFT4G_Create_Genomic_DB-master
perl make-SIFT-db-all.pl --config $config


fi

