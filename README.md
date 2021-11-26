# Populus-Genomic-Consequences-of-Hybridization
Scripts for Liu et al (2021) Demographic history and natural selection shape patterns of deleterious mutation load and barriers to introgression across Populus genome. submitted

# Documentation of Scripts

# 1. Population structure and phylogenetic analyses

snpable.ptremula.sh - Scripts to use SNPable to mask the genome

populus227.pop_struc.sh - Scripts to perform population structure (fastSTRUCTURE and pca)

outgroup_merge_vcf.populus227.ld_pruned.sh - Scripts to integrate outgroup species and construct the NJ tree 


# 2. Classification of parental and hybrids individuals 

hybrid_index_heterozygosity.sh - Scripts to extract the fixed differences SNPs between pairs of species and estimate the hybrid index and heterozygosity for both parental species and the hybrid species

hybrid_index.ancestral_proportion.py - Python scripts to calculate the hybrid index and ancestral proportion

chro.populus227.phylogenetic.sh - Scripts to construct the chloroplast NJ tree

mito.populus227.phylogenetic.sh - Scripts to construct the mitochondrial NJ tree 


# 3. Estimation of genetic variation, LD and species demographic history

pixy.pi_fst_dxy.100kb.sh - Scripts to estimate the nucleotide diversity (pi) and genetic divergence (Fst and dxy) over 100 Kbp non-overlapping windows

ld.singleton.tajimad.sh - Scripts to perform LD, singleton and Tajima's D estimates


# 4. Assessment of genetic load among species

pixy.pi.0_4fold.intergenic.sh - Scripts to estimate the heterozygosity at zero-fold, four-fold and intergenic sites per individual

dfe_0_4_fold.sh - Scripts to perform DFE analysis

sift4g_database.pte2.sh - Scripts to build SIFT4G database for reference genome

sift4g_mutation_load.sh - Scripts to estimate the relative proportion of homozygous and heterozygous derived alleles for LoF, deleterious, tolerated and synonymous sites per individual



# 5. Detection of past introgression across the genome

twisst.populus.sh - Scripts to use TWISST to perform analysis over 10Kbp windows

twisst.populus.sites1000.sh - Scripts to use TWISST to perform analysis over 1000 SNPs

dsuite.sh - Scripts to use Dsuite to detect calculate D and f4 statistics

introgressed_islands.dsuite.sh - Scripts to identify introgressed regions in each trio and also compare the observed values of various statistics to 1000 randomizations

volcanofinder.populus.new.sh - Scripts to use VolcanoFinder to detect signals of adaptive introgression sweeps

selscan.populus7species.sh - Scripts to use selscan to calculate the iHH12 statistic



