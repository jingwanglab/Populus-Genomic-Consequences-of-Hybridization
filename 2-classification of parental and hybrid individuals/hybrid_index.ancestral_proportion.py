#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##################################################################
##the script used to calculate hybrid index and ancestral proportion based on parental species genotype and the genotype of all hybrids

import pandas as pd
import numpy as np
import sys
import getopt
import os

args=sys.argv

os.chdir(args[1])

#####read in the genotype file
joint_set=pd.read_csv(args[2],sep='\t')

n_col=joint_set.columns.size

columns=list(joint_set.columns)
hybrids=columns[4:]
species1=columns[2]
species2=columns[3]

species1_prop=[]  ###hybrid index when assuming species1 as 1 and speceis2 as 0

######calculating hybrid index (ancestral proportion of species1)
for i in range(len(hybrids)):
    species1_hom=[]
    for j in range(len(joint_set)):
        if joint_set[hybrids[i]][j] == joint_set[species1][j]:
            species1_hom.append(1.0)
        elif joint_set[hybrids[i]][j] == joint_set[species2][j]:
            species1_hom.append(0.0)
        elif joint_set[hybrids[i]][j] == "0/1":
            species1_hom.append(0.5)
        elif joint_set[hybrids[i]][j] == "./.":
            continue
    species1_prop.append((sum(species1_hom))/(len(species1_hom)))

#####calculating the heterozygosity of the potential hybrids
hets=[]

for i in range(len(hybrids)):
    het=[]
    for j in range(len(joint_set)):
        if joint_set[hybrids[i]][j] == "0/1":
            het.append(1.0)
        elif joint_set[hybrids[i]][j] == "./.":
            continue
        else:
            het.append(0.0)
    hets.append((sum(het))/(len(het)))            


######summarizing the above estimates into the table

ance_tab_ind=pd.DataFrame({'ind':hybrids,(species1+'_anc'):species1_prop,'het_prop':hets})
ance_tab_ind.to_csv((species1+'.'+species2+'.hybrids.'+species1+'_anc.hetero.txt'),index=False,sep='\t')





