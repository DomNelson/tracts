# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 14:33:24 2015

@author: dominic
"""
import matplotlib.pylab as pylab

import tracts_ped as ped
#import hmm_struct as struct
#import hmm_struct as struct
#import hmm_pedigree_old as old
import os

# MigrantProps = [0.2, 0.05] # Proportion of pedigree that will be new migrants
MigPropMat = [[15, 0, 0.5]]
SourceProps = [1, 0]
DemeSwitch = 0.1 # Chance of child leaving deme of parents
rho = 1 # Recombination rate
RecentMigrants = True
##@@ AllAncestry = True is assumed to assign parents ie no ancestry means you
## are not a leaf
AllAncestry = True
GlobalStepSize = 0.01
outfile = os.path.expanduser('test.jpg')

ChromLengths = [2.865747830, 2.64751457082595, 2.23363180733515, 
                2.15492839808593, 2.04089356863902, 1.92039918028429, 
                1.87852676459211, 1.68003441747308, 1.78206001355185, 
                1.81366917101923, 1.58218649890248, 1.74679023161126,
                1.26778791112187, 1.20202583329567, 1.39297570875973, 
                1.340377262456, 1.2849052927734, 1.17708922675517, 
                1.07733846085975, 1.08266933913055, 0.627864782064372, 
                0.741095623349923]
                

sim_ped = ped.Pedigree(sampleind = None, 
                 Gens = 15, 
                 RecentMigrants = RecentMigrants, 
                 AllAncestry = AllAncestry,
                 DemeSwitch = DemeSwitch,
                 SourceProps = SourceProps,
                 MigPropMat = MigPropMat)
                 
print "Creating individuals..."
#a.MakePedigree(DemeSwitch = DemeSwitch, MigrantProps = MigrantProps, 
    # RecentMigrants = RecentMigrants)
print "Sorting individuals..."
sim_ped.SortLeafNode()
# inds = sim_ped.indlist

print "Building chromosomes..."
sim_ped.MakeGenomes(ChromLengths = ChromLengths, rho = rho, smoothed = True,
                     Gamete = True)
                     
samp_ind = sim_ped.indlist[1]
samp_chrom1 = samp_ind.chromosomes['M0']
for tract in samp_chrom1.tracts:
    print tract.label
#struct.PlotChrom(samp_ind.chromosomes.values())
#gamete = sim_ped.indlist[0]
#struct.PlotChrom(gamete.chromosomes.values())



