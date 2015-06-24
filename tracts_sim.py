# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 14:33:24 2015

@author: dominic
"""
import matplotlib.pylab as pylab
import numpy as np
import tracts_ped as ped
#import hmm_struct as struct
#import hmm_struct as struct
#import hmm_pedigree_old as old
import os
import tracts
import sys
import cPickle

# MigrantProps = [0.2, 0.05] # Proportion of pedigree that will be new migrants
# MigPropMat = [[8, 0.1, 0], [12, 0, 0.1]]
DemeSwitch = 0.1 # Chance of child leaving deme of parents
rho = 1 # Recombination rate
##@@ AllAncestry = True is assumed to assign parents ie no ancestry means you
## are not a leaf
GlobalStepSize = 0.01
colordict = {0:'red', 1:'blue', 2:'green'}


ChromLengths = [2.865747830, 2.64751457082595, 2.23363180733515, 
                2.15492839808593, 2.04089356863902, 1.92039918028429, 
                1.87852676459211, 1.68003441747308, 1.78206001355185, 
                1.81366917101923, 1.58218649890248, 1.74679023161126,
                1.26778791112187, 1.20202583329567, 1.39297570875973, 
                1.340377262456, 1.2849052927734, 1.17708922675517, 
                1.07733846085975, 1.08266933913055, 0.627864782064372, 
                0.741095623349923]
                

migfile = sys.argv[1]
plotoutfile = os.path.expanduser(sys.argv[2])
numinds = int(sys.argv[3])
popoutfile = os.path.expanduser(sys.argv[4])
bedpath = os.path.expanduser(sys.argv[5])

migmat = np.genfromtxt(migfile)

indlist = []
for i in range(numinds):
    sim_ped = ped.Pedigree(sampleind = None, 
                     DemeSwitch = DemeSwitch,
                     MigPropMat = migmat)
    print "Simulating individual", i, "of", numinds                     
#    print "Creating individuals..."
    #a.MakePedigree(DemeSwitch = DemeSwitch, MigrantProps = MigrantProps, 
        # RecentMigrants = RecentMigrants)
#    print "Sorting individuals..."
    # sim_ped.SortLeafNode()
    # inds = sim_ped.indlist
    
#    print "Building chromosomes..."
    sim_ped.MakeGenomes(ChromLengths = ChromLengths, rho = rho, smoothed = True,
                         Gamete = False)
                         
    samp_ind = sim_ped.indlist[0]
    tracts_ind = samp_ind.to_tracts_indiv()
    indlist.append(tracts_ind)
    
    if bedpath != "None":
        outfile = os.path.join(bedpath, "IND" + str(i))
        samp_ind.to_bed_file(outfile)
        outfile = os.path.join(bedpath, "tracts_IND" + str(i))
        ped.tracts_ind_to_bed(tracts_ind, outfile)
#    samp_chrom1 = samp_ind.chromosomes['M0']
#    for tract in samp_chrom1.tracts:
#        print tract.label
    #struct.PlotChrom(samp_ind.chromosomes. values())
    #gamete = sim_ped.indlist[0]
    #struct.PlotChrom(gamete.chromosomes.values())
    
    
pop = tracts.population(list_indivs = indlist)
pop.plot_global_tractlengths(colordict, outfile = plotoutfile)

if popoutfile != "None":
    with open(popoutfile, 'wb') as f:
        cPickle.dump(pop, f, cPickle.HIGHEST_PROTOCOL)
        

