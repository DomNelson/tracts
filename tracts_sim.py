# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 14:33:24 2015

@author: dominic
"""
import matplotlib.pylab as pylab
import numpy as np
import tracts_ped as ped
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
labels = ['EUR', 'NAT', 'AFR']
colordict = {'EUR':'red', 'NAT':'blue', 'AFR':'green'}
#colordict = {'0':'red', '1':'blue', '2':'green'}
#colordict = {'AFR':'red', 'EUR':'blue'}#, 2:'green'}



#ChromLengths = [2.865747830, 2.64751457082595, 2.23363180733515, 
#                2.15492839808593, 2.04089356863902, 1.92039918028429, 
#                1.87852676459211, 1.68003441747308, 1.78206001355185, 
#                1.81366917101923, 1.58218649890248, 1.74679023161126,
#                1.26778791112187, 1.20202583329567, 1.39297570875973, 
#                1.340377262456, 1.2849052927734, 1.17708922675517, 
#                1.07733846085975, 1.08266933913055, 0.627864782064372, 
#                0.741095623349923]
                

cM_ChromLengths = [ 277.6846783 ,  263.4266571 ,  224.5261258 ,  212.8558223 ,
                203.9634184 ,  192.9822446 ,  186.9212679 ,  170.2156421 ,
                168.2431216 ,  179.0947462 ,  159.5132079 ,  172.8505693 ,
                126.9025447 ,  116.3957107 ,  131.405539  ,  134.9600594 ,
                129.2943145 ,  119.0324459 ,  107.8670432 ,  108.0521931 ,
                61.46827149,   72.68689882]

ChromLengths = [length / 100. for length in cM_ChromLengths]

migfile = sys.argv[1]
migmat = np.genfromtxt(migfile)
numinds = int(sys.argv[2])
method = sys.argv[3]
bedpath = os.path.expanduser(sys.argv[4])
try:
    popoutfile = os.path.expanduser(sys.argv[5])
    plotoutfile = os.path.expanduser(sys.argv[6])
except IndexError:
    popoutfile = "None"
    plotoutfile = "None"


indlist = []
for i in range(numinds):
    print "Simulating individual", i, "of", numinds                     
    if method == "forward":
        sim_ped = ped.Pedigree(sampleind = None, 
                         DemeSwitch = DemeSwitch,
                         MigPropMat = migmat,
                         labels = labels)
        sim_ped.MakeGenomes(ChromLengths = ChromLengths, rho = rho, smoothed = True,
                             Gamete = False)
        samp_ind = sim_ped.indlist[0]
        tracts_ind = samp_ind.to_tracts_indiv()
    elif method == "PSMC":
        P = ped.Pedigree(migmat, labels = labels)
        P.SortLeafNode()
        P.BuildTransMatrices()
        tracts_ind = P.PSMC_ind(ChromLengths)

    indlist.append(tracts_ind)
    
    if bedpath != "None":
        if not os.path.exists(bedpath):
            os.makedirs(bedpath)
#        outfile = os.path.join(bedpath, "IND" + str(i + 1))
#        samp_ind.to_bed_file(outfile)
        outfile = os.path.join(bedpath, "IND" + str(i + 1))
        ped.tracts_ind_to_bed(tracts_ind, outfile, conv = "M->cM")
    
pop = tracts.population(list_indivs = indlist)
pop.plot_global_tractlengths(colordict, outfile = plotoutfile)

if popoutfile != "None":
    popoutpath = os.path.dirname(popoutfile)
    if not os.path.exists(popoutpath):
        os.makedirs(popoutpath)
    with open(popoutfile, 'wb') as f:
        cPickle.dump(pop, f, cPickle.HIGHEST_PROTOCOL)
        

