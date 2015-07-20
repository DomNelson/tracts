# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 14:33:24 2015

@author: dominic
"""
#import matplotlib.pylab as pylab
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
labels = ['EUR', 'NAT', 'AFR']
#colordict = {'EUR':'red', 'NAT':'blue', 'AFR':'green'}
#colordict = {'EUR':'yellow', 'NAT':'green'}
colordict = {0:'red', 1:'blue', 2:'green'}
#colordict = {'AFR':'red', 'EUR':'blue'}#, 2:'green'}



#ChromLengths = [2.865747830, 2.64751457082595, 2.23363180733515, 
#                2.15492839808593, 2.04089356863902, 1.92039918028429, 
#                1.87852676459211, 1.68003441747308, 1.78206001355185, 
#                1.81366917101923, 1.58218649890248, 1.74679023161126,
#                1.26778791112187, 1.20202583329567, 1.39297570875973, 
#                1.340377262456, 1.2849052927734, 1.17708922675517, 
#                1.07733846085975, 1.08266933913055, 0.627864782064372, 
#                0.741095623349923]
                

#cM_ChromLengths = [ 277.6846783 ,  263.4266571 ,  224.5261258 ,  212.8558223 ,
#                203.9634184 ,  192.9822446 ,  186.9212679 ,  170.2156421 ,
#                168.2431216 ,  179.0947462 ,  159.5132079 ,  172.8505693 ,
#                126.9025447 ,  116.3957107 ,  131.405539  ,  134.9600594 ,
#                129.2943145 ,  119.0324459 ,  107.8670432 ,  108.0521931 ,
#                61.46827149,   72.68689882]
#
#ChromLengths = [length / 100. for length in cM_ChromLengths]

ChromLengths = [30]

migfile = sys.argv[1]
if migfile != "None":
    migmat = np.genfromtxt(migfile)
else:
    migmat = None
pedfile = sys.argv[2]
if pedfile == "None":
    pedfile = None
ancfile = sys.argv[3]
if ancfile == "None":
    ancfile = None
numinds = int(sys.argv[4])
method = sys.argv[5]
outdir = os.path.expanduser(sys.argv[6])
if not os.path.exists(outdir):
    print "Output path does not exist"
    sys.exit()
bed_dir = os.path.join(outdir + 'BED/')
if not os.path.exists(bed_dir):
    os.makedirs(bed_dir)
popname = sys.argv[7]
#try:
#    popoutfile = sys.argv[7]
#    plotoutfile = sys.argv[8]
#except IndexError:
#    popoutfile = "None"
#    plotoutfile = "None"

indlist = []
for i in range(numinds):
    print "Simulating individual", i, "of", numinds                     
    if method == "forward":
        P = ped.Pedigree(sampleind = None, 
                 DemeSwitch = DemeSwitch,
                 MigPropMat = migmat,
                 pedfile = pedfile,
                 ancfile = ancfile,
                 labels = labels,
                 split_parents = False)
        leaflist, nodelist = P.SortLeafNode(P.indlist)
#        TMat = P.BuildTransMatrices(leaflist, nodelist)
        P.MakeGenomes(ChromLengths = ChromLengths, smoothed = True,
                             Gamete = False)
        ##@@ This could be slow in a large pedigree
        samp_ind = [ind for ind in P.indlist if ind.depth == 0]
        if len(samp_ind) == 1:
            samp_ind = samp_ind[0]
        else:
            print "Depth error: mutiple roots to the pedigree"
            print samp_ind
            sys.exit()
        tracts_ind = samp_ind.to_tracts_indiv()
    ## We split the pedigree into maternal/paternal sides when simulating with
    ## PSMC
    elif method == "PSMC":
        P = ped.Pedigree(sampleind = None, 
                 DemeSwitch = DemeSwitch,
                 MigPropMat = migmat,
                 pedfile = pedfile,
                 ancfile = ancfile,
                 labels = labels,
                 split_parents = True)
        M_leaflist, M_nodelist = P.SortLeafNode(P.mother_indlist)
        F_leaflist, F_nodelist = P.SortLeafNode(P.father_indlist)
        M_TMat = P.BuildTransMatrices(M_leaflist, M_nodelist)
        F_TMat = P.BuildTransMatrices(F_leaflist, F_nodelist)
        
        tracts_ind = P.PSMC_ind(M_TMat, F_TMat, M_leaflist, F_leaflist, ChromLengths)
    else:
        print "Unknown simulation method"
        sys.exit
    ## Save simulated individual to list
    indlist.append(tracts_ind)
    ## Write simulated individuals to BED files
    if bed_dir != "None":
        outfile = os.path.join(bed_dir, "IND" + str(i + 1))
        ped.tracts_ind_to_bed(tracts_ind, outfile, conv = "M->cM")

## Plot tracts distribution for simulated population
pop = tracts.population(list_indivs = indlist)
(bins, data) = pop.get_global_tractlengths(npts=50)
#outdir = "./out"
if migmat is None:
    migmat, ancestries = P.ped_to_migmat(P.indlist)
D = tracts.demographic_model(mig=migmat)

with open(outdir + popname + "_bins", 'w') as fbins:
    fbins.write("\t".join(map(str, bins)))

with open(outdir + popname + "_dat", 'w') as fdat:
    for label in data.keys():
        fdat.write("\t".join(map(str, data[label])) + "\n")

with open(outdir + popname + "_mig", 'w') as fmig:
    for line in D.mig:
        fmig.write("\t".join(map(str, line)) + "\n")

with open(outdir + popname + "_pred", 'w') as fpred:
    for popnum in range(len(data)):
        fpred.write(
            "\t".join(map(
                str,
                pop.nind * np.array(D.expectperbin(ChromLengths, popnum, bins))))
            + "\n")
            
            
#if plotoutfile != "None":
#plotoutfile = os.path.join(outdir, plotoutfile)
#plotoutfile = os.path.join(outdir, popname, '_plot.png')            
#pop.plot_global_tractlengths(colordict, outfile = plotoutfile)

## Option to write population instance to file
#if popoutfile != "None":
#    os.path.join(outdir, popoutfile)
#    popoutpath = os.path.dirname(popoutfile)
#    if not os.path.exists(popoutpath):
#        os.makedirs(popoutpath)
#    with open(popoutfile, 'wb') as f:
#        cPickle.dump(pop, f, cPickle.HIGHEST_PROTOCOL)
        

