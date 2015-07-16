# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:15:29 2015

@author: dominic
"""
import numpy as np
import sys, os
sys.path.append(os.path.abspath('../'))
import tracts_ped as ped
import tracts

labels = ['EUR', 'NAT', 'AFR']
#colordict = {'EUR':'red', 'NAT':'blue', 'AFR':'green'}
colordict = {0:'yellow', 1:'green'}
#colordict = {0:'red', 1:'blue', 2:'green'}

cM_ChromLengths = [ 277.6846783 ,  263.4266571 ,  224.5261258 ,  212.8558223 ,
                203.9634184 ,  192.9822446 ,  186.9212679 ,  170.2156421 ,
                168.2431216 ,  179.0947462 ,  159.5132079 ,  172.8505693 ,
                126.9025447 ,  116.3957107 ,  131.405539  ,  134.9600594 ,
                129.2943145 ,  119.0324459 ,  107.8670432 ,  108.0521931 ,
                61.46827149,   72.68689882]
ChromLengths = [length / 100. for length in cM_ChromLengths]
migmat = [[0,0],[0,0],[0.75, 0.25]]
#migmat = [[0,0],[0.75, 0.25]]

pedtype = sys.argv[1]
outpath = os.path.expanduser(sys.argv[2])
if not os.path.exists(outpath):
    print "Output path does not exist"
    sys.exit()
bedpath = outpath + 'BED/'
if not os.path.exists(bedpath):
    os.makedirs(bedpath)

if pedtype == '0001':
    P = ped.Pedigree(migmat, split_parents=True)
    leaflist = P.SortLeafNode(P.indlist)[0]
    
    ## If this sum is one we know that there is only one leaf with ancestry 1
    while np.sum([leaf.ancestry for leaf in leaflist]) != 1:
        P = ped.Pedigree(migmat, split_parents=True)
        leaflist = P.SortLeafNode(P.indlist)[0]
elif pedtype == '1001':
    P = ped.Pedigree(migmat, split_parents=True)
    M_leaflist, M_nodelist = P.SortLeafNode(P.mother_indlist)
    F_leaflist, F_nodelist = P.SortLeafNode(P.father_indlist)
    
    ## We want one different ancestor on each side of the pedigree
    while (np.sum([leaf.ancestry for leaf in M_leaflist]) != 1 or 
            np.sum([leaf.ancestry for leaf in F_leaflist]) != 1):
        P = ped.Pedigree(migmat, split_parents=True)
        M_leaflist, M_nodelist = P.SortLeafNode(P.mother_indlist)
        F_leaflist, F_nodelist = P.SortLeafNode(P.father_indlist)
else:
    print "Unknown pedigree type"
    sys.exit()
    

M_leaflist, M_nodelist = P.SortLeafNode(P.mother_indlist)
F_leaflist, F_nodelist = P.SortLeafNode(P.father_indlist)
M_TMat = P.BuildTransMatrices(M_leaflist, M_nodelist)[0]
F_TMat = P.BuildTransMatrices(F_leaflist, F_nodelist)[0]
    
#print P.TMat
#for leaf in P.leaflist:
#    print leaf.ancestry
# bedpath = os.path.expanduser('~/project/tracts/sims/results/PSMC/simple')


### PSMC simulations
PSMCinds = []
for i in range(1000):
    print "PSMC", i
    PSMCind = P.PSMC_ind(M_TMat, F_TMat, M_leaflist, F_leaflist, ChromLengths)
    PSMCinds.append(PSMCind)
    outfile = os.path.join(bedpath, "PSMC_IND" + str(i + 1))
    ped.tracts_ind_to_bed(PSMCind, outfile, conv = "M->cM")
    
# plotoutfile = os.path.expanduser('~/project/tracts/sims/results/PSMC/simple/PSMC0001.png')
# plotoutpath = os.path.expanduser(sys.argv[2])
if pedtype == '0001':
    plotoutfile = outpath + "PSMC0001.png"
elif pedtype == '1001':
    plotoutfile = outpath + "PSMC1001.png"
pop = tracts.population(list_indivs = PSMCinds)
pop.plot_global_tractlengths(colordict, outfile = plotoutfile)


## Forward simulations
#if pedtype == '0001':
#    P = ped.Pedigree(migmat, split_parents=False)
#    leaflist = P.SortLeafNode(P.indlist)[0]
#    
#    ## If this sum is one we know that there is only one leaf with ancestry 1
#    while np.sum([leaf.ancestry for leaf in leaflist]) != 1:
#        P = ped.Pedigree(migmat, split_parents=False)
#        leaflist = P.SortLeafNode(P.indlist)[0]
#elif pedtype == '1001':
#    P = ped.Pedigree(migmat, split_parents=False)
#    leaflist, nodelist = P.SortLeafNode(P.indlist)
#    
#    ## We want one different ancestor on each side of the pedigree
#    while (np.sum([leaf.ancestry for leaf in leaflist[:2]]) != 1 or 
#            np.sum([leaf.ancestry for leaf in leaflist[2:]]) != 1):
#        P = ped.Pedigree(migmat, split_parents=False)
#        leaflist, nodelist = P.SortLeafNode(P.indlist)
#else:
#    print "Unknown pedigree type"
#    sys.exit()
#    
#
##M_leaflist, M_nodelist = P.SortLeafNode(P.mother_indlist)
##F_leaflist, F_nodelist = P.SortLeafNode(P.father_indlist)
##leaflist, nodelist = P.SortLeafNode(P.indlist)
#TMat = P.BuildTransMatrices(leaflist, nodelist)[0]
##F_TMat = P.BuildTransMatrices(F_leaflist, F_nodelist)[0]
#
#
##P = ped.Pedigree(migmat, split_parents=False)
##leaflist, nodelist = P.SortLeafNode(P.indlist)
##
#### If this sum is one we know that there is only one leaf with ancestry 1
##while np.sum([leaf.ancestry for leaf in leaflist]) != 1:
##    P = ped.Pedigree(migmat, split_parents=False)
##    leaflist, nodelist = P.SortLeafNode(P.indlist)
##
##TMat = P.BuildTransMatrices(leaflist, nodelist)
#
#forwardinds = []
#for i in range(1000):
#   print "Forward", i
#   P.MakeGenomes(TMat, ChromLengths = ChromLengths, smoothed = True,
#                            Gamete = False)
#   forwardind = P.indlist[0].to_tracts_indiv()
#   forwardinds.append(forwardind)
#   outfile = os.path.join(bedpath, "forward_IND" + str(i + 1))
#   ped.tracts_ind_to_bed(forwardind, outfile, conv = "M->cM")
#   
## plotoutfile = os.path.expanduser('~/project/tracts/sims/results/PSMC/simple/forward0001.png')
#if pedtype == '0001':
#    plotoutfile = outpath + "forward0001.png"
#elif pedtype == '1001':
#    plotoutfile = outpath + "forward1001.png"
#pop = tracts.population(list_indivs = forwardinds)
#pop.plot_global_tractlengths(colordict, outfile = plotoutfile)