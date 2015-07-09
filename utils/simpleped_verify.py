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
#colordict = {0:'yellow', 1:'green'}
colordict = {0:'red', 1:'blue', 2:'green'}

cM_ChromLengths = [ 277.6846783 ,  263.4266571 ,  224.5261258 ,  212.8558223 ,
                203.9634184 ,  192.9822446 ,  186.9212679 ,  170.2156421 ,
                168.2431216 ,  179.0947462 ,  159.5132079 ,  172.8505693 ,
                126.9025447 ,  116.3957107 ,  131.405539  ,  134.9600594 ,
                129.2943145 ,  119.0324459 ,  107.8670432 ,  108.0521931 ,
                61.46827149,   72.68689882]
ChromLengths = [length / 100. for length in cM_ChromLengths]
migmat = [[0,0],[0,0],[0.75, 0.25]]
#migmat = [[0,0],[0.75, 0.25]]

P = ped.Pedigree(migmat)
P.SortLeafNode()
P.BuildTransMatrices()

## If this sum is one we know that there is only one leaf with ancestry 1
while np.sum([leaf.ancestry for leaf in P.leaflist]) != 1:
    P = ped.Pedigree(migmat)
    P.SortLeafNode()
    P.BuildTransMatrices()
    
#print P.TMat
#for leaf in P.leaflist:
#    print leaf.ancestry
bedpath = os.path.expanduser('~/project/tracts/sims/results/PSMC/simple')

PSMCinds = []
for i in range(1000):
    print "PSMC", i
    PSMCind = P.PSMC_ind(ChromLengths)
    PSMCinds.append(PSMCind)
    outfile = os.path.join(bedpath, "PSMC_IND" + str(i + 1))
    ped.tracts_ind_to_bed(PSMCind, outfile, conv = "M->cM")
    
plotoutfile = os.path.expanduser('~/project/tracts/sims/results/PSMC/simple/PSMC0001.png')
pop = tracts.population(list_indivs = PSMCinds)
pop.plot_global_tractlengths(colordict, outfile = plotoutfile)

#forwardinds = []
#for i in range(1000):
#    print "Forward", i
#    P.MakeGenomes(ChromLengths = ChromLengths, rho = 1, smoothed = True,
#                             Gamete = False)
#    forwardind = P.indlist[0].to_tracts_indiv()
#    forwardinds.append(forwardind)
#    outfile = os.path.join(bedpath, "forward_IND" + str(i + 1))
#    ped.tracts_ind_to_bed(forwardind, outfile, conv = "M->cM")
#    
#plotoutfile = os.path.expanduser('~/project/tracts/sims/results/PSMC/simple/forward0001.png')
#pop = tracts.population(list_indivs = forwardinds)
#pop.plot_global_tractlengths(colordict, outfile = plotoutfile)