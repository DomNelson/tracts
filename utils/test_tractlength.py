# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 11:37:28 2015

@author: dominic
"""
import numpy as np

def PSMC_chromosome(chromlength, rho, depth):
    tractlist = []
    startpnt = 0
#    Lambda = (1. + rho * chromlength * (depth - 1)) / chromlength
#    Lambda = chromlength / (1. + rho * chromlength * (depth - 1))
    Lambda = rho * (depth - 1)
    endpnt = np.random.exponential(1. / Lambda)
    ## Fill in all tracts up to the last one, which is done after
    while endpnt < chromlength:
        tractlist.append([startpnt, endpnt])
        ## Now build the next tract
        startpnt = endpnt
#        Lambda = (1. + rho * chromlength * (depth - 1)) / chromlength
#        Lambda = chromlength / (1. + rho * chromlength * (depth - 1))
        Lambda = rho * (depth - 1)
        endpnt = endpnt + np.random.exponential(1. / Lambda)
    ## Fill in the last tract and build chromosome
    tractlist.append([startpnt, chromlength])
    
    return tractlist, len(tractlist)
    
explens = []
poissonlens = []
generation = 6
length = 3
rho = 1
for i in range(10000):
    explens.append(PSMC_chromosome(length, 1, generation)[1])
    poissonlens.append(1 + np.sum([np.random.poisson(length) for i in range(generation - 1)]))

print 1. /( (1. + rho * chromlength * (generation - 1)) / length)
print 1. / (rho * (generation - 1))
print "Ancestor depth:", generation
print "Exponential tracts:", np.mean(explens), np.var(explens)
print "Poisson tracts:", np.mean(poissonlens), np.var(poissonlens)