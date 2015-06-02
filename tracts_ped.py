# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 15:32:35 2014

@author: Dom
"""

import numpy as np
import scipy.sparse as sp
import PIL.ImageDraw as ImageDraw, PIL.Image as Image, PIL.ImageFont as ImageFont
#import hmm_struct as struct
import sys
import tracts
from collections import OrderedDict

def weighted_choice(weights):
    totals = []
    running_total = 0

    for w in weights:
        running_total += w
        totals.append(running_total)

    rnd = np.random.random() * running_total
    for i, total in enumerate(totals):
        if rnd < total:
            return i

## We want to create a pedigree using the following steps:
##
## 1) create two parents for leaf node
## 2) chance to switch demes (parents are from different deme)
## 3) chance that the individual is a migrant (leaf of pedigree)

class Pedigree:
    def __init__(self, sampleind, Gens, RecentMigrants, AllAncestry,
                 DemeSwitch, SourceProps, MigPropMat):
        if sampleind is None:
            self.sampleind = indiv()
        else:
            self.sampleind = sampleind
        self.indlist = [self.sampleind]
        self.nextgenlist = []
        self.Gens = Gens
        self.RecentMigrants = RecentMigrants
        self.AllAncestry = AllAncestry
        self.SourceProps = SourceProps
        self.MigPropMat = MigPropMat#self.expand_migmat(MigPropMat)
        self.currentgenlist = self.indlist
        
        for i in range(self.Gens):
            
            ## Create parents for most recently created generation
            for ind in self.currentgenlist:
                if ind.ancestry is None:
                    
                    ## If not a migrant, check for deme switch - affects parents
                    demecheck = np.random.random()
                    if demecheck < DemeSwitch:
                        self.parentdeme = (ind.deme + 1) % 2
                    
                    ## Create mother
                    # mother_ancestry = self.SetAncestry(ind.depth + 1)
                    mother_ancestry = self.SetAncestry_matrix(ind.depth + 1)
                    ind_mother = indiv(depth = ind.depth + 1, 
                                           deme = ind.parentdeme, 
                                           position = ind.position + [0], 
                                           ancestry = mother_ancestry)
                    self.nextgenlist.append(ind_mother)
                    ind.parentlist.append(ind_mother)
                    
                    
                    ## Create father
                    # father_ancestry = self.SetAncestry(ind.depth + 1)
                    father_ancestry = self.SetAncestry_matrix(ind.depth + 1)
                    ind_father = indiv(depth = ind.depth + 1, 
                                           deme = ind.parentdeme, 
                                           position = ind.position + [1], 
                                           ancestry = father_ancestry)
                    self.nextgenlist.append(ind_father)
                    ind.parentlist.append(ind_father)
                    
                    ##@@ Track descendents for use in transition matrices. Will
                    ## eventually have to do this with separate function so
                    ## it can handle arbitrary pedigrees.
                    ind_mother.descendants += ind.descendants
                    ind_mother.descendants.append(ind)
                    ind_father.descendants += ind.descendants
                    ind_father.descendants.append(ind)
                    
            self.indlist += self.nextgenlist
            self.currentgenlist = self.nextgenlist
            self.nextgenlist = []
                            

    def SetAncestry_matrix(self, depth):
        if depth > len(self.MigPropMat):
            print "Error: No migrant proportions for generation", depth
            sys.exit()
        if np.random.random() < np.sum(self.MigPropMat[depth]):
            migindex = weighted_choice(self.MigPropMat[depth])
            return migindex
        else:
            if depth == self.Gens:
                sourceindex = weighted_choice(self.SourceProps)
                return sourceindex
            # else:
            #     return None
    def expand_migmat(self, migmat):
        if len(migmat) == self.Gens:
            print "Matrix already proper size: skipping expansion"
            return migmat

        new_migmat = [[] for gen in range(self.Gens)]
        ## Generation of migration is stored in the first entry
        ## of every row
        for row in migmat:
            try:
                ## Subtract 1 so columns in migmat are numbered
                ## from 1 not 0
                new_migmat[row[0] - 1].extend(row[1:])
            except IndexError:
                print "Migration in generation beyond pedigree"
                sys.exit
            numancs = len(row) - 1
        for row in new_migmat:
            if len(row) == 0:
                for i in range(numancs):
                    row.append(0)
        return new_migmat
            

    ## Sorts individuals into lists of leaves and nodes
    def SortLeafNode(self):
        self.leaflist = []
        self.nodelist = []
        
        ## Individuals are leaves if they are from the oldest generation, or 
        ## if they are a migrant from an earlier generation
        for ind in self.indlist:
            if ind.ancestry != None or ind.depth == self.Gens:
                self.leaflist.append(ind)
            else:
                self.nodelist.append(ind)
        
        for i in range(len(self.leaflist)):
            self.leaflist[i].label = i
        for i in range(len(self.nodelist)):
            self.nodelist[i].label = i
    
    def MakeGenomes(self, ChromLengths, rho, smoothed = True, Gamete = False):
        self.Gamete = Gamete
        self.rho = rho
        ##@@ Make sure new chromosomes plot correctly
        for ind in self.indlist:
            if ind.ancestry is not None:
                for i in range(len(ChromLengths)):
                    ## Tracks donor individual using their position in the
                    ## pedigree.
                    first_tract = tracts.tract(0, ChromLengths[i], ind.ancestry)
                    chrom0 = tracts.chrom(ls=ChromLengths[i], 
                                               tracts = [first_tract])
                    chrom1 = tracts.chrom(ls=ChromLengths[i], 
                                               tracts = [first_tract])
                    ind.chromosomes['M' + str(i)] = chrom0
                    ind.chromosomes['F' + str(i)] = chrom1

        for i in range(1, self.Gens + 1):
            for ind in self.indlist:
                if ind.depth == self.Gens - i and ind.ancestry is None:
                    for k in range(len(ChromLengths)):
                        try:
                            m0 = ind.parentlist[0].chromosomes['M' + str(k)]
                            m1 = ind.parentlist[0].chromosomes['F' + str(k)]
                            mom_chroms = tracts.chropair([m0, m1])
                            chrom0 = mom_chroms.recombine()
                            chrom0._smooth()
                            
                            f0 = ind.parentlist[1].chromosomes['M' + str(k)]
                            f1 = ind.parentlist[1].chromosomes['F' + str(k)]
                            dad_chroms = tracts.chropair([f0, f1])
                            chrom1 = dad_chroms.recombine()
                            chrom1._smooth()
                            
                            ind.chromosomes['M' + str(k)] = chrom0
                            ind.chromosomes['F' + str(k)] = chrom1
                        except KeyError:
                            print "Chromosome not found: depth", self.Gens - i
                            print ind
                            print ind.depth
                            print ind.label
                            sys.exit()

        for ind in self.indlist:
            if len(ind.chromosomes) != 2 * len(ChromLengths):
                print "Warning: not all chromosomes created!"
                print len(ind.chromosomes), len(ChromLengths)
                print ind.position
        if self.Gamete is True:
            hapind = indiv(depth = -1, chromosomes = None)
            for i in range(len(ChromLengths)):
                m0 = self.indlist[0].chromosomes['M' + str(i)]
                m1 = self.indlist[0].chromosomes['F' + str(i)]
                chroms = tracts.chropair([m0, m1])
                newchrom = chroms.recombine()
                newchrom._smooth()                
                hapind.chromosomes[str(i)] = newchrom
            ## Once done, add the haploid individual to the beginning of
            ## the indlist
            self.indlist = [hapind] + self.indlist

        
    def PlotPedigree(self, indlabels = True, tractlabels = False, 
                     showgamete = False, outfile = None):
        ##@@ Gamete is created in plot function... this shouldn't be the case
        if self.Gens >= 9:
            proceed = raw_input("Warning: Can't plot 9 or more generations \
                                clearly. Continue? [y/n]")
            if proceed != "y":
                return
        x_range = (40 * len(self.leaflist) + 5) * self.Gens
        if self.Gamete is True and showgamete is True:
            y_range = 300 * (self.Gens + 2)
        else:
            y_range = 300 * (self.Gens + 1)
        if indlabels is True or tractlabels is True:

            font = ImageFont.truetype("/Users/Dominic/project/PSMC_tracts/bin/Futura.ttc",25)
        im = Image.new("RGB", (x_range, y_range))
        draw = ImageDraw.Draw(im)
#        plotlist = self.indlist
#        if showgamete is True:
#            self.gamete = recombine(self.indlist[0].chromosomes[0], 
#                               self.indlist[0].chromosomes[1],
#                               smoothed = False)
#            hapind = indiv(depth = -1, chromosomes = [self.gamete])
#            plotlist = [hapind] + plotlist
        for ind in self.indlist:
            if self.Gamete is True and showgamete is True:
                ystart = (ind.depth + 1) * 300
            else:
                ystart = ind.depth * 300
            for i in range(len(ind.chromosomes)):
                for tract in ind.chromosomes[i].tracts:
                    if tract.label ==  0:
                        colour = 255
                    elif tract.label == 1:
                        colour = 170
                    xstart = x_range / 2
                    for j in range(len(ind.position)):
                        if ind.position == []:
                            xstart = x_range / 2
                        elif ind.position[j] == 0:
                            xstart -= x_range / (2 ** (j + 2))
                        elif ind.position[j] == 1:
                            xstart += x_range / (2 ** (j + 2))
                        else:
                            print "Invalid position"
                    draw.rectangle(((xstart + 5 + i * 15, ystart + tract.start * 100 + 5),
                                    (xstart + 15 + i * 15, ystart + tract.end * 100 + 5)),
                                    fill=colour)
                if indlabels is True:
                    if ind.depth >= 0:
                        draw.text((xstart - 30, ystart + 50), 
                                  str(ind.label), 
                                     (255, 255, 255), font = font)
                    else:
                        draw.text((xstart - 100, ystart + 50), 
                                  "Gamete", 
                                     (255, 255, 255), font = font)
        if outfile is not None:
            im.save(outfile)
        else:
            im.show()
            
            
class indiv:
    def __init__(self, depth = 0, deme = 0, position = None, ancestry = None,
                 chromosomes = None):
        self.depth = depth
        self.deme = deme
        self.ancestry = ancestry
        self.parentdeme = deme
        self.parentlist = []
        self.descendants = []
        self.label = None
        if position is None:
            self.position = []
        else:
            self.position = position
        if chromosomes is None:
            ## Two empty lists, one maternal one paternal
            self.chromosomes = OrderedDict()
        else:
            if type(chromosomes) == list:
                raise
            self.chromosomes = chromosomes
        
    ## Turns individual position to binary number then returns it as an integer
    def PosAsBinary(self):
        ## First check if individual is sample, since 'int' can't handle an
        ## empty position.
        ## Duplicate values are handled later, when creating the transition
        ## matrices.
            
        if self.position != []:
            s = ''.join(map(str, self.position))
            return int(s, 2)
        else:
            return 0
            
    def to_tracts_indiv(self):
        chromlengths = [chrom.get_len() for chrom in self.chromosomes.values()]
        tractsind = tracts.indiv(Ls = chromlengths, label = self.ancestry)
        tractsind.chroms = self.chromosomes.values()
        return tractsind