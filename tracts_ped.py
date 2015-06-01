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
        self.MigPropMat = self.expand_migmat(MigPropMat)
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
        ## Migrant matrix uses depth - 1 since the sequenced individual
        ## at depth = 0 can't be a migrant, reducing rows by 1
        if np.random.random() < np.sum(self.MigPropMat[depth - 1]):
            migindex = weighted_choice(self.MigPropMat[depth - 1])
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
                            chrom0 = recombine(ind.parentlist[0].chromosomes['M' + str(k)], 
                                                   ind.parentlist[0].chromosomes['F' + str(k)],
                                                   smoothed)
                            chrom1 = recombine(ind.parentlist[1].chromosomes['M' + str(k)], 
                                                   ind.parentlist[1].chromosomes['F' + str(k)],
                                                   smoothed)
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
                self.gamete = recombine(self.indlist[0].chromosomes['M' + str(i)], 
                                           self.indlist[0].chromosomes['F' + str(i)],
                                           smoothed = True)
                hapind.chromosomes[str(i)] = self.gamete
            ## Once done, add the haploid individual to the beginning of
            ## the indlist
            self.indlist = [hapind] + self.indlist
            
    
    def BuildTransMatrices(self, stepsize):
        self.GlobalStepSize = stepsize
        self.SortLeafNode()
        ## Create 'biggest possible' matrices, even though some branches may
        ## terminate early at recent migrants. Two rows for each possible node
        ## (maternal and paternal side), and one for each possible leaf.
        ##@@ Should clean up empty rows/columns, or define as sparse matrix
        self.SparseLtN = sp.csr_matrix((2 * len(self.nodelist), len(self.leaflist)))
        self.SparseNtL = sp.csr_matrix((2 * len(self.nodelist), len(self.leaflist)))
        self.NonSelf = self.rho * self.GlobalStepSize * self.Gens
        self.Self = 1 - self.NonSelf
        self.LtN = np.zeros((2 * len(self.nodelist), len(self.leaflist)))
        self.NtL = np.zeros((2 * len(self.nodelist), len(self.leaflist)))        
        
        for i in range(len(self.leaflist)):
            ## Transition probability depends on leaf depth, which is also 
            ## equal to its number of descendants
            TransitionProb = 1. / self.leaflist[i].depth
            for descendant in self.leaflist[i].descendants:
                ## Check if we are coming from maternal (0) or paternal (1) 
                ## side. Rows are determined by adding the individuals binary
                ## position (ie 1, 0, 1 = 5) to the maximum possible number of
                ## individuals from all previous generations.
                if self.leaflist[i].position[descendant.depth] == 0:
                    ##@@ Why? print "Leaflist", i, ", Descendant depth", descendant.depth
                    LtNRow = 2 * (descendant.label)
                elif self.leaflist[i].position[descendant.depth] == 1:
                    LtNRow = 2 * (descendant.label) + 1
                ##@@ Why? print 2 * (2 ** descendant.depth - 1 + descendant.PosAsBinary())
                try:
                    self.LtN[LtNRow][i] = TransitionProb
                except IndexError:
                    print "Out of bounds: Row", i, "does not exist"
        self.SparseLtN = sp.csr_matrix(self.LtN)
        
        for i in range(len(self.leaflist)):
            for descendant in self.leaflist[i].descendants:
                ## Check if we are coming from maternal (0) or paternal (1) 
                ## side. Rows are determined by adding the individuals binary
                ## position (ie 1, 0, 1 = 5) to the maximum possible number of
                ## individuals from all previous generations. Select opposite
                ## maternal/paternal row compared to LtN above.
                if self.leaflist[i].position[descendant.depth] == 0:
                    NtLRow = 2 * (descendant.label) + 1
                elif self.leaflist[i].position[descendant.depth] == 1:
                    NtLRow = 2 * (descendant.label)
                self.NtL[NtLRow][i] = 1. / 2 ** (self.leaflist[i].depth - descendant.depth - 1)
        self.SparseNtL = sp.csr_matrix(self.NtL)
               
    def TransMatrixCombine(self):
        Mult = np.dot(np.transpose(self.LtN), self.NtL)

        for i in range(len(self.leaflist)):
            
            NonSelf = self.rho * self.GlobalStepSize * self.Gens
                
            assert NonSelf < 1, "Self-Transition probability > 1 %f" % (1 - NonSelf)
            assert NonSelf > 0, "Self-Transition probability < 0 %f" % (1 - NonSelf)
            
            ## Scaling of transition probabilities
            Mult[:,i] *= NonSelf
            ## Self-transition probability
            Mult[i][i] = 1 - NonSelf
        
        #print "Self-transition probability: ", 1 - NonSelf
        return Mult
        
    def BuildEmissionMatrix(self):
        ##@@ So far only 2 different ancestries/emitted symbols
        EmissionMatrix = np.zeros([len(self.leaflist), 2])
        for i in range(len(self.leaflist)):
            EmissionMatrix[i][self.leaflist[i].ancestry] = 1
        
        return EmissionMatrix
        
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


def recombine(chromosome0, chromosome1, smoothed = True):
    assert chromosome0.get_len() == chromosome1.get_len(), "Mismatched chromosomes"
    events = 0
    while events == 0:
        events = np.random.poisson(lam = chromosome0.get_len())
    #print "Number of recombination events:", events
    switchpoints = []
    switchpoints.extend(np.random.uniform(0, chromosome0.get_len(), events))
    switchpoints.sort()
    newtracts = []
    index = [0, 0] 
    chromo_pair = [chromosome0.tracts, chromosome1.tracts]
    currentchrom = np.random.random_integers(0,1)
    
    for i in range(events):
        ##@@ We should also check for the (unlikely) case when a recombination
        ##@@ lands exactly on a tract switch. I'm sure this could also be 
        ##@@ cleaned up a bit
        ## Copy full tracts from the current chromosome until one containing a 
        ## switchpoint is reached, saving the location on both chromosomes, and
        ## entering the split tract for the next step
        if i == 0:
            while chromo_pair[0][index[0]].end < switchpoints[i]:
                if currentchrom == 0:
                    newtracts.append(chromo_pair[0][index[0]])
                index[0] += 1        
            while chromo_pair[1][index[1]].end < switchpoints[i]:
                if currentchrom == 1:
                    newtracts.append(chromo_pair[1][index[1]])
                index[1] += 1
        else:
            while chromo_pair[0][index[0]].end < switchpoints[i]:
                if currentchrom == 0 and chromo_pair[0][index[0]].start > switchpoints[i - 1]:
                    newtracts.append(chromo_pair[0][index[0]])
                index[0] += 1        
            while chromo_pair[1][index[1]].end < switchpoints[i]:
                if currentchrom == 1 and chromo_pair[1][index[1]].start > switchpoints[i - 1]:
                    newtracts.append(chromo_pair[1][index[1]])
                index[1] += 1
        
        ## Now that we have found the tracts, copy up to the current switch
        ## point, from either the start of the tract or the last switch,
        ## whichever is closer.
        if i == 0:
            start = chromo_pair[currentchrom][index[currentchrom]].start
            end = switchpoints[i]
        else:
            start = np.max([chromo_pair[currentchrom][index[currentchrom]].start, 
                            switchpoints[i - 1]])
            end = switchpoints[i]
        label = chromo_pair[currentchrom][index[currentchrom]].label    
        if start > end:
            print "Bad order, ", start, end
        
        newtracts.append(tracts.tract(start, end, label))
        
        ## Switch chromosomes following recombination
        if currentchrom == 0:
            currentchrom = 1
        elif currentchrom == 1:
            currentchrom = 0
        else:
            print "oops!"        
        
        ## Need to catch the last section of a tract after the recombinations
        ## it contains, otherwise we would skip to the start of the next tract
        if i < (events - 1):
            if chromo_pair[currentchrom][index[currentchrom]].end < switchpoints[i + 1]:
                start = switchpoints[i]
                end = chromo_pair[currentchrom][index[currentchrom]].end
                label = chromo_pair[currentchrom][index[currentchrom]].label
                newtracts.append(tracts.tract(start, end, label))

    ## Add remaining partial tract following last recombination event, then add all
    ## following complete tracts.
    ##@@ Should this be switchpoints[events - 1]?
    start = switchpoints[-1]
    end = chromo_pair[currentchrom][index[currentchrom]].end
    label = chromo_pair[currentchrom][index[currentchrom]].label
    newtracts.append(tracts.tract(start, end, label))

    for i in range(index[currentchrom] + 1, len(chromo_pair[currentchrom])):
        newtracts.append(chromo_pair[currentchrom][i])
    
    newchrom = tracts.chrom(ls = chromosome0.get_len(), tracts = newtracts)
    if smoothed is True:
        newchrom._smooth()
    return newchrom
