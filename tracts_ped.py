# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 15:32:35 2014

@author: Dom
"""

import numpy as np
import scipy.sparse
import PIL.ImageDraw as ImageDraw, PIL.Image as Image, PIL.ImageFont as ImageFont
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
    def __init__(self, MigPropMat, sampleind=None, DemeSwitch=0, labels=None,
                    split_parents = False):
        if sampleind is None:
            self.sampleind = indiv()
        else:
            self.sampleind = sampleind
        self.indlist = [self.sampleind]
        self.nextgenlist = []
        self.DemeSwitch = DemeSwitch
        ## Subtract one from length of migration matrix because it
        ## includes the zeroth generation, ie the sampled individual
        self.Gens = len(MigPropMat) - 1
        self.MigPropMat = MigPropMat
        self.currentgenlist = self.indlist
        if labels is None:
            self.labels = range(len(MigPropMat[0]))
        else:
            self.labels = labels
            
        self.sampleind.ancestry = self.SetAncestry_matrix(0)
        if self.sampleind.ancestry is not None:
            print "The proband is an ancestor"
            sys.exit()
        
        ## If true, we build separate pedigrees on the maternal and paternal
        ## sides.
        if split_parents is True:
            print "Building separate maternal/paternal pedigrees"
            sampleparents = self.set_parents(self.sampleind)
            ## Remove descendants of parents so they become the root of their
            ## respective pedigrees.
            for parent in sampleparents:
                parent.descendants = []
                parent.position = []
            ## We set the individual and their parents, so the remaining pedigrees
            ## are 2 generations shorter than the migration matrix
            self.mother_indlist = self.build_ped(len(self.MigPropMat) - 2,
                                                 sampleind = sampleparents[0])
            self.father_indlist = self.build_ped(len(self.MigPropMat) - 2,
                                                 sampleind = sampleparents[1])
            self.indlist = [self.sampleind] + self.mother_indlist + self.father_indlist
            ## Now adjsut the depth of all individuals to account for split
            for ind in self.indlist:
                ind.depth -= 1

        else:                                             
            self.indlist = self.build_ped(len(self.MigPropMat) - 1,
                                            sampleind = self.sampleind)
        
        
    def set_parents(self, ind):
        ## If not a migrant, check for deme switch - affects parents
        ##@@ Poorly implemented, only works for two demes
        demecheck = np.random.random()
        if demecheck < self.DemeSwitch:
            ind.parentdeme = (ind.deme + 1) % 2
        ## Create mother
        mother_ancestry = self.SetAncestry_matrix(ind.depth + 1)
        ind_mother = indiv(depth = ind.depth + 1, 
                               deme = ind.parentdeme, 
                               position = ind.position + [0], 
                               ancestry = mother_ancestry)
        ind.parentlist.append(ind_mother)
        ## Create father
        father_ancestry = self.SetAncestry_matrix(ind.depth + 1)
        ind_father = indiv(depth = ind.depth + 1, 
                               deme = ind.parentdeme, 
                               position = ind.position + [1], 
                               ancestry = father_ancestry)
        ind.parentlist.append(ind_father)
        ##@@ Track descendents for use in transition matrices. Will
        ## eventually have to do this with separate function so
        ## it can handle arbitrary pedigrees.
        ind_mother.descendants += ind.descendants
        ind_mother.descendants.append(ind)
        ind_father.descendants += ind.descendants
        ind_father.descendants.append(ind)
        
        return [ind_mother, ind_father]

            
    def build_ped(self, gens, depth = 0, sampleind = None, DemeSwitch = 0.1):
        if sampleind is None:
            sampleind = indiv(depth = depth)
        else:
            sampleind = sampleind
        indlist = [sampleind]
        currentgenlist = indlist
        nextgenlist = []
        for i in range(gens):
            ## Create parents for most recently created generation
            for ind in currentgenlist:
                if ind.ancestry is None:
                    nextgenlist.extend(self.set_parents(ind))
            indlist += nextgenlist
            currentgenlist = nextgenlist
            nextgenlist = []
        return indlist
                            

    def SetAncestry_matrix(self, depth, AllAncestry = True):
        if depth > len(self.MigPropMat):
            print "Error: No migrant proportions for generation", depth
            sys.exit()
        if depth == self.Gens and AllAncestry is True:
            sourceindex = weighted_choice(self.MigPropMat[depth])
            return self.labels[sourceindex]
        if np.random.random() < np.sum(self.MigPropMat[depth]):
            migindex = weighted_choice(self.MigPropMat[depth])
            return self.labels[migindex]
        else:
            return None
            
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
                sys.exit()
            numancs = len(row) - 1
        for row in new_migmat:
            if len(row) == 0:
                for i in range(numancs):
                    row.append(0)
        return new_migmat
            

    ## Sorts individuals into lists of leaves and nodes
    def SortLeafNode(self, indlist):
        leaflist = []
        nodelist = []
        
        ## Individuals are leaves if they are from the oldest generation, or 
        ## if they are a migrant from an earlier generation
        for ind in indlist:
            ##@@ Does depth check still work for sub-pedigrees?
            if ind.ancestry != None or ind.depth == self.Gens:
                leaflist.append(ind)
            else:
                nodelist.append(ind)
        
        for i in range(len(leaflist)):
            leaflist[i].label = i
        for i in range(len(nodelist)):
            nodelist[i].label = i

        return leaflist, nodelist

    def BuildTransMatrices(self, leaflist, nodelist, rho = 1):
        ## Create 'biggest possible' matrices, even though some branches may
        ## terminate early at recent migrants. Two rows for each possible node
        ## (maternal and paternal side), and one for each possible leaf.
        ##@@ Should clean up empty rows/columns, or define as sparse matrix
        ##@@ Use sparse matrices properly for performance improvements
        self.SparseLtN = scipy.sparse.csr_matrix((2 * len(nodelist), len(leaflist)))
        self.SparseNtL = scipy.sparse.csr_matrix((2 * len(nodelist), len(leaflist)))
        LtN = np.zeros((2 * len(nodelist), len(leaflist)))
        NtL = np.zeros((2 * len(nodelist), len(leaflist)))        
        
        for i in range(len(leaflist)):
            ## Transition probability depends on leaf depth, which is also 
            ## equal to its number of descendants
            TransitionProb = 1. / leaflist[i].depth
            for descendant in leaflist[i].descendants:
                ## Check if we are coming from maternal (0) or paternal (1) 
                ## side. Rows are determined by adding the individuals binary
                ## position (ie 1, 0, 1 = 5) to the maximum possible number of
                ## individuals from all previous generations.
                if leaflist[i].position[descendant.depth] == 0:
                    ##@@ Why? print "Leaflist", i, ", Descendant depth", descendant.depth
                    LtNRow = 2 * (descendant.label)
                elif leaflist[i].position[descendant.depth] == 1:
                    LtNRow = 2 * (descendant.label) + 1
                ##@@ Why? print 2 * (2 ** descendant.depth - 1 + descendant.PosAsBinary())
                try:
                    LtN[LtNRow][i] = TransitionProb
                except IndexError:
                    print "Out of bounds: Row", i, "does not exist"
        self.SparseLtN = scipy.sparse.csr_matrix(LtN)
        
        for i in range(len(leaflist)):
            for descendant in leaflist[i].descendants:
                ## Check if we are coming from maternal (0) or paternal (1) 
                ## side. Rows are determined by adding the individuals binary
                ## position (ie 1, 0, 1 = 5) to the maximum possible number of
                ## individuals from all previous generations. Select opposite
                ## maternal/paternal row compared to LtN above.
                ##@@ This could be made more robust
                if leaflist[i].position[descendant.depth] == 0:
                    NtLRow = 2 * (descendant.label) + 1
                elif leaflist[i].position[descendant.depth] == 1:
                    NtLRow = 2 * (descendant.label)
                NtL[NtLRow][i] = 1. / 2 ** (leaflist[i].depth - descendant.depth - 1)
        self.SparseNtL = scipy.sparse.csr_matrix(NtL)
        TMat = np.dot(np.transpose(LtN), NtL)
        ##@@ Think about whether this is a safe threshold
        if np.abs(np.sum(TMat[0]) - 1) > 1e-6:
            print "Transition probabilities do not sum to one. Matrix may be transposed"
            print TMat
            print LtN
            print NtL
            sys.exit()

        return TMat, LtN, NtL

    # def build_tmat(self):
    #     for ind in self.indlist:


    
    def MakeGenomes(self, TMat, ChromLengths, smoothed = True, Gamete = False):
        self.Gamete = Gamete
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
                    ind.chromosomes[i] = tracts.chropair([chrom0, chrom1])
                    # ind.chromosomes['F' + str(i)] = chrom1

        for i in range(1, self.Gens + 1):
            for ind in self.indlist:
                if ind.depth == self.Gens - i and ind.ancestry is None:
                    for k in range(len(ChromLengths)):
                        try:
                            # m0 = ind.parentlist[0].chromosomes['M' + str(k)]
                            # m1 = ind.parentlist[0].chromosomes['F' + str(k)]
                            # mom_chroms = tracts.chropair([m0, m1])
                            mom_chroms = ind.parentlist[0].chromosomes[k]
                            chrom0 = mom_chroms.recombine()
                            chrom0._smooth()
                            
                            # f0 = ind.parentlist[1].chromosomes['M' + str(k)]
                            # f1 = ind.parentlist[1].chromosomes['F' + str(k)]
                            # dad_chroms = tracts.chropair([f0, f1])
                            dad_chroms = ind.parentlist[1].chromosomes[k]
                            chrom1 = dad_chroms.recombine()
                            chrom1._smooth()
                            
                            ind.chromosomes[k] = tracts.chropair([chrom0, chrom1])
                            # ind.chromosomes['M' + str(k)] = chrom0
                            # ind.chromosomes['F' + str(k)] = chrom1
                        except KeyError:
                            print "Chromosome not found: depth", self.Gens - i
                            print ind
                            print ind.depth
                            print ind.label
                            sys.exit()

        for ind in self.indlist:
            if len(ind.chromosomes) != len(ChromLengths):
                print "Warning: not all chromosomes created!"
                print len(ind.chromosomes), len(ChromLengths)
                print ind.position
        if self.Gamete is True:
            hapind = indiv(depth = -1, chromosomes = None)
            for i in range(len(ChromLengths)):
                # m0 = self.indlist[0].chromosomes['M' + str(i)]
                # m1 = self.indlist[0].chromosomes['F' + str(i)]
                # chroms = tracts.chropair([m0, m1])
                chroms = self.indlist[0].chromosomes[i]
                newchrom = chroms.recombine()
                newchrom._smooth()                
                hapind.chromosomes[i] = newchrom
            ## Once done, add the haploid individual to the beginning of
            ## the indlist
            self.indlist = [hapind] + self.indlist

    def PSMC_chromosome(self, TMat, leaflist, chromlength, rho):
        tractlist = []
        leafindex = np.random.randint(len(leaflist))
        leaf = leaflist[leafindex]
        startpnt = 0
#        Lambda = (1. + rho * chromlength * (leaf.depth - 1)) / chromlength
        Lambda = rho * (leaf.depth)
        endpnt = np.random.exponential(1. / Lambda)
        ## Fill in all tracts up to the last one, which is done after
        while endpnt < chromlength:
            tractlist.append(tracts.tract(startpnt, endpnt, leaf.ancestry))
            ## Now build the next tract
            startpnt = endpnt
#            Lambda = (1. + rho * chromlength * (leaf.depth - 1)) / chromlength
            Lambda = rho * (leaf.depth)
            endpnt = endpnt + np.random.exponential(1. / Lambda)
            transprobs = TMat[leafindex]
            ##@@ This could be sped up using weighted_choice function
            leafindex = np.random.choice(range(len(leaflist)), p=transprobs)
            leaf = leaflist[leafindex]
        ## Fill in the last tract and build chromosome
        tractlist.append(tracts.tract(startpnt, chromlength, leaf.ancestry))
        chrom = tracts.chrom(tracts = tractlist)
        chrom._smooth()
        
        return chrom
        
    def PSMC_chropair(self, TMat, leaflist, chromlength, rho):
        chroms = []
        for i in range(2):
            chroms.append(self.PSMC_chromosome(chromlength, TMat, leaflist, rho))
        chropair = tracts.chropair(chroms)
        
        return chropair        

    def PSMC_ind(self, M_TMat, F_TMat, M_leaflist, F_leaflist, 
                chromlengths, rho=1.):
        indiv = tracts.indiv(Ls = chromlengths, label = "None")
        chropairs = []
        for length in chromlengths:
            mother_chrom = self.PSMC_chromosome(M_TMat, M_leaflist, length, rho)
            father_chrom = self.PSMC_chromosome(F_TMat, F_leaflist, length, rho)
            chropairs.append(tracts.chropair([mother_chrom, father_chrom]))
        indiv.chroms = chropairs
        
        return indiv
        
    def PlotPedigree(self, indlabels = True, tractlabels = False, 
                     showgamete = False, outfile = None):
        ##@@ Gamete is created in plot function... this shouldn't be the case
        if self.Gens >= 9:
            proceed = raw_input("Warning: Can't plot 9 or more generations \
                                clearly. Continue? [y/n]")
            if proceed != "y":
                return
        x_range = (40 * len(self.leaflist) + 5) * self.Gens
        # if self.Gamete is True and showgamete is True:
        if showgamete is True:
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
            # if self.Gamete is True and showgamete is True:
            if showgamete is True:
                ystart = (ind.depth + 1) * 300
            else:
                ystart = ind.depth * 300
            for i in range(len(ind.chromosomes)):
                for k in range(2):
                    for tract in ind.chromosomes[i].chroms[k].tracts:
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
            self.chromosomes = OrderedDict()
        else:
            if type(chromosomes) == list:
                raise
            else:
                self.chromosomes = chromosomes
        ## We can't have more than two parents
        assert len(self.parentlist) <= 2
        
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
        chromlengths = [chropair.applychrom(tracts.chrom.get_len)[0] for chropair in self.chromosomes.values()]
        tractsind = tracts.indiv(Ls = chromlengths, label = self.ancestry)
        tractsind.chroms = self.chromosomes.values()
        return tractsind
        
    def to_bed_file(self, outfile):
        header = "Chr\tStart(bp)\tEnd(bp)\tVit\tStart(cM)\tEnd(cM)\tFbk\tConfidence\tPosterior\n"
        for phase in range(2):
            if phase == 0:
                newoutfile = outfile + "_A.bed"
            elif phase == 1:
                newoutfile = outfile + "_B.bed"
            with open(newoutfile, 'w') as f:
                f.write(header)            

            for i in range(len(self.chromosomes)):
                with open(newoutfile, 'a') as f:
                    for tract in self.chromosomes[i].copies[phase].tracts:
                        line = ""
                        line += str(i + 1) + "\t"
                        line += str(0) + "\t"
                        line += str(0) + "\t"
                        line += str(tract.label) + "\t"
                        line += str(tract.start) + "\t"
                        line += str(tract.end) + "\t"
                        line += "\n"
                        
                        f.write(line)

def tracts_ind_to_bed(ind, outfile, conv = None):
    header = "Chr\tStart(bp)\tEnd(bp)\tVit\tStart(cM)\tEnd(cM)\tFbk\tConfidence\tPosterior\n"
    for phase in range(2):
        if phase == 0:
            newoutfile = outfile + "_A.bed"
        elif phase == 1:
            newoutfile = outfile + "_B.bed"
        with open(newoutfile, 'w') as f:
            f.write(header)

        if conv == "M->cM":
            conv_fact = 100
        elif conv == "cM->M":
            conv_fact = 0.01
        else:
            conv_fact = 1

        for i in range(len(ind.chroms)):
            with open(newoutfile, 'a') as f:
                for tract in ind.chroms[i].copies[phase].tracts:
                    line = ""
                    line += str(i + 1) + "\t"
                    line += str(0) + "\t"
                    line += str(0) + "\t"
                    line += str(tract.label) + "\t"
                    line += str(tract.start * conv_fact) + "\t"
                    line += str(tract.end * conv_fact) + "\t"
                    line += "\n"
                    
                    f.write(line)    

        
#migmat = [[0,0], [0,0],[0,0], [0,0], [0,0],[0,0],[0,0], [0,0],[0.5,0], [0.5,0.5]]
##migmat = [[0,0], [0,0],[0,0],[0.5,0.5]]
##migmat = [[0,0], [0.5,0.5]]
#                
#P = Pedigree(migmat)
#P.SortLeafNode()
#P.BuildTransMatrices()
#
##print P.TMat
#
#ChromLengths = [2.865747830, 2.64751457082595, 2.23363180733515, 
#                2.15492839808593, 2.04089356863902, 1.92039918028429, 
#                1.87852676459211, 1.68003441747308, 1.78206001355185, 
#                1.81366917101923, 1.58218649890248, 1.74679023161126,
#                1.26778791112187, 1.20202583329567, 1.39297570875973, 
#                1.340377262456, 1.2849052927734, 1.17708922675517, 
#                1.07733846085975, 1.08266933913055, 0.627864782064372, 
#                0.741095623349923]
#
#a = P.PSMC_ind(ChromLengths)
#for tract in a.chroms[0].copies[0].tracts:
#    print tract.start, tract.end, tract.label
                
                
                
                