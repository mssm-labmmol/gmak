from math import sqrt 
import runcmd
from random import choice
import time
import os
import re
from vbga import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Cartesian Position
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class CartesianPosition (object):

    def __init__ (self, r):
        self.dim = len(r)
        self.r   = r

    def __eq__ (self, other):
        return ((self.dim == other.dim) and (self.r == other.r))

    def __repr__ (self):
        out_str = "["
        for i in range(self.dim):
            out_str += "%10.3e" % self.r[i]
        out_str += "]"
        return out_str 

    def __getitem__ (self,i):
        return self.r[i]

    def distance_to (self, other):
        out_sum = 0.0
        for i in range(self.dim):
            out_sum += (self.r[i] - other.r[i]) ** 2
        return sqrt(out_sum)

    def shift (self,other):
        for i in range(self.dim):
            self.r[i] += other.r[i]

    def travel_to_another_dimension (self):
        self.r.append(0)
        self.dim += 1
        
    def dimensional_stack (self, other):
        output = CartesianPosition (self.r + other.r)
        return output

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# VolumeElement
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class VolumeElement (object):

    def __init__ (self, position, volume):
        self.position = position
        self.volume   = volume
        self.dimension = position.dim

    def __eq__ (self,other):
        return ((self.volume == other.volume) and (self.position == other.position))

    def __repr__ (self):
        return self.position.__repr__()

    def __getitem__ (self,i):
        return self.position[i]

    def get_volume (self):
        return self.volume

    def distance_to (self, other):
        return self.position.distance_to(other.position)

    def is_within_distance_of (self, other, dist):
        return (self.distance_to(other) <= dist)

    def shift (self, displacement):
        self.position.shift(displacement)

    def travel_to_another_dimension (self, lengthInNewDirection):
        self.position.travel_to_another_dimension()
        self.volume = self.volume * lengthInNewDirection

    def dimensional_stack (self, other):
        newElement = VolumeElement (self.position.dimensional_stack(other.position), self.volume*other.volume)
        return newElement

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# VolumeRegion
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class VolumeRegion (object):
    
    def __init__ (self, dimension):
        self.elements = []
        self.dimension = dimension

    def __getitem__ (self,i):
        return self.elements[i]

    def __add__ (self, other):
        newRegion = VolumeRegion(self.dimension)
        newRegion.copy_from(self)
        for element in other.elements:
            newRegion.add_element(element)
        return newRegion

    def __repr__ (self):
        out_str = "[ region \n"
        for el in self.elements:
            out_str += "\t" + el.__repr__()  + "\n"
        out_str += "]\n"
        return out_str

    def copy_from (self,other):
        self.elements = [el for el in other.elements]
        self.dimension = other.dimension

    def add_element (self, element):
        if element not in self.elements:
            if (element.dimension == self.dimension):
                self.elements.append(element)
            else:
                print ("Error: element and region dimensions don't match.")
                exit()

    def get_volume (self):
        out = 0
        for element in self.elements:
            out += element.get_volume()
        return out

    def get_extension (self):
        listOfX = [r[0] for r in self.elements]
        listOfY = [r[1] for r in self.elements]
        sizeX = max(listOfX) - min(listOfX)
        sizeY = max(listOfY) - min(listOfY)
        return sizeX*sizeY

    def distance_to_region (self, other):
        mindist = -1
        for el_i in self.elements:
            for el_j in other.elements:
                thisdist = el_i.distance_to(el_j)
                if (mindist == -1):
                    mindist = thisdist
                else:
                    if (thisdist < mindist):
                        mindist = thisdist
                    if (mindist == 0):
                        return mindist
        return mindist

    def clone_radius_around_element (self, radius, i):
        newRegion = VolumeRegion (self.dimension)
        for el in self.elements:
            if (el.is_within_distance_of(self.elements[i], radius)):
                newRegion.add_element(el)
        return newRegion

    def shift (self, displacement):
        for el in self.elements:
            el.shift(displacement)

    def travel_to_another_dimension (self, lengthInNewDirection):
        for el in self.elements:
            el.travel_to_another_dimension(lengthInNewDirection)

    def cartesian_product (self, other):
        newRegion = VolumeRegion (self.dimension + other.dimension)
        for el_i in self.elements:
            for el_j in other.elements:
                newRegion.add_element(el_i.dimensional_stack(el_j))
        return newRegion

    def size (self):
        return len(self.elements)

    def has_element (self, el):
        for myel in self.elements:
            if (el == myel):
                return True
        return False
        

# CubicgridRegion: a special type of VolumeRegion
class CubicgridRegion (VolumeRegion,object):

    def __init__ (self, nBoxesAtSide, boxLength, dimension):
        if (dimension == 1):
            super(CubicgridRegion,self).__init__(1)
            for i in range(nBoxesAtSide):
                self.add_element(VolumeElement(CartesianPosition([i*boxLength]), boxLength))
        else:
            self.sideLength = nBoxesAtSide * boxLength
            self.__init__(nBoxesAtSide, boxLength, 1)
            product = self.cartesian_product(CubicgridRegion(nBoxesAtSide, boxLength, dimension - 1))
            # now translate to self, because cartesian_product returns a new element
            self.copy_from(product)

    def get_size_len (self):
        return self.sideLength


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# VolumePrinter (mainly for debugging)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class TikzInterface (object):

    def __init__ (self, scaleX, scaleY):
        self.scale_x = scaleX
        self.scale_y = scaleY

    def start (self, stream):
        stream.write("\\documentclass{standalone}\n")
        stream.write("\\usepackage{tikz}\n")
        stream.write("\\begin{document}\n")
        stream.write("\\begin{tikzpicture}[xscale=%.f,yscale=%.f]\n" % (self.scale_x,self.scale_y))

    def end (self, stream):
        stream.write("\\end{tikzpicture}\n")
        stream.write("\\end{document}\n")

    def push_label_to (self, stream, label, position, color):
        stream.write("\\path [%s] (%f,%f) node {%s};\n" % (color,position[0],position[1],label))

    def push_element_to (self, stream, element, lineWidth, lineColor, fillColor, fillOpacity, gradOpacity):
        if (element.dimension == 2):
            xstart = element[0]
            ystart = element[1]
            xend   = xstart + sqrt(element.get_volume())
            yend   = ystart + sqrt(element.get_volume())
            stream.write ("\\filldraw [fill=%s!%d!white,fill opacity=%f,draw=%s] (%f,%f) rectangle (%f,%f);\n" % (fillColor,gradOpacity,fillOpacity/100.0,lineColor,xstart,ystart,xend,yend))
        else:
            print ("TikzInterface was designed to deal with 2 dimensions!\n")
            exit()

    def compile_command (self):
        return "pdflatex"

class VolumePrinter (object):

    def __init__ (self, interface):
        self.fn = "/tmp/volumeprintertmpfile"
        self.stream = open(self.fn, "w")
        self.interface = interface
        self.interface.start(self.stream)

    def drawLabel (self, label, position, color):
        self.interface.push_label_to(self.stream, label, position, color)

    def drawElement (self, element, lineWidth, lineColor, fillColor, fillOpacity, gradOpacity):
        self.interface.push_element_to(self.stream, element, lineWidth, lineColor, fillColor, fillOpacity, gradOpacity)

    def drawRegion (self, region, lineWidth, lineColor, fillColor, fillOpacity, gradOpacity):
        for el in region.elements:
            self.drawElement(el, lineWidth, lineColor, fillColor, fillOpacity, gradOpacity)

    def savefig (self, output):
        self.interface.end(self.stream)
        self.stream.close()
        runcmd.run(self.interface.compile_command() + " " + self.fn + "> /dev/null 2>&1")
        runcmd.run("mv volumeprintertmpfile.pdf " + output)
        runcmd.run("rm volumeprintertmpfile.*")


class coverageIndividualBase (object):

    def __init__ (self, size, sideLength, dimension, samples, influenceRadius, hardness):

        self.influenceRadius = influenceRadius
        self.size = size
        self.basicGrid = CubicgridRegion (size, sideLength, dimension)
        self.samplesList = samples
        self.hardness = hardness
        self.regions = []
        for sample in samples:
            influenceRegion = self.basicGrid.clone_radius_around_element(influenceRadius,sample)
            self.regions.append(influenceRegion)

    def addSample (self, sample):
        self.samplesList.append(sample)
        newInfluenceRegion = self.basicGrid.clone_radius_around_element(self.influenceRadius,sample)
        self.regions.append(newInfluenceRegion)

    def getSize (self):
        return int(self.size)

class coverageIndividual (object):

    def __init__ (self, coverageBase, maxNumberOfNewSamples):
        self.volumeCalculated = -1
        self.maxNumberOfNewSamples = maxNumberOfNewSamples
        self.influenceRadius = coverageBase.influenceRadius
        self.hardness = coverageBase.hardness
        self.basicGrid = coverageBase.basicGrid
        # dictionary, keys are the numbers of the samples on which the regions are based
        self.oldRegions = {}
        self.newRegions = {}
        # define oldRegions
        for i, sampleIndex in enumerate(coverageBase.samplesList):
            self.oldRegions[sampleIndex] = coverageBase.regions[i]

    def randomize (self):
        self.newRegions = {}
        numberOfNewRegions = choice(range(0, self.maxNumberOfNewSamples + 1))
        for _ in range(numberOfNewRegions):
            newSample = -1
            while (newSample in self.oldRegions) or (newSample in self.newRegions) or (newSample == -1):
                newSample = choice(range(self.basicGrid.size()))
            self.newRegions[newSample] = self.basicGrid.clone_radius_around_element(self.influenceRadius,newSample)

    def has_element (self, element):
        for region in self.oldRegions.values():
            if region.has_element(element):
                return True
        for region in self.newRegions.values():
            if region.has_element(element):
                return True
        return False

    def get_grid_volume (self):
        volume = self.basicGrid.get_volume()
        return volume

    def get_volume (self):
        total = 0
        if (self.volumeCalculated == -1):
            for el in self.basicGrid:
                if (self.has_element(el)):
                    total += el.get_volume()
            self.volumeCalculated = total
        else:
            total = self.volumeCalculated
        return total

    def get_extension (self):
        combinedRegion = VolumeRegion (self.basicGrid.dimension)
        for region in list(self.oldRegions.values()) + list(self.newRegions.values()):
            combinedRegion += region
        return combinedRegion.get_extension()

    def get_unoccupied_volume (self):
        return (self.get_grid_volume() - self.get_volume())

    def get_unoccupied_extension (self):
        return (self.get_grid_volume() - self.get_extension())
        
    def fitness (self):
        if (self.get_new_samples() == []):
            return 5000000000000000
        idx = list(self.newRegions.keys())[0]
        energy = 0
        for idx2 in list(self.oldRegions.keys()):
            dist = self.basicGrid[idx].distance_to(self.basicGrid[idx2])
            if (dist >= 2*self.influenceRadius):
                energy += (dist - 2*self.influenceRadius)**2
            if (dist < 2*self.influenceRadius):
                energy += self.hardness*(dist - 2*self.influenceRadius)**2
        return energy

    def print_to_file (self, interface, fn):
        printer = VolumePrinter (interface)
        printer.drawRegion (self.basicGrid, 1.0, 'gray', 'white', 100, 100)
        for oldRegion in self.oldRegions:
            printer.drawRegion(self.oldRegions[oldRegion], 1.0, 'black', 'blue', 40, 100)
        for newRegion in self.newRegions:
            printer.drawRegion(self.newRegions[newRegion], 1.0, 'black', 'red', 40, 100)
        printer.savefig(fn)

    def get_new_samples (self):
        return list(self.newRegions.keys())

    def get_total_number_of_samples (self):
        return len(list(self.newRegions.keys()) + list(self.oldRegions.keys()))

    def unset_volume (self):
        self.volumeCalculated = -1


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# GA Functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def SelectionMethod(population, selectionSize=10):
    population.sort(key=lambda x: x.fitValue, reverse=False)
    return population[:selectionSize]   #return 0,1,2....selectionSize

def CrossoverMethod(father, mother, childOne, childTwo):
    # select a new region at random in father and mother and swap them
    # forming their childs
    if (father.get_new_samples() == []) or (mother.get_new_samples() == []):
        childOne = father
        childTwo = mother
        childOne.unset_volume()
        childTwo.unset_volume()
        return [childOne, childTwo]
    fatherChoice = choice(list(father.newRegions.keys()))
    motherChoice = choice(list(mother.newRegions.keys()))
    copyFatherRegion = father.newRegions[fatherChoice]
    copyMotherRegion = mother.newRegions[motherChoice]
    childOne = father
    childTwo = mother
    del(childOne.newRegions[fatherChoice])
    childOne.newRegions[motherChoice] = copyMotherRegion
    del(childTwo.newRegions[motherChoice])
    childTwo.newRegions[fatherChoice] = copyFatherRegion
    childOne.unset_volume()
    childTwo.unset_volume()
    return [childOne, childTwo]

def MutationMethod(mutant):
    choices = ['add', 'delete']
    chosen = choice(choices)
    if (chosen == 'add'):
        newSample = -1
        while (newSample in mutant.oldRegions) or (newSample in mutant.newRegions) or (newSample == -1):
            newSample = choice(range(mutant.basicGrid.size()))
        mutant.newRegions[newSample] = mutant.basicGrid.clone_radius_around_element(mutant.influenceRadius, newSample)
    if (chosen == 'delete'):
        if (mutant.get_new_samples() != []):
            deleteSample = choice(list(mutant.newRegions.keys()))
            del(mutant.newRegions[deleteSample])
    mutant.unset_volume()
    return mutant

def Fitness(individual):
    fitness = individual.fitness()
    return fitness


class coverInterface (object):

    def __init__ (self, gridSize=0, previousSamples=[], maxSamples=0, nRuns=0, volFrac=0.85, hardness=500.00, radius=4.80, popSize=200, nGens=500):
        self.maxSamples = maxSamples
        self.nRuns = nRuns
        self.volFrac = volFrac
        self.hardness = hardness
        self.radius = radius
        self.popSize = popSize
        self.nGens = nGens
        self.previousSamples = previousSamples[:]
        self.minvolFrac = 0.50
        self.noStart = False
        # base
        self.covBase = coverageIndividualBase (gridSize, 1, 2, previousSamples[:], self.radius, self.hardness)

    def clear (self):
        self.previousSamples = []

    def prepareForRun (self, samplesList, gridSize):
        self.clear()
        self.previousSamples = samplesList[:]
        self.covBase = coverageIndividualBase (gridSize, 1, 2, samplesList[:], self.radius, self.hardness)

    def run (self):
        chosenSamples = []
        bestRun = -1
        sumTimes = 0

        print ("BEGIN GACOVER RUN")
        for irun in range(self.nRuns):
            self.prepareForRun(self.previousSamples, self.covBase.getSize())
            startTime = time.time()
            bestVolFrac = -1
            thisVolFrac = -1
            numberOfSamples = len(self.previousSamples)
            addedSamples = []
            while (thisVolFrac < self.volFrac) and ((thisVolFrac < self.minvolFrac) or (numberOfSamples < self.maxSamples)):
                # Run GA and get best volFrac and numberOfNewSamples.
                myGA = VBGA(coverageIndividual, self.covBase, 1, Fitness, SelectionMethod, CrossoverMethod, 10, MutationMethod, 5, self.popSize)
                myGA.run(self.nGens)
                # Get best individual.
                optind = myGA.getBest()
                newSamples = optind.get_new_samples()
                numberOfSamples = optind.get_total_number_of_samples()
                numberOfNewSamples = len(newSamples)
                thisVolFrac = float(optind.get_volume()) / float(optind.get_grid_volume())
                if numberOfNewSamples > 0:
                    addedSamples += newSamples
                    self.covBase.addSample (newSamples[0])
            thisTime = (time.time() - startTime)
            sumTimes += thisTime
            if (thisVolFrac > bestVolFrac):
                bestVolFrac = thisVolFrac
                chosenSamples = addedSamples[:]
                bestRun = irun
            try:
                tikz = TikzInterface (1.0, 1.0)
                optind.print_to_file(tikz, "gacover_%d.pdf" % (irun+1))
            except:
                print("Warning: An exception occurred while printing a cover grid to a file.")
            print ("\t---> GACOVER RUN %d --------> Time spent: %d s" % (irun+1, int(thisTime)))
            print ("\t                    --------> Volume fraction covered: %f" % (thisVolFrac))
            print ("\t                    --------> Optimal samples: " + str(addedSamples))
        print ("Total time spent: %d s" % sumTimes)
        print ("Best volume fraction covered: %f" % bestVolFrac)
        print ("Best optimal samples: " + str(chosenSamples))
        print ("END GACOVER RUN")

        return chosenSamples

    # stream is at a '$gacover' line
    def readFromStream (self, stream):
        for line in stream:
            if line[0] == '#':
                continue
            if (re.match(r"^\$end.*",line)):
                return
            splitted = line.split()
            if (splitted[0] == 'maxsamples'):
                self.maxSamples = int(splitted[1])
            if (splitted[0] == 'nruns'):
                self.nRuns = int(splitted[1])
            if (splitted[0] == 'volfrac'):
                self.volFrac = float(splitted[1])
            if (splitted[0] == 'radius'):
                self.radius = float(splitted[1])
            if (splitted[0] == 'popsize'):
                self.popSize = int(splitted[1])
            if (splitted[0] == 'ngens'):
                self.nGens = int(splitted[1])
            if (splitted[0] == 'nostart'):
                self.noStart = True
        
if __name__ == "__main__":
    samples = [544]
    covInt = coverInterface (33, previousSamples=samples, maxSamples=15, nRuns=10, volFrac=0.85, hardness=500.00, radius=4.80, popSize=100, nGens=20)
    covInt.run()
