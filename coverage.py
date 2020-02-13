from math import sqrt 
from random import choice
import os

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
        os.system(self.interface.compile_command() + " " + self.fn + "> /dev/null 2>&1")
        os.system("mv volumeprintertmpfile.pdf " + output)
        os.system("rm volumeprintertmpfile.*")


class coverageIndividualBase (object):

    def __init__ (self, size, sideLength, dimension, samples, influenceRadius):

        self.influenceRadius = influenceRadius
        self.basicGrid = CubicgridRegion (size, sideLength, dimension)
        self.samplesList = samples
        self.regions = []
        for sample in samples:
            influenceRegion = self.basicGrid.clone_radius_around_element(influenceRadius,sample)
            self.regions.append(influenceRegion)

class coverageIndividual (object):

    def __init__ (self, coverageBase, numberOfNewSamples):
        self.numberOfNewSamples = numberOfNewSamples
        self.influenceRadius = coverageBase.influenceRadius
        self.basicGrid = coverageBase.basicGrid
        # dictionary, keys are the numbers of the samples on which the regions are based
        self.oldRegions = {}
        self.newRegions = {}
        # define oldRegions
        for i, sampleIndex in enumerate(coverageBase.samplesList):
            self.oldRegions[sampleIndex] = coverageBase.regions[i]

    def randomize (self):
        # clear
        self.newRegions = {}
        for _ in range(self.numberOfNewSamples):
            newSample = -1
            while (newSample in self.oldRegions) or (newSample in self.newRegions) or (newSample == -1):
                newSample = choice(range(self.basicGrid.size()))
            self.newRegions[newSample] = self.basicGrid.clone_radius_around_element(self.influenceRadius,newSample)

    def get_grid_volume (self):
        return self.basicGrid.get_volume()

    def get_volume (self):
        combinedRegion = VolumeRegion (self.basicGrid.dimension)
        for region in list(self.oldRegions.values()) + list(self.newRegions.values()):
            combinedRegion += region
        return combinedRegion.get_volume()

    def get_unoccupied_volume (self):
        return (self.get_grid_volume() - self.get_volume())
        
    def fitness (self):
        #volumeTerm = self.get_unoccupied_volume() / self.get_grid_volume()
        #distanceTerm = 0
        #selfDistanceTerm = 0
        #count = 0
        #for j in self.newRegions:
        #    for i in self.oldRegions:
        #        distanceTerm += self.oldRegions[i].distance_to_region(self.newRegions[j])
        #        count += 1
        #    for i in self.newRegions:
        #        if (i != j):
        #            distanceTerm += self.newRegions[i].distance_to_region(self.newRegions[j])
        #            count +=1 
        #distanceTerm /= count
        #distanceTerm /= self.basicGrid.get_size_len()
        #print ("%.4f %.4f %.10f" % (volumeTerm, distanceTerm, volumeTerm + distanceTerm))
        #return (volumeTerm + distanceTerm)
        wallTerm = 0
        for i in self.newRegions:
            wallTerm += 0.25 / (self.basicGrid[i][0] + self.influenceRadius)
            wallTerm += 0.25 / (self.basicGrid[i][1] + self.influenceRadius)
            wallTerm += 0.25 / (self.basicGrid.get_size_len() - self.basicGrid[i][0] + self.influenceRadius)
            wallTerm += 0.25 / (self.basicGrid.get_size_len() - self.basicGrid[i][1] + self.influenceRadius)
        wallTerm /= len(self.newRegions)
        regionsTerm = 0
        count = 0
        for i in self.newRegions:
            for j in self.oldRegions:
                regionsTerm += 1.0 / (self.basicGrid[i].distance_to(self.basicGrid[j]))
                count += 1
            for j in self.newRegions:
                if (j != i):
                    regionsTerm += 1.0 / (self.basicGrid[i].distance_to(self.basicGrid[j]))
                    count += 1
        regionsTerm /= count
        return wallTerm + regionsTerm - 0.1*len(self.newRegions)
        

    def print_to_file (self, interface, fn):
        printer = VolumePrinter (interface)
        printer.drawRegion (self.basicGrid, 1.0, 'gray', 'white', 100, 100)
        for oldRegion in self.oldRegions:
            printer.drawRegion(self.oldRegions[oldRegion], 1.0, 'black', 'blue', 40, 100)
        for newRegion in self.newRegions:
            printer.drawRegion(self.newRegions[newRegion], 1.0, 'black', 'red', 40, 100)
        printer.savefig(fn)
     

if __name__ == '__main__':

    tikz = TikzInterface(1.0, 1.0)
    base = coverageIndividualBase (33, 1, 2, [0,100,544,918], 3.3)
    for i in range(100):
        individual = coverageIndividual (base, choice([1,2,3,4,5,6,7,8,9,10]))
        individual.randomize()
        individual.print_to_file(tikz, "aha_%f.pdf" % individual.fitness())
        
