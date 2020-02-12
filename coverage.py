from math import sqrt 
from os   import system

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
        self.elements = other.elements
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

# CubicgridRegion: a special type of VolumeRegion
class CubicgridRegion (VolumeRegion,object):

    def __init__ (self, nBoxesAtSide, boxLength, dimension):
        if (dimension == 1):
            super(CubicgridRegion,self).__init__(1)
            for i in range(nBoxesAtSide):
                self.add_element(VolumeElement(CartesianPosition([i*boxLength]), boxLength))
        else:
            self.__init__(nBoxesAtSide, boxLength, 1)
            product = self.cartesian_product(CubicgridRegion(nBoxesAtSide, boxLength, dimension - 1))
            # now translate to self, because cartesian_product returns a new element
            self.copy_from(product)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# VolumePrinter (mainly for debugging)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class VolumePrinter (object):

    def __init__ (self):

        self.fn = "/tmp/volumeprintertmpfile.tex"
        self.stream = open(self.fn, "w")



    def savefig (self, output):

        self.stream.close()



if __name__ == '__main__':
    # test
    mygrid = CubicgridRegion(33, 1.0, 2)
    radiusaround = mygrid.clone_radius_around_element(3.3, 0)
    radiusaround2 = mygrid.clone_radius_around_element(3.3, 544)
    print(mygrid)
    print(radiusaround)
    print(mygrid.get_volume())
    print(radiusaround.get_volume())
    print(radiusaround.distance_to_region(radiusaround2))

