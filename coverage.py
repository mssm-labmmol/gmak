from math import sqrt 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Cartesian Position
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class CartesianPosition:

    def __init__ (self, r):
        self.dim = len(r)
        self.r   = r

    def __eq__ (self, other):
        return ((self.dim == other.dim) and (self.r == other.r))

    def __repr__ (self):
        out_str = "["
        for i in range(self.dim):
            out_str += "%10.3e" % self.r[i]
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

class VolumeElement:

    def __init__ (self, position, volume):
        self.position = position
        self.volume   = volume

    def __eq__ (self,other):
        return ((self.volume == other.volume) and (self.position == other.position))

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
        newElement = VolumeElement (dimensional_stack(self.r,other.r), self.volume*other.volume)
        return newElement

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# VolumeRegion
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class VolumeRegion:
    
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

    def clone_radius_around_element (self, other, radius, i):
        newRegion = VolumeRegion (self.dimension)
        for el in other.elements:
            if (el.is_within_distance_of(other.elements[i], radius)):
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
            for el_j in self.elements:
                newRegion.add_element(el_i.dimensional_stack(el_j))
        return newRegion


# CubicgridRegion: a special type of VolumeRegion
class CubicgridRegion (VolumeRegion):

    def __init__ (self, nBoxesAtSide, boxLength, dimension):

        if (dimension == 1):
            self = VolumeRegion (1)
            for i in range(nBoxesAtSide):
                self.add_element(VolumeElement(CartesianPosition([i*boxLength]), boxLength))
            return self
        else:
            self = CubicgridRegion(nBoxesAtSide, boxLength, 1)
            self.cartesian_product(CubicgridRegion(nBoxesAtSide, boxLength, dimension - 1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# VolumePrinter (mainly for debugging)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


