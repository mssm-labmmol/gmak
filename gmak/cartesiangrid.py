import numpy as np
import itertools
import gmak.runcmd as runcmd

class CartesianGrid:

    def __init__(self, lens):
        self.lens = tuple(lens)
        self.dim  = len(lens)

    def getDim(self):
        return self.dim

    def getLens(self):
        return self.lens

    def getVolume(self):
        return np.prod(self.lens)

    def getCenterAsLinear(self):
        tot = self.getVolume()
        if (tot % 2) == 0:
            return int(tot/2)
        else:
            return int((tot - 1)/2)

    def getCenterAsTuple(self):
        return self.linear2tuple(self.getCenterAsLinear())

    def getDisplacement(self, i, j):
        tuple_i = self.linear2tuple(i)
        tuple_j = self.linear2tuple(j)
        displ   = tuple([y - x for x,y in zip(tuple_i, tuple_j)])
        return displ

    def linear2tuple(self, linpos):
        """Returns tuple position for linear position."""
        count = 0
        for n in np.ndindex(self.lens):
            if (count == linpos):
                return n
            count += 1
        raise ValueError("Linear position {} not found.".format(linpos))

    def tuple2linear(self, pos):
        """Returns linear position for tuple position."""
        linpos = 0
        for n in np.ndindex(self.lens):
            if (n == pos):
                return linpos
            linpos += 1
        raise ValueError("Position {} not found.".format(pos))

    def isTupleInGrid(self, tup):
        for i,val in enumerate(tup):
            if (val < 0) or (val >= self.lens[i]):
                return False
        return True

    def getCornersAsLinear(self):
        iter_arg = [(0, n-1) for n in self.lens]
        iter_obj = itertools.product(*iter_arg)
        points   = [p for p in iter_obj]
        expanded_points = [self.tuple2linear(p) for p in points]
        return expanded_points

    def getTupleIndexes(self):
        """
        Returns a list of the tuple indexes associated with the grid.
        """
        return list(np.ndindex(self.lens))

def flat2tuple(gridshape, idx):
    return CartesianGrid(gridshape).linear2tuple(idx)

class CartesianGridIterator:

    def __init__(self, cartesianGrid):
        self.src = cartesianGrid
        self.curr = -1

    def __iter__(self):
        return self

    def __next__(self):
        self.curr += 1
        if self.curr < self.src.getVolume():
            return self.src.linear2tuple(self.curr)
        raise StopIteration

class CartesianGridLinearIterator:

    def __init__(self, cartesianGrid):
        self.src = cartesianGrid
        self.curr = -1

    def __iter__(self):
        return self

    def __next__(self):
        self.curr += 1
        if self.curr < self.src.getVolume():
            return self.curr
        raise StopIteration
    
class CartesianGridMask:

    def __init__(self, srcGrid, destGrid, shiftTuple):
        if srcGrid.getDim() != destGrid.getDim():
            raise ValueError("Can't mask grids of different dimensions.")
        self.srcMask   = np.zeros(tuple(srcGrid.getLens()), dtype=int)
        self.destMask  = np.zeros(tuple(destGrid.getLens()), dtype=int)

        # build pre-destMask
        for i,pos in enumerate(CartesianGridIterator(destGrid)):
            self.destMask[pos] = -1

        # build srcMask
        for pos in CartesianGridIterator(srcGrid):
            posDestCoords = tuple([pos[i] - shiftTuple[i] for i in range(srcGrid.getDim())])
            if (destGrid.isTupleInGrid(posDestCoords)):
                fillWith = destGrid.tuple2linear(posDestCoords)
                self.destMask[posDestCoords] = srcGrid.tuple2linear(tuple(pos))
            else:
                fillWith = -1
            self.srcMask[pos] = fillWith

    def getSourceMask(self):
        return self.srcMask

    def getSourceMaskLinear(self):
        return self.srcMask.flatten()

    def getDestMask(self):
        return self.destMask

    def getDestMaskLinear(self):
        return self.destMask.flatten()

