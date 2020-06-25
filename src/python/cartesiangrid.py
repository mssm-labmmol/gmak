import numpy as np

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

    def getDestMask(self):
        return self.destMask
