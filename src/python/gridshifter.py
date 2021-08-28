import os
import runcmd
from logger import *
from parameters import *
from cartesiangrid import *
from gridbase import *

# # creates (soft) symbolic link from target to link_name
# def create_symbolic_link (target, link_name):
#     if not (os.path.isdir(link_name)):
#         print ("Linking: ln -s %s %s" % (target, link_name))
#         runcmd.run("ln -s %s %s" % (target, link_name))
#     else:
#         print ("Warning: Link %s already exists.\n" % link_name)

class GridShifter:

    def __init__ (self, grid, maxshifts, margins, ncut, keepsamples):

        # Margins delimiting the border region.  self.margins = {}
        self.grid      = grid
        self.maxshifts = maxshifts
        self.nshifts   = 0
        # # Default values means no shifting
        # self.margins = [[0.0, 1.0] for i in range(grid.get_dim())]
        self.margins = margins

        # Fraction of the best points which, if inside the border region,
        # trigger the shifting of the grid.
        self.ncut = ncut

        # Do we keep samples from old grids?
        self.keepSamples = keepsamples

        # Stores the CG.
        self.cg = tuple([0 for i in range(grid.get_dim())])

    def getCurrentNumberOfShifts(self):
        return self.nshifts

    def getCGasLinear(self):
        return self.grid.tuple2linear(self.cg)
            
    # Calculates CG of the best ncut percent points.
    # Returns it as a tuple.
    # Store in internal variables.
    def calcCG (self, optimizer):

        thr = int(self.ncut * self.grid.get_linear_size())
        cg  = [0 for i in range(len(self.grid.get_size()))]

        for k in range(thr):
            idx = optimizer.stateScores[k][0]
            point = self.grid.linear2tuple(idx)
            for i, x in enumerate(cg):
                cg[i] += point[i]

        for i, x in enumerate(cg):
            cg[i] = int(cg[i]/thr)

        globalLogger.putMessage('MESSAGE: CG position is {}'.format(cg))
        self.cg = tuple(cg)
        return self.cg

    # With CG calculated, verify if we need to shift.
    def doWeShift (self, optimizer):
        if (self.nshifts >= self.maxshifts):
            return False
        self.calcCG(optimizer)
        grid_size = self.grid.get_size()
        for i in range(len(grid_size)):
            dim_min = int(self.margins[i][0] * grid_size[i])
            dim_max = int(self.margins[i][1] * grid_size[i])
            cg_i    = self.cg[i]
            if (cg_i < dim_min) or (cg_i > dim_max):
                return True
        return False

    def checkShift (self, optimizer):
        if not (self.doWeShift(optimizer)):
            globalLogger.putMessage('MESSAGE: The grid *WAS NOT* shifted to CG = {} (linear = {})'.format(self.cg, self.getCGasLinear()), dated=True)
            return False
        globalLogger.putMessage('MESSAGE: The grid was shifted to CG = {} (linear = {})'.format(self.cg, self.getCGasLinear()), dated=True)
        # self.shift(grid, paramLoop, old_workdir, new_workdir)
        self.nshifts += 1
        return True

    def shift(self, optimizer):
        if not(self.checkShift(optimizer)):
            return False

        linearCG = self.getCGasLinear()
        tupleCG  = self.cg
        cartesianGrid = self.grid.getCartesianGrid()
        
        # Alter parameter space
        self.grid.setNewCenterForParameters(linearCG)
        
        # Write new topologies 
        self.grid.incrementPrefixOfTopologies()
        self.grid.writeTopologies()

        # Set center as sample.
        self.grid[linearCG].set_as_sample()

        shiftTuple = cartesianGrid.getDisplacement(cartesianGrid.getCenterAsLinear(), linearCG)
        
        # New gridpoints.
        newGridpoints = [None for i in range(self.grid.get_linear_size())]
        mask = CartesianGridMask(self.grid.getCartesianGrid(), self.grid.getCartesianGrid(), shiftTuple)
        maskArray = mask.getDestMaskLinear()

        # Set new gridpoints.
        for i, m in enumerate(maskArray):
            if (m != -1) and (self.keepSamples):
                newGridpoints[i] = self.grid[m]
                newGridpoints[i].resetWithNewId(i)
            else:
                newGridpoints[i] = GridPoint(self.grid, i)

        # Put new gridpoints in grid
        self.grid.setGridpoints(newGridpoints)
        return True

    @staticmethod
    def createFromGridAndDict (grid, dictargs):
        if dictargs is not None:
            maxshifts = int(dictargs['maxshifts'][0])
            margins   = [[float(dictargs['margins'][2*i]), float(dictargs['margins'][2*i+1])] for i in range(grid.get_dim())]
            ncut      = float(dictargs['ncut'][0])
            keepsamples = True if dictargs['keepsamples'] == 'yes' else False
        else:
            maxshifts = 0
            margins   = None
            ncut      = 0.0
            keepsamples = True
        return GridShifter(grid, maxshifts, margins, ncut, keepsamples)

class EmptyGridShifter(GridShifter):
    def __init__(self):
        return
    def shift(self):
        raise NotImplementedError("Using EmptyGridShifter for shifting.")
    
