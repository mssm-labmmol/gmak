import os
import numpy as np
import gmak.runcmd as runcmd
import gmak.logger as logger
from abc import ABC, abstractmethod
from gmak.parameter_space import *
from gmak.cartesiangrid import *
from gmak.gridbase import *
from copy import deepcopy

# # creates (soft) symbolic link from target to link_name
# def create_symbolic_link (target, link_name):
#     if not (os.path.isdir(link_name)):
#         print ("Linking: ln -s %s %s" % (target, link_name))
#         runcmd.run("ln -s %s %s" % (target, link_name))
#     else:
#         print ("Warning: Link %s already exists.\n" % link_name)


class IndexTransform:
    """
    This class transforms a tuple index that is specific to a grid-shift
    iteration into a unique and general tuple index.
    """
    def __init__(self, dim):
        self.origin = tuple([0 for i in range(dim)])

    def update(self, origin):
        if self.origin is None:
            self.origin = -np.array(origin, dtype=int)
        else:
            self.origin -= np.array(origin, dtype=int)

    def transform(self, tuple_index):
        return tuple(np.array(tuple_index, dtype=int) - self.origin)


class BaseGridShifter(ABC):
    """
    Attributes
    ----------
    
    grid : ParameterGrid
        The ParameterGrid object of the run.

    index_transform : IndexTransform
        An IndexTransform instance.

    nshifts : int
        Current shift iteration.
    
    maxshifts : int
        Maximum number of shifts.

    keepsamples : boolean
        Indicates if data for previously simulated points is kept
        after shifting.

    """

    @abstractmethod
    def calc_new_origin(self,
                        tuple_indexes,
                        scores,
                        propnames,
                        averages,
                        uncertainties):
        """
        Returns the tuple or linear index for the new origin, or None
        if shifting is not desired.
        
        Parameters
        ----------
        tuple_indexes : list of size N
            Flat list of tuple indexes with size D along linear indexes.
        
        scores : list of size N
            Flat list of scores along linear indexes.

        propnames:  list of size P
            List of property names.

        averages: np.ndarray of shape (N,P)
            Matrix with average values estimated for each combination
            of linear index and property.

        uncertainties: np.ndarray of shape (N,P)
            Matrix with uncertainties estimated for each combination
            of linear index and property.
        
        Returns
        -------
        index : int, tuple of size D, or None
             Linear or tuple index of the new origin, or None not to shift.
        """
        pass


    def get_current_number_of_shifts(self):
        return self.nshifts


    def shift(self, optimizer):
        """
        Template for shifting.
        Returns True if shifting occurred, False otherwise.
        """
        if (self.nshifts >= self.maxshifts):
            logger.globalLogger.putMessage('MESSAGE: The grid *WAS NOT* shifted because maxshifts was reached.', dated=True)
            return False

        # get cartesian grid
        cartesian_grid = self.grid.getCartesianGrid()

        # prepare data to calculate new origin
        scores = deepcopy(optimizer.stateScores)
        scores.sort(key=lambda x: x[0])
        scores = [s[1] for s in scores]
        tuple_indexes = cartesian_grid.getTupleIndexes()
        prop_names = optimizer.getProperties()
        averages = np.array([[gp.get_property_estimate(prop) for prop in prop_names] for gp in self.grid])
        uncertainties = np.array([[gp.get_property_err(prop) for prop in prop_names] for gp in self.grid])

        # new origin
        new_origin = self.calc_new_origin(tuple_indexes, scores, prop_names, averages, uncertainties)

        if new_origin is None:
            logger.globalLogger.putMessage('MESSAGE: The grid *WAS NOT* shifted.', dated=True)
            return False

        # increment shifts
        self.index_transform.update(new_origin)
        logger.globalLogger.putMessage(f"MESSAGE: The grid *WAS* shifted to index = {new_origin}.", dated=True)
        self.nshifts += 1

        # set new origin for variations (note that cartesian grid does
        # not change)
        self.grid.setNewOriginForParameters(new_origin)

        # New gridpoints.
        new_gridpoints = [None for i in range(self.grid.get_linear_size())]
        mask = CartesianGridMask(cartesian_grid,
                                 cartesian_grid,
                                 new_origin)
        mask_array = mask.getDestMaskLinear()

        # Set new gridpoints.
        for i, m in enumerate(mask_array):
            if (m != -1) and (self.keepsamples):
                new_gridpoints[i] = self.grid[m]
                new_gridpoints[i].resetWithNewId(i)
            else:
                new_gridpoints[i] = GridPoint(self.grid, i)

        # Put new gridpoints in grid
        self.grid.setGridpoints(new_gridpoints)

        return True

    def merge(self, other):
        self.maxshifts = other.maxshifts



class GridShifter(BaseGridShifter):
    """
    Default GridShifter.
    """
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
        self.keepsamples = keepsamples

        self.index_transform = IndexTransform(self.grid.get_dim())

    def merge(self, other):
        self.maxshifts = other.maxshifts
        self.margins = other.margins
        self.ncut = other.ncut

    @classmethod
    def create_from_grid_and_dict(cls, grid, dictargs):
        if dictargs is not None:
            maxshifts = int(dictargs['maxshifts'][0])
            margins   = [[float(dictargs['margins'][2*i]), float(dictargs['margins'][2*i+1])] for i in range(grid.get_dim())]
            ncut      = float(dictargs['ncut'][0])
            keepsamples = False #True if dictargs['keepsamples'][0] == 'yes' else False
        else:
            maxshifts = 0
            margins   = None
            ncut      = 0.0
            keepsamples = True
        return cls(grid, maxshifts, margins, ncut, keepsamples)

    def calc_new_origin(self,
                        tuple_indexes,
                        scores,
                        propnames,
                        averages,
                        uncertainties):
        # Calculate CG.
        sorted_indexes = np.argsort(scores)
        thr = int(self.ncut * self.grid.get_linear_size())
        cg  = [0 for i in range(len(self.grid.get_size()))]
        for k in range(thr):
            idx = sorted_indexes[k]
            point = self.grid.linear2tuple(idx)
            for i, x in enumerate(cg):
                cg[i] += point[i]
        for i, x in enumerate(cg):
            cg[i] = int(cg[i]/thr)
        # Convert CG to a tuple.
        cg = tuple(cg)
        logger.globalLogger.putMessage('MESSAGE: CG position is {}'.format(cg))
        # Check if it lies in the inner or outer part of the grid.
        grid_size = self.grid.get_size()
        for i in range(len(grid_size)):
            dim_min = int(self.margins[i][0] * grid_size[i])
            dim_max = int(self.margins[i][1] * grid_size[i])
            if (cg[i] < dim_min) or (cg[i] > dim_max):
                # It lies in the outer part -> Shift!
                cartesian_grid = self.grid.getCartesianGrid()
                current_center = cartesian_grid.getCenterAsLinear()
                cg_as_linear   = cartesian_grid.tuple2linear(cg)
                shift_tuple    = cartesian_grid.getDisplacement(current_center, cg_as_linear)
                return shift_tuple
        # It lies in the inner part -> Don't shift!
        return None


class EmptyGridShifter(BaseGridShifter):
    def __init__(self):
        return
    def calc_new_origin(self,
                        tuple_indexes,
                        scores,
                        propnames,
                        averages,
                        uncertainties):
        raise NotImplementedError("Using EmptyGridShifter for shifting.")


class CustomGridShifter(BaseGridShifter):
    def __init__(self, calculator):
        self.calculator = calculator
        self.nshifts = 0
        self.maxshifts = 0
        self.keepsamples = False

    def calc_new_origin(self,
                        tuple_indexes,
                        scores,
                        propnames,
                        averages,
                        uncertainties):
        return self.calculator(tuple_indexes,
                               scores,
                               propnames,
                               averages,
                               uncertainties)


    @classmethod
    def from_dict(cls, bd, grid):
        for key in ['type', 'maxshifts']:
            if (key not in bd.keys()):
                raise Exception(f"CustomGridShifter must specify a '{key}'.")
        bd['maxshifts'][0] = int(bd['maxshifts'][0])
        out = CustomGridShifterFactory.create(bd['type'][0])
        for k in bd.keys():
            if len(bd[k]) == 1:
                setattr(out, k, bd[k][0])
            else:
                setattr(out, k, bd[k])
        out.grid = grid
        out.index_transform = IndexTransform(out.grid.get_dim())
        return out


class CustomGridShifterFactory:

    # this is a dict 'type' : calculate_new_origin function
    ptable = {}

    @classmethod
    def add_custom_gridshifter(cls, type_name, calculator):
        cls.ptable[type_name] = calculator

    @classmethod
    def create(cls, type_name):
        return CustomGridShifter(cls.ptable[type_name])


def create_gridshifter(grid, bd):
    if ('type' not in bd.keys()) or (bd['type'][0] == 'default'):
        # if type is not specified or is 'default', use default
        return GridShifter.create_from_grid_and_dict(grid, bd)
    # otherwise, it is a custom gridshifter
    try:
        return CustomGridShifter.from_dict(bd, grid)
    except KeyError:
        raise NotImplementedError(f"Can't create GridShifter of type {bd['type'][0]}.")
