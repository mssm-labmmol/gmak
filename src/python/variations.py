import numpy as np
from abc import ABC, abstractmethod, abstractproperty
from copy import deepcopy
from cartesiangrid import *

class AbstractVariationFunction(ABC):
    """
    Contract attributes: domain_dim, image_dim
    Contract methods: apply, set_new_center
    """
    @abstractmethod
    def apply(self, list_of_args): pass

    @abstractmethod
    def set_new_center(self, i): pass

class VariationFunctionReadFromFile(AbstractVariationFunction):
    def __init__(self, dim, size, fn):
        self.fn  = fn
        self.domain_dim = dim
        self.image_dim = dim
        self.size = size
        self._data = np.loadtxt(self.fn)
        self._data = self._data.reshape(self._data.shape[0], -1)
        if (dim != self._data.shape[1]):
            raise ValueError("Wrong number of dimensions: expected {} in file {}, but got {}.".format(dim, self.fn, self._data.shape[1]))
        if (size != self._data.shape[0]):
            raise ValueError("Wrong size: expected {} in file {}, but got {}.".format(size, self.fn, self._data.shape[0]))
    def apply(self, args):
        return self._data[args]
    def set_new_center(self, i):
        raise NotImplementedError("Can't shift VariationFunctionReadFromFile.")

class VariationFunctionConstant(AbstractVariationFunction):
    def __init__(self, constants):
        self.domain_dim = len(constants)
        self.image_dim  = len(constants)
        self.constants  = constants
    def apply(self, args):
        return self.constants
    def set_new_center(self, i):
        return

class VariationFunctionScale(AbstractVariationFunction):
    def __init__(self, factors):
        self.domain_dim = len(factors)
        self.image_dim  = len(factors)
        self.factors    = np.array(factors)
    def apply(self, args):
        return np.array(args) * self.factors
    def set_new_center(self, i):
        raise NotImplementedError("Can't shift VariationFunctionScale.")

class VariationFunctionCartesian(AbstractVariationFunction):
    def __init__(self, starts, steps, lens):
        self.starts = starts
        self.steps  = steps
        self.lens   = lens
        self.domain_dim = len(starts)
        self.image_dim  = len(starts)
        self._cartesianGrid = CartesianGrid(lens)
        self._core_calcs()
    def _core_calcs(self):
        self._raveled = np.zeros((self.domain_dim, self._cartesianGrid.getVolume()))
        for i, pt in enumerate(CartesianGridIterator(self._cartesianGrid)):
            arr = np.array(pt)
            self._raveled[:,i] = np.array(self.starts) + self.steps * arr
    def apply(self, args):
        return self._raveled[:,args]
    def set_new_center(self, i):
        currentCenter      = self._cartesianGrid.getCenterAsLinear()
        centerDisplacement = self._cartesianGrid.getDisplacement(currentCenter, i)
        newStarts          = np.array(self.starts) + self.steps * np.array(centerDisplacement)
        # update data
        self.starts        = newStarts
        self._core_calcs()
        
class AbstractVariation(ABC):
    """
    Interface:
    ----------
    
    int dim()
    int size()
    2D (size x dim) np.ndarray gen_data()
    void set_new_center(int)

    """
    def dim(self): return self.dim

    def size(self): return self.size

    @abstractmethod
    def gen_data(self): pass

    @abstractmethod
    def set_new_center(self, i): pass

class VariationFromFunction(AbstractVariation):
    
    def __init__(self, dim, size, func):
        self.dim  = dim
        self.size = size
        self.func = func

    def gen_data(self):
        _data = np.zeros((self.size, self.dim))
        for i in range(self.size):
            _data[i,:] = self.func.apply(i)
        return _data

    def set_new_center(self, i):
        self.func.set_new_center(i)

class VariationFromFile(VariationFromFunction):
    def __init__(self, dim, size, fn):
        self.dim = dim
        self.size = size
        self.fn = fn
        self.func = VariationFunctionReadFromFile(dim, size, fn)

class VariationConstant(VariationFromFunction):
    def __init__(self, dim, size, constants):
        self.dim = dim
        self.size = size
        self.func = VariationFunctionConstant(constants)

class VariationCartesian(VariationFromFunction):
    def __init__(self, dim, size, starts, steps, lens):
        self.dim = dim
        self.size = size
        self.func = VariationFunctionCartesian(starts, steps, lens)

class VariationFromVariation(AbstractVariation):
    """
    Function is applied to output of source variation.
    """
    def __init__(self, dim, size, variation, func):
        self.dim = dim
        self.size = size
        self.variation = deepcopy(variation)
        self.func = func
        if (self.func.domain_dim != self.variation.dim) or (self.func.image_dim != self.dim) or (self.dim != variation.dim):
            raise ValueError("Incompatible sizes in VariationFromVariation.")

    def gen_data(self):
        _data = np.zeros((self.size, self.dim))
        for i in range(self.size):
            _data[i, :] = self.func.apply( self.variation.func.apply(i) )
        return _data

    def set_new_center(self, i):
        self.variation.set_new_center(i)

class VariationFromVariationScale(VariationFromVariation):
    def __init__(self, dim, size, variation, factors):
        func = VariationFunctionScale(factors)
        super().__init__(dim, size, variation, func)

class DomainSpace:
    """
    Attributes:
    -----------

        2D np.ndarray data
        AbstractVariation generators[]

    Methods:
    """
    def __init__(self, generators):
        self.generators = generators
        self.update_data()

    def update_data(self):
        sparse_data     = []
        for gen in self.generators:
            sparse_data.append( gen.gen_data() ) 
        self.data = np.concatenate(sparse_data, axis=1)
        self.data = self.data.reshape(self.data.shape[0], -1)

    def get(self, i):
        return self.data[i,:]

    def get_linear_size(self):
        return data.shape[0]

    def get_dim(self):
        return data.shape[1]

    def set_new_center(self, i):
        for gen in self.generators:
            gen.set_new_center(i) 
        self.update_data()

    def write_to_file(self, fn):
        np.savetxt(fn, self.data)

