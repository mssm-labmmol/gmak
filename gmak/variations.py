import numpy as np
import gmak.runcmd as runcmd
from abc import ABC, abstractmethod, abstractproperty
from copy import deepcopy
from gmak.cartesiangrid import *


class AbstractVariationFunction(ABC):
    """
    Contract attributes: domain_dim, image_dim
    Contract methods: apply, set_new_center, rescale
    """
    @abstractmethod
    def apply(self, list_of_args): pass

    @abstractmethod
    def set_new_center(self, i): pass

    @abstractmethod
    def set_new_origin(self, i): pass

    @abstractmethod
    def rescale(self, factors): pass

    @abstractmethod
    def get_sizes(self): pass

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
    def set_new_origin(self, i):
        raise NotImplementedError("Can't shift VariationFunctionReadFromFile.")
    def rescale(self, factors):
        raise NotImplementedError("Can't rescale VariationFunctionReadFromFile.")
    def get_sizes(self):
        return None

class VariationFunctionConstant(AbstractVariationFunction):
    def __init__(self, constants):
        self.domain_dim = len(constants)
        self.image_dim  = len(constants)
        self.constants  = constants
    def apply(self, args):
        return self.constants
    def set_new_center(self, i):
        return
    def set_new_origin(self, i):
        return
    def rescale(self, factors):
        return
    def get_sizes(self):
        return None

class VariationFunctionScale(AbstractVariationFunction):
    def __init__(self, factors):
        self.domain_dim = len(factors)
        self.image_dim  = len(factors)
        self.factors    = np.array(factors)
    def apply(self, args):
        return np.array(args) * self.factors
    def set_new_center(self, i):
        raise NotImplementedError("Can't shift VariationFunctionScale.")
    def set_new_origin(self, i):
        raise NotImplementedError("Can't shift VariationFunctionScale.")
    def rescale(self, factors):
        return
    def get_sizes(self):
        return None

class VariationFunctionFromString(AbstractVariationFunction):
    def __init__(self, domain_dim, func_strings):
        self.domain_dim = domain_dim
        self.image_dim = len(func_strings)
        self.func_strings = func_strings
        super().__init__()

    @staticmethod
    def _1d_corefunc(x, string):
        import math
        allowed_chars = "xeE0123456789+-*()./[]"
        # Evaluating exp(...) and log(...) is also allowed.
        rep_string = string.replace("exp(", "").replace("log(", "")
        for char in rep_string:
            if char not in allowed_chars:
                raise Exception("Unsafe math expression")
        return eval(string,
                    {'__builtins__': None},
                    {'log': math.log,
                     'exp': math.exp,
                     'log10': math.log10,
                     'x': x})

    def _corefunc(self, x):
        return np.array([self._1d_corefunc(x, s) for s in self.func_strings])

    def apply(self, args):
        return self._corefunc(args)

    def set_new_center(self, i):
        raise NotImplementedError

    def set_new_origin(self, i):
        raise NotImplementedError

    def rescale(self, factors):
        raise NotImplementedError

    def get_sizes(self):
        raise NotImplementedError

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
    def set_new_origin(self, i):
        displacement = self._cartesianGrid.getDisplacement(0, i)
        newStarts          = np.array(self.starts) + self.steps * np.array(displacement)
        # update data
        self.starts        = newStarts
        self._core_calcs()
    def rescale(self, factors):
        if (len(factors) != self.domain_dim):
            raise ValueError("({} != {}) len(factors) != self.domain_dim".format(len(factors), self.domain_dim))
        for f in factors:
            if not isinstance(f, int):
                raise ValueError("Subgrid factors must be integers.")
        # update defining data
        newLens  = [f*(l-1) + 1 for f,l in zip(factors, self.lens)]
        newSteps = [s/float(f) for s,f in zip(self.steps, factors)]
        self.__init__(self.starts, newSteps, newLens)

    def get_sizes(self):
        return self.lens

class AbstractVariation(ABC):
    """
    Interface:
    ----------

    int dim()
    int size()
    2D (size x dim) np.ndarray gen_data()
    void set_new_center(int)
    void set_new_origin(int)
    void rescale(list<int>)
    list<int> get_sizes()

    """
    def dim(self): return self.dim

    def size(self): return self.size

    @abstractmethod
    def gen_data(self): pass

    @abstractmethod
    def set_new_center(self, i): pass

    @abstractmethod
    def set_new_origin(self, i): pass

    @abstractmethod
    def rescale(self, factors): pass

    @abstractmethod
    def get_sizes(self): pass

    def get_type(self): return self.type_string

    def write_block_to_stream(self, stream): raise NotImplementedError

class VariationFromFunction(AbstractVariation):

    def __init__(self, dim, size, func):
        self.dim  = dim
        self.size = size
        self.func = func
        self.type_string = "function"

    def gen_data(self):
        _data = np.zeros((self.size, self.dim))
        for i in range(self.size):
            _data[i,:] = self.func.apply(i)
        return _data

    def set_new_center(self, i):
        self.func.set_new_center(i)

    def set_new_origin(self, i):
        self.func.set_new_origin(i)

    def rescale(self, factors):
        self.func.rescale(factors)
        self.size = np.prod(self.func.get_sizes())

    def get_sizes(self):
        return self.func.get_sizes()

class VariationFromFile(VariationFromFunction):
    def __init__(self, dim, size, fn):
        self.dim = dim
        self.size = size
        self.fn = fn
        self.func = VariationFunctionReadFromFile(dim, size, fn)
        self.type_string = "file"

class VariationConstant(VariationFromFunction):
    def __init__(self, dim, size, constants):
        self.dim = dim
        self.size = size
        self.func = VariationFunctionConstant(constants)
        self.type_string = "constant"

class VariationCartesian(VariationFromFunction):
    def __init__(self, dim, size, starts, steps, lens):
        self.dim = dim
        self.size = size
        self.func = VariationFunctionCartesian(starts, steps, lens)
        self.type_string = "cartesian"

    def write_block_to_stream(self, stream):
        stream.write("type cartesian\n")
        stream.write("start " + " ".join(map(str, self.func.starts)) + "\n")
        stream.write("step "  + " ".join(map(str, self.func.steps)) + "\n")
        stream.write("size "  + " ".join(map(str, self.func.lens)) + "\n")

class VariationExplicit(AbstractVariation):
    """
    Parameter values are specified explicitly in batches of size `dim'.
    """
    def __init__(self, dim, values):
        """
        dim(int): number of dimensions
        values(list of float): list of size `dim' * `size', where `size'
                               is the number of parameters
        """
        self.size = int(len(values)/dim)
        self.dim  = dim
        self.type_string = "explicit"
        self._data = np.reshape(np.array(values), (self.size, self.dim))

    def gen_data(self):
        return self._data

    def set_new_center(self, i):
        raise ValueError("Can't set new center for explicit variation.")

    def set_new_origin(self, i):
        raise ValueError("Can't set new origin for explicit variation.")

    def rescale(self, factors):
        raise ValueError("Can't rescale explicit variation.")

    def get_sizes(self):
        return (self.size,)

class VariationFromVariation(AbstractVariation):
    """
    Function is applied to output of source variation.
    """
    def __init__(self, dim, size, variation, func):
        self.dim = dim
        self.size = size
        self.variation = deepcopy(variation)
        self.func = func
        self.type_string = "from_variation"
        if (self.func.domain_dim != self.variation.dim) or (self.func.image_dim != self.dim):
            raise ValueError("Incompatible sizes in VariationFromVariation.")

    def gen_data(self):
        _inner_data = self.variation.gen_data()
        _data = np.zeros((self.size, self.dim))
        for i in range(self.size):
            _data[i, :] = self.func.apply(_inner_data[i,:])
        return _data

    def rescale(self, factors):
        self.variation.rescale(factors)
        self.size = self.variation.size

    def set_new_center(self, i):
        self.variation.set_new_center(i)

    def set_new_origin(self, i):
        self.variation.set_new_origin(i)

    def get_sizes(self):
        return self.variation.get_sizes()

class FunctionDecoratedVariation(VariationFromVariation):
    def __init__(self, variation, func_strings):
        size = variation.size
        domain_dim = variation.dim
        func_call = VariationFunctionFromString(domain_dim, func_strings)
        super().__init__(func_call.image_dim, size, variation, func_call)
        self.func_strings = func_strings

    def write_block_to_stream(self, stream):
        self.variation.write_block_to_stream(stream)
        full_function_string = " ".join(self.func_string)
        stream.write(f"function {full_function_string}\n")


class VariationFromVariationScale(VariationFromVariation):
    def __init__(self, dim, size, variation, factors):
        func = VariationFunctionScale(factors)
        super().__init__(dim, size, variation, func)
        self.type_string = "scale"

class AbstractVariationFactory:

    @staticmethod
    def parseTillEnd(stream):
        options_dict = {}
        for line in stream:
            if line[0] == '#':
                continue
            if line.rstrip() == '$end':
                break
            splittedLine             = line.split()
            identifier               = splittedLine[0]
            options                  = splittedLine[1:]
            options_dict[identifier] = options
        return options_dict

    @staticmethod
    def _createCartesian(stream, used):
        options_dict = AbstractVariationFactory.parseTillEnd(stream)
        if 'function' not in options_dict.keys():
            options_dict['function'] = None
        dim          = len(options_dict['size'])
        starts       = np.array(options_dict['start'], dtype=float)
        steps        = np.array(options_dict['step'], dtype=float)
        lens         = np.array(options_dict['size'], dtype=int)
        size         = np.prod(lens)
        return VariationCartesian(dim, size, starts, steps, lens), options_dict['function']

    @staticmethod
    def _createScale(stream, used):
        options_dict = AbstractVariationFactory.parseTillEnd(stream)
        if 'function' not in options_dict.keys():
            options_dict['function'] = None
        variation    = used.generators[0]
        factors      = np.array(options_dict['factors'], dtype=float)
        dim          = len(options_dict['factors'])
        size         = used.get_linear_size()
        return VariationFromVariationScale(dim, size, variation, factors), options_dict['function']

    @staticmethod
    def _createExplicit(stream, used):
        options_dict = AbstractVariationFactory.parseTillEnd(stream)
        if 'function' not in options_dict.keys():
            options_dict['function'] = None
        dim = int(options_dict['dim'][0])
        # values are read in batches of size dim
        values = [float(x) for x in options_dict['values']]
        return VariationExplicit(dim, values), options_dict['function']

    @staticmethod
    def _createTie(stream, used):
        options_dict = AbstractVariationFactory.parseTillEnd(stream)
        if 'function' not in options_dict.keys():
            # for this case it is necessary
            raise ValueError("A 'coupled' variation must specify a function.")
        try:
            variation    = used.generators[0]
        except AttributeError:
            raise Exception(f"In a 'coupled' variation, make sure the "
                            f"line \"using <base_variation>\" appears before the "
                            f"line \"type  coupled\" in the input file.")
        size         = used.get_linear_size()
        dim          = len(options_dict['function'])
        # no need to return func_strings below, because I am already creating a
        # FunctionDecoratedVariation
        return (FunctionDecoratedVariation(variation, options_dict['function']),
                None)

    @staticmethod
    def readFromTypeAndStream(typestring, stream, used):
        # all of them can be decorated
        func_dict = {
            'cartesian' : AbstractVariationFactory._createCartesian,
            'coupled'       : AbstractVariationFactory._createTie,
            'explicit'  : AbstractVariationFactory._createExplicit,
        }
        try:
            base, func_strings = func_dict[typestring](stream, used)
            if func_strings is not None:
                return FunctionDecoratedVariation(base, func_strings)
            else:
                return base
        except KeyError:
            raise NotImplementedError("VariationFactory for {} is not implemented.".format(typestring))

# ----------------------------------------------------------------------
# DomainSpace
# ----------------------------------------------------------------------

class DomainSpace:
    """
    Attributes:
    -----------

        2D np.ndarray data of shape (linear_size, dim)
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
        return self.data.shape[0]

    def get_dim(self):
        return self.data.shape[1]

    def get_sizes(self):
        """This method returns the sizes along each dimension when this makes
        sense.  For a cartesian 33x33 grid, for instance, it will
        return [33, 33].  Note that the lenght of the return list is
        not necessarily equal to the dimension.  For instance, we can
        have a N-lines file representing an arbitrary change in two
        parameters where get_sizes() = [N]."""
        for gen in self.generators:
            thisSizes = gen.get_sizes()
            if thisSizes is not None:
                return thisSizes
        raise ValueError("Could not find generator with appropriate get_sizes() method.")

    def set_new_center(self, i):
        for gen in self.generators:
            gen.set_new_center(i)
        self.update_data()

    def set_new_origin(self, i):
        for gen in self.generators:
            gen.set_new_origin(i)
        self.update_data()

    def rescale(self, factors):
        for gen in self.generators:
            gen.rescale(factors)
        self.update_data()

    def write_to_stream(self, stream):
        m, n = self.data.shape
        for i in range(m):
            for j in range(n):
                stream.write("%18.7e" % self.data[i,j])
            stream.write('\n')

    def write_to_file(self, fn):
        np.savetxt(fn, self.data)

    def write_block_to_stream(self, stream):
        for gen in self.generators:
            gen.write_block_to_stream(stream)

    def get_data(self):
        return self.data

class DomainSpaceFactory:

    @staticmethod
    def readFromTypeAndStream(typestring, stream, used):
        variation = AbstractVariationFactory.readFromTypeAndStream(typestring, stream, used)
        return DomainSpace([variation])
