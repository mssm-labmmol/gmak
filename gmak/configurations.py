from abc import ABC, abstractmethod
from shutil import copyfile
from copy import deepcopy
import os
import gmak.runcmd as runcmd
from gmak.config import ConfigVariables
import gmak.custom_attributes
import warnings

class ConfigurationError(RuntimeError):
    pass

class BoxFactory:
    @classmethod
    def from_file(cls, path):
        if path.endswith('.gro'):
            return cls.from_gro(path)
        elif path.endswith('.pdb'):
            return cls.from_pdb(path)
        else:
             warnings.warn(f"Unknown file extension for {path};"
                             " will use a dummy box.")
             return Box()

    @classmethod
    def from_gro(cls, path):
        with open(path, "r") as fp:
            line = ""
            for line in fp:
                continue
            box = line.split()
            edges = tuple([float(b) for b in box])
            if len(edges) == 3:
                edges = tuple([float(b) for b in box])
                return RectangularBox(*edges)
            else:
                warnings.warn(f"File {path} does not contain a "
                            "rectangular box; a dummy box will be used",
                            RuntimeWarning)
                return Box()

    @classmethod
    def from_pdb(cls, path):
        found = False
        with open(path, 'r') as fp:
            for line in fp:
                splitted = line.split()
                if splitted[0] == 'CRYST1':
                    a, b, c, alpha, beta, gamma = map(float, splitted[1:7])
                    found = True
                    break
        if not found:
            warnings.warn(f"File {path} does not contain a "
                        "box; a dummy box will be used",
                        RuntimeWarning)
            return Box()
        else:
            if (alpha != 90.00) or (beta != 90.00) or (gamma != 90.00):
                warnings.warn(f"File {path} does not contain a "
                            "rectangular box; a dummy box will be used",
                            RuntimeWarning)
                return Box()
            else:
                return RectangularBox(0.1 * a, 0.1 * b, 0.1 * c)
                
        
    @classmethod
    def from_string_list(cls, init_list):
        type = init_list[0]
        if type == 'rectangular':
            edges = (float(init_list[1]), float(init_list[2]), float(init_list[3]))
            return RectangularBox(*edges)
        elif type == 'cubic':
            edges = (float(init_list[1]), float(init_list[1]), float(init_list[1]))
            return RectangularBox(*edges)
        else: 
            raise ValueError(f"Can't construct a non-rectangular box: {init_list}.")

class Box:
    """
    So far, no general box behavior, but this may be needed in the
    future.
    """
    def __init__(self):
        pass

class RectangularBox(Box):

    def __init__(self, lx, ly, lz):
        self.lx = lx
        self.ly = ly
        self.lz = lz

    def get_edge_lengths(self):
        return (self.lx, self.ly, self.lz)

    def extend_axis(self, factor, which='z'):
        edges = (self.lx, self.ly, self.lz)
        if which == 'x':
            extended_edges = (factor * edges[0], edges[1], edges[2])
        elif which == 'y':
            extended_edges = (edges[0], factor * edges[1], edges[2])
        elif which == 'z':
            extended_edges = (edges[0], edges[1], factor * edges[2])
        else:
            raise ValueError(f"Axis {which} was not recognized.")
        return RectangularBox(*extended_edges)


class BaseConfiguration(ABC):
    """
    Base configuration class.
    """

    def get_path(self):
        """
        Returns the path of the configuration.
        """
        if os.path.isfile(self.path):
            return self.path
        else:
            raise FileNotFoundError(f"File {self.path} does not exist or is not a file.")

    def get_box(self):
        if self.box is not None:
            return self.box
        else:
            raise ValueError("Configuration has no box.")


    def copy(self, new_path, tmp=False):
        new_conf = deepcopy(self)
        copyfile(self.get_path(), new_path)
        new_conf.path = new_path
        new_conf.tmp  = tmp
        return new_conf


    def __del__(self):
        if self.tmp:
            #print(f"Removing temporary configuration {self.path}.")
            if os.path.isfile(self.path):
                os.remove(self.path)


class Configuration(BaseConfiguration,
                    gmak.custom_attributes.CustomizableAttributesMixin):

    def __init__(self, path, box=None, tmp=False):
        self.path = os.path.abspath(path)
        self.box = box
        self.tmp = tmp

    def extend_axis(self, out_path, factor, which='z', tmp=False):
        # fix periodicity
        fix_path = os.path.join(os.path.dirname(out_path), 'fix-periodicity.gro')
        runcmd.run(
            "echo 0 | %s trjconv -f %s -s %s -o %s -pbc whole" % (
                ConfigVariables.gmx,
                self.path,
                self.path,
                fix_path))
        fix_conf = ConfigurationFactory.from_file(fix_path, tmp=True)
        # then extend
        extended_box = self.get_box().extend_axis(factor, which)
        nx, ny, nz = extended_box.get_edge_lengths()
        runcmd.run(
            "%s editconf -f %s -box %f %f %f -o %s" % (ConfigVariables.gmx,
                                                       fix_conf.get_path(),
                                                       nx,
                                                       ny,
                                                       nz,
                                                       out_path))
        return Configuration(out_path, extended_box, tmp)

class ConfigurationFactory(ABC,
                           gmak.custom_attributes.CustomizableAttributesMixin):

    @classmethod
    def from_file(cls, input_file, tmp=False):
        if not os.path.isfile(input_file):
            raise FileNotFoundError(f"File {input_file} does not exist or is not a file.")
        box = BoxFactory.from_file(input_file)
        return Configuration(input_file, box, tmp)
    
    @abstractmethod
    def _construct_configuration(self, path, tmp=False):
        pass

    def construct_configuration(self, path, tmp=False):
        conf = self._construct_configuration(path, tmp)
        self.clone_custom_attributes(conf)
        return conf


class FullCoordsConfigurationFactory(ConfigurationFactory):
    """
    ConfigurationFactory for fully-constructed configurations.
    """
    def __init__(self, configuration):
        self.configuration = configuration

    def _construct_configuration(self, path, tmp=False):
        return self.configuration.copy(path, tmp)


class FollowProtocolConfigurationFactory(ConfigurationFactory):
    """
    ConfigurationFactory for recycling the last frame of a protocol.
    
    *NOT* parallelizable.
    """
    def __init__(self, protocol, grid):
        self.protocol = protocol
        self.grid = grid
        self._current_gridpoint = -1
        self._samples = None

    def _construct_configuration(self, path, tmp=False):
        if (self._current_gridpoint == -1):
            # then update samples
            self._samples = self.grid.get_samples_id()
            self._current_gridpoint = 0
        original_conf = self.protocol.get_last_frame(
            self.grid[self._samples[self._current_gridpoint]])
        new_conf = original_conf.copy(path, tmp)
        self._current_gridpoint += 1
        if self._current_gridpoint >= len(self._samples):
            self._current_gridpoint = -1
        return new_conf


class GmxLiquidConfigurationFactory(ConfigurationFactory):

    def __init__(self, conf_molec, n_molec, box):
        self.conf = conf_molec
        self.n = n_molec
        self.box = box

    def _construct_configuration(self, path, tmp=False):
        apath = os.path.abspath(path)
        edges = self.box.get_edge_lengths()
        runcmd.gmx_insert_molecules([
            '-ci', self.conf.get_path(),
            '-nmol', str(self.n),
            '-box', str(edges[0]), str(edges[1]), str(edges[2]),
            '-o', apath,
        ])
        return Configuration(apath, self.box, tmp)


GmxGasConfigurationFactory = FullCoordsConfigurationFactory


class GmxSlabConfigurationFactory(ConfigurationFactory):

    def __init__(self, conf_molec, n_molec, inner_box, axis=None, factor=None):
        self.conf = conf_molec
        self.n = n_molec
        if axis is not None:
            self.axis = axis
        else:
            self.axis = 'z'
        self.inner_box = inner_box
        if factor is not None:
            self.factor = factor
        else:
            self.factor = 5.0

    def _construct_configuration(self, path, tmp=False):
        # First create a Liquid, then expand.
        apath = os.path.abspath(path)
        liquid_path = os.path.join(os.path.dirname(apath), 'pre-liquid.gro')
        liquid_fac = GmxLiquidConfigurationFactory(self.conf, self.n, self.inner_box)
        liquid_conf = liquid_fac.construct_configuration(liquid_path, True)
        # Extend
        new_conf = liquid_conf.extend_axis(apath, self.factor, self.axis, tmp)
        return new_conf
        

class GmxSlabFollowExtendConfigurationFactory(ConfigurationFactory):
    """
    Combine a FollowProtocolConfigurationFactory with extending the
    axis of the output configuration.
    
    Differ from FollowProtocolConfigurationFactory on type string by
    using 'extend BASE_PROTOCOL AXIS FACTOR' instead of 'follow
    BASE_PROTOCOL'.
    """
    def __init__(self, protocol, grid, axis=None, factor=None):
        self.protocol = protocol
        if axis is not None:
            self.axis = axis
        else:
            self.axis = 'z'
        if factor is not None:
            self.factor = factor
        else:
            self.factor = 5.0
        self.grid = grid
        self.follow_fac = FollowProtocolConfigurationFactory(self.protocol, self.grid)

    def _construct_configuration(self, path, tmp=False):
        apath = os.path.abspath(path)
        follow_path = os.path.join(os.path.dirname(apath), 'pre-liquid.gro')
        follow_conf = self.follow_fac.construct_configuration(follow_path, True)
        # extend
        return follow_conf.extend_axis(apath, self.factor, which=self.axis, tmp=tmp)
        

class GmxSolvationConfigurationFactory(ConfigurationFactory):

    def __init__(self, conf_solute, n_solute, conf_solv, n_solv, box):
        self.conf_solute = conf_solute
        self.n_solute = n_solute
        self.conf_solv = conf_solv
        self.n_solv = n_solv
        self.box = box

    def _construct_configuration(self, path, tmp=False):
        apath = os.path.abspath(path)
        path_solv = os.path.join(os.path.dirname(apath), 'pre-solv.gro')
        # put solvent in box with temp conf
        fac_solv = GmxLiquidConfigurationFactory(self.conf_solv, self.n_solv, self.box)
        conf_solv = fac_solv.construct_configuration(path_solv, tmp=True)
        # put solute for final conf
        runcmd.gmx_insert_molecules([
            '-f', conf_solv.get_path(),
            '-ci', self.conf_solute.get_path(),
            '-nmol', str(self.n_solute),
            '-o', apath
        ])
        return Configuration(apath, self.box, tmp)

