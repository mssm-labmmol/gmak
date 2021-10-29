from abc import ABC, abstractmethod
from shutil import copyfile
from copy import deepcopy
import os

class Box:
    pass

class RectangularBox(Box):

    def __init__(self, lx, ly, lz):
        self.lx = lx
        self.ly = ly
        self.lz = lz

    def get_edge_lengths(self):
        return (self.lx, self.ly, self.lz)

    @classmethod
    def from_ext(self, path):
        pass

    @classmethod
    def from_gro(self, path):
        pass

    @classmethod
    def from_pdb(self, path):
        pass


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


    def copy(self, new_path):
        new_conf = deepcopy(self)
        shuftil.copyfile(self.get_path(), new_path)
        new_conf.path = new_path
        return new_conf


class Configuration(BaseConfiguration):

    def __init__(self, path, box=None):
        self.path = path
        self.box = box

    @classmethod
    def from_ext(cls, input_file):
        pass

    @classmethod
    def from_gro(cls, input_file):
        pass

    @classmethod
    def from_pdb(cls, input_file):
        pass

class GmxConfiguration(BaseConfiguration):

    def extend_axis(self, out_path, factor, which='z'):
        pass
    

class ConfigurationFactory(ABC):
    
    @abstractmethod
    def construct_configuration(self, path):
        pass


class FullCoordsConfigurationFactory(ConfigurationFactory):
    """
    ConfigurationFactory for fully-constructed configurations.
    """
    def __init__(self, configuration):
        self.configuration = configuration

    def construct_configuration(self, path):
        self.configuration.copy()


class FollowProtocolConfigurationFactory(ConfigurationFactory):
    """
    ConfigurationFactory for recycling the last frame of a protocol.
    """
    def __init__(self, protocol, grid):
        self.protocol = protocol
        self.grid = grid
        self._current_gridpoint = 0

    def construct_configuration(self, path):
        original_conf = self.grid[self._current_gridpoint].getLastFrame(self.protocol.name)
        self._current_gridpoint += 1
        if self._current_gridpoint >= self.grid.get_linear_size():
            self._current_gridpoint = 0


class GmxLiquidConfigurationFactory(ConfigurationFactory):
    pass

class GmxGasConfigurationFactory(ConfigurationFactory):
    pass

class GmxSlabConfigurationFactory(ConfigurationFactory):
    pass

class GmxSlabFollowExtend(ConfigurationFactory):
    """
    Combine a FollowProtocolConfigurationFactory with extending the
    axis of the output configuration.
    
    Differ from FollowProtocolConfigurationFactory on type string by
    using 'extend BASE_PROTOCOL AXIS FACTOR' instead of 'follow
    BASE_PROTOCOL'.
    """
    pass

class GmxSolvationConfigurationFactory(ConfigurationFactory):
    pass
