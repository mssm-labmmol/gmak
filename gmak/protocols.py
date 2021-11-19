from enum import Enum, auto
from abc import ABC, abstractmethod
import gmak.runcmd as runcmd
from gmak.simulate import *
from shutil import copyfile
from gmak.surrogate_model import *
from gmak.custom_attributes import CustomizableAttributesMixin
import warnings
import re
import copy
import os
import gmak.logger as logger
from gmak.mdputils import *
import gmak.component_properties as component_properties
from gmak.traj_ana import get_box
from gmak.configurations import ConfigurationFactory, FollowProtocolConfigurationFactory

class Ensemble(Enum):
    NVT = 0
    NPT = 1

    @classmethod
    def from_string(cls, string):
        return cls(['nvt', 'npt'].index(string))


class DefaultExtendMixin:
    """
    Assumes the following attributes:
        - self.mdps
        - self.maxSteps
        - self.minFactor
    """
    @classmethod
    def _ext_from_dict(cls, bd):
        try:
            maxSteps = int(bd['maxsteps'][0])
        except KeyError:
            maxSteps = None
        try:
            minFactor = float(bd['minextend'][0])
        except KeyError:
            minFactor = None
        return maxSteps, minFactor

    def calc_initial_len(self):
        """
        Returns the initial length of the simulations.
        """
        return int(read_nsteps_from_mdp(self.mdps[-1]))


    def calc_extend(self, gridpoint, optimizer, protocolsHash):
        """
        Returns the new number of steps in the protocol, or None if there
        is no need to extend the protocol.

        The new number of steps will always be less than self.maxSteps
        to avoid very long simulations. Also, it will always be at
        least self.minFactor * (the current number of steps), where,
        by default, self.minFactor = 1.1.
        """
        errs_tols = gridpoint.get_errs_tols(optimizer, self, protocolsHash)
        last_lenspec = gridpoint.getProtocolLength(self)
        maxSteps = self.maxSteps
        lenspec = None
        for prop in errs_tols.keys():
            err = errs_tols[prop]['err']
            tol = errs_tols[prop]['tol']
            if err > tol:
                scalingFactor = (err / tol) ** 2
                if (lenspec is None) or (scalingFactor > lenspec):
                    lenspec = scalingFactor
        if lenspec is not None:
            if lenspec < self.minFactor:
                logger.globalLogger.putMessage(f"Adjusted extension factor to minimum "
                                               f"value {self.minFactor}.")
                lenspec = self.minFactor
            lenspec = int(last_lenspec * lenspec)
            if lenspec > maxSteps:
                logger.globalLogger.putMessage(f"Tried to set length to {lenspec},"
                                               f" but capped it to {maxSteps}.")
                lenspec = maxSteps
            if lenspec == last_lenspec:
                return # None
            else:
                return lenspec
        else:
            return # None



class BaseProtocol(ABC):
    
    def get_filtering_properties(self):
        at_properties = self.get_component_properties()
        # filter those that are timeseries
        standard_properties = []
        for name, obj in at_properties.items():
            if obj.is_timeseries:
                standard_properties.append(name)
                
        if (('potential' not in standard_properties) and
            (self.get_mbar_model() is not None)):
            standard_properties.append('potential')
            
        if self.has_pv():
            if (('pV' not in standard_properties) and
                (self.get_mbar_model() is not None)):
                standard_properties.append('pV')
        return standard_properties


    def get_nonfiltering_properties(self):
        raw_properties = self.get_properties()
        filter_properties = self.get_filtering_properties()
        return [p for p in raw_properties if p not in
                filter_properties]


    def get_properties (self):
        return self.properties


    def get_component_properties(self):
        return self.component_properties


    def add_component_property(self, ap_obj):
        self.component_properties[ap_obj.name] = ap_obj


    def add_property(self, prop: str):
        try:
            if prop not in self.properties:
                self.properties.append(prop)
        except AttributeError:
            self.properties = []
            if prop not in self.properties:
                self.properties.append(prop)


    def get_temperature(self):
        return self.get_custom_attributes().temperature


    def add_surrogate_model (self, surrogate_model_string,
                             component_property, bool_legacy):
        """Add a (surrogate_model, property) tuple to the protocol. """
        model = init_surrogate_model_from_string(
            surrogate_model_string, bool_legacy)
        try:
            self.surrogate_models.append((model, component_property))
        except AttributeError:
            # this probably means the list has not been initialized
            self.surrogate_models = []
            self.surrogate_models.append((model, component_property))


    def get_mbar_model (self):
        """Returns the first instance of an MBAR model for this protocol."""
        for m, p in self.surrogate_models:
            if (m.kind == 'mbar'):
                return m


    def get_interp_models (self):
        """Returns a list of instances of Interpolation."""
        out_models = []
        for m, p in self.surrogate_models:
            if not m.requiresReweight():
                out_models.append(m)
        return out_models


    def get_interp_models_props (self):
        """Returns a list of tuples (instance of Interpolation, property)."""
        out = []
        for m, p in self.surrogate_models:
            if not m.requiresReweight():
                out.append((m,p))
        return out


    def get_all_models_props(self):
        """Returns a list of tuples (instance of Surrogate Model, property)."""
        out = []
        for m, p in self.surrogate_models:
                out.append((m,p))
        return out


    def get_avg_err_estimate_of_property (self, prop, kind):
        """
        Returns a tuple (EA_k, dEA_k) for atomic property 'prop' estimated
        via 'kind' surrogate model for all states.

        """
        if (kind == 'mbar'):
            model = self.get_mbar_model()
            index = self.get_reweighting_properties().index(prop)
        else:
            index = 0
            for m, p in self.surrogate_models:
                if (p == prop) and (m.kind == kind):
                    model = m
                    break
        return model.EA_pk[index,:], model.dEA_pk[index,:]

    def get_avg_err_estimate_of_property_for_gridpoint(self, gridpoint, prop, kind):
        EA, dEA = self.get_avg_err_estimate_of_property(prop, kind)
        return EA[gridpoint.id], dEA[gridpoint.id]
    
    def requires_corners (self):
        """Checks if this protocol requires corners."""
        for m, p in self.surrogate_models:
            if m.requiresCorners():
                return True
        return False

    
    def requires_reweight (self):
        """Checks if this protocol requires reweight."""
        for m, p in self.surrogate_models:
            if m.requiresReweight():
                return True
        return False


    def get_reweighting_properties (self):
        """Based on the surrogate models for each property, return a list
        with the atomic properties that need reweighting."""
        if self.get_mbar_model() is None:
            return []
        else:
            output_list = ['potential']
            if (self.has_pv()):
                output_list.append('pV')
            for m, p in self.surrogate_models:
                if m.requiresReweight():
                    if p not in output_list:
                        output_list.append(p)
            return output_list


    def get_nonreweighting_properties (self):
        """Based on the surrogate models for each property, return a list
        with the atomic properties that DO NOT need reweighting."""
        output_list = []
        for m, p in self.surrogate_models:
            if not (m.requiresReweight()):
                if p not in output_list:
                    output_list.append(p)
        return output_list


    def has_pv(self):
        return self._has_pv

    @abstractmethod
    def get_last_frame(self, gridpoint):
        pass


class GmxBaseProtocol(BaseProtocol):
    """
    Gromacs-specific functions.
    """
    
    # Override
    def get_temperature(self):
        temp = 0.0
        fp = open(self.mdps[-1], "r")
        for line in fp:
            if re.match(r"^ref[_-]t", line):
                temp = float(line.split()[2])
        fp.close()
        return temp

    def get_last_frame(self, gridpoint):
        return ConfigurationFactory.from_file(
            gridpoint.protocol_outputs[self.name]['gro'])

    
class GmxCoreAlchemicalProtocol(GmxBaseProtocol):

    def __init__ (self,
                  name,
                  system,
                  conf_fac,
                  mdps,
                  properties):
        self.name = name
        self.system = system
        self.mdps = mdps
        # ConfigurationFactory is FullCooords... on the previous
        # lambda for lambda != 0 and GmxSolvation... for lambda == 0.
        self.conf_fac = conf_fac
        self.properties = properties
        self.surrogate_models = []
        self.parent = None
        self.output_files = None

    def run_gridpoint_at_dir (self, gridpoint, workdir, simulate=True):
        labels = [str(x) for x in range(len(self.mdps))]
        nsteps = gridpoint.getProtocolLength(self.parent)
        conf   = self.conf_fac.construct_configuration(
            os.path.join(workdir, "solvation.gro"))
        out = simulate_protocol_solvation(
            conf.get_path(),
            gridpoint.getTopologyPath(self.system),
            self.mdps,
            nsteps,
            labels,
            workdir,
            simulate)
        self.output_files = out
        gridpoint.add_protocol_output(self, out)
        return out

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        if self.requires_reweight():
            return self.run_gridpoint_at_dir(gridpoint, workdir, simulate=False)

    @classmethod
    def from_dict(cls, bd, coordinates, grid):
        # read number of lambda points from some mdp file
        imu = mdpUtils()
        imu.parse_file(bd['mdps'][-1])
        nlambdas = imu.get_nlambdas()
        # force calc-lambda-neighbors to be -1
        if (int(imu.get_option('calc-lambda-neighbors')[0]) != -1):
            warnings.warn('calc-lambda-neighbors will be set to -1')
        outputProtocols = []
        for i in range(nlambdas):
            # write mdp files for each state, setting lambda appropriately
            mdps = []
            for m in bd['mdps']:
                prefix = os.path.abspath(m)[:-4]
                fn = prefix + '_{}.mdp'.format(i)
                mdps.append(fn)
                imu.clean()
                imu.parse_file(m)
                imu.set_lambda_state(i)
                imu.set_option('calc-lambda-neighbors', ['-1'])
                imu.write_to_file(fn)
            # first state follows nothing, but succeeding states follow their
            # respective previous state
            if (i == 0):
                conf_fac = coordinates[bd['coords'][0]]
            else:
                conf_fac = FollowProtocolConfigurationFactory(
                    outputProtocols[i-1], grid)
            # create protocol for each state and append it to the list
            # of states returned by this factory
            outputProtocols.append(
                cls(bd['name'][0] + '-{}'.format(i),
                    bd['system'][0],
                    conf_fac,
                    mdps,
                    []))
        return outputProtocols


class GmxAlchemicalProtocol(GmxBaseProtocol,
                            DefaultExtendMixin,
                            CustomizableAttributesMixin):

    def __init__(self,
                 name,
                 core_protocols,
                 mdps,
                 properties,
                 ensemble,
                 maxSteps=None,
                 minFactor=None):
        self.name = name
        self.core_protocols = core_protocols
        self.ensemble = ensemble
        # Set _has_pv just to avoid breaking existing code.
        if self.ensemble == Ensemble.NPT:
            self._has_pv = True
        else:
            self._has_pv = False
        self.mdps = mdps
        self.properties = properties
        self.surrogate_models = []
        if maxSteps is None:
            raise ValueError("GmxAlchemicalProtocol must define maxsteps.")
        else:
            self.maxSteps = maxSteps
        if minFactor is None:
            self.minFactor = 1.1
        else:
            self.minFactor = minFactor
        # This is initialized later on.
        self.component_properties = {}


    @classmethod
    def from_dict(cls, bd, systems, coordinates, protocols, grid):
        maxSteps, minFactor = cls._ext_from_dict(bd)
        # first create the core protocols
        core_protocols = GmxCoreAlchemicalProtocol.from_dict(bd, coordinates, grid)
        if 'ensemble' in bd.keys():
            ensemble = Ensemble.from_string(bd['ensemble'][0])
        else:
            # adopt NVT in case nothing is said.
            # this should only affect including or not pV.
            ensemble = Ensemble.NVT
        # then create the object
        out = cls(bd['name'][0],
                  core_protocols,
                  bd['mdps'],
                  [],
                  ensemble,
                  maxSteps=maxSteps,
                  minFactor=minFactor)
        out.set_custom_attributes_from_blockdict(bd, systems, coordinates,
                                                 protocols)
        return out

    # overriden
    def run_gridpoint_at_dir(self, gridpoint, workdir, simulate=True):
        out = {}
        for i, p in enumerate(self.core_protocols):
            # replace workdir name
            core_workdir = gridpoint.makeSimudir(
                gridpoint.baseGrid.makeProtocolSimudir(p))
            this_out = p.run_gridpoint_at_dir(gridpoint,
                                              core_workdir,
                                              simulate=simulate)
            for k,v in this_out.items():
                try:
                    out[k].append(v)
                except KeyError:
                    out[k] = [v]
        gridpoint.add_protocol_output(self, out)

    # overriden
    def prepare_gridpoint_at_dir(self, gridpoint, workdir):
        if self.requires_reweight():
            out = {}
            for i, p in enumerate(self.core_protocols):
                # replace workdir name
                core_workdir = gridpoint.makeSimudir(
                    gridpoint.baseGrid.makeProtocolSimudir(p))
                this_out = p.prepare_gridpoint_at_dir(gridpoint, core_workdir)
                for k,v in this_out.items():
                    try:
                        out[k].append(v)
                    except KeyError:
                        out[k] = [v]
            gridpoint.add_protocol_output(self, out)


class GmxProtocol(GmxBaseProtocol,
                  DefaultExtendMixin,
                  CustomizableAttributesMixin):
    """
    A Gromacs-based simulation protocol for general usage.

    
    Attributes
    ----------

    name : string

    system : string

    conf_fac : ConfigurationFactory
    
    ensemble : Ensemble

    mdps : list of str/os.PathLike

    properties : list of str
        Property types (e.g., density, potential, ...).
        Usually initialized with the empty list.

    component_properties : dict of BaseAtomicProperty
        Keys are the types of the properties.
        This is not initialized together with the Protocol.

    surrogate_models: list of (SurrogateModel, str)
        The string is the property type associated with the
        surrogate model (e.g., density, potential, ...).

    ..from DefaultExtendMixin:

    maxSteps : int
        Maximum number of steps of the production run.

    minFactor : float (optional)
        Minimum factor to multiply by the number of steps when
        the production run is extended. If not given, it is set
        to 1.1.
    """
    
    def __init__(self, name, system, conf_fac, mdps, ensemble, maxSteps,
                 minFactor=None):
        """
        Parameters
        ----------
        
        name : string

        system : string

        conf_fac : ConfigurationFactory

        ensemble : Ensemble

        mdps : list of str/os.PathLike

        properties : list of str
            Property types (e.g., density, potential, ...).
            Usually initialized with the empty list.

        maxSteps : int
            Maximum number of steps of the production run.

        minFactor : float (optional)
            Minimum factor to multiply by the number of steps when
            the production run is extended. If not given, it is set
            to 1.1.

        """
        from gmak.read_input import InputError
        
        self.name = name
        self.system = system
        self.conf_fac = conf_fac
        self.ensemble = ensemble
        self.mdps = mdps
        # Set _has_pv just to avoid breaking existing code.
        if self.ensemble == Ensemble.NPT:
            self._has_pv = True
        else:
            self._has_pv = False
        if maxSteps is None:
            raise InputError("GmxProtocol must specify maxsteps.")
        else:
            self.maxSteps = maxSteps
        if minFactor is None:
            self.minFactor = 1.1
        else:
            self.minFactor = minFactor
        # This is initialized later on.
        self.component_properties = {}
        
    @classmethod
    def from_dict(cls, bd, systems, coordinates, protocols):
        maxSteps, minFactor = cls._ext_from_dict(bd)
        if 'ensemble' in bd.keys():
            ensemble = Ensemble.from_string(bd['ensemble'][0])
        else:
            # adopt NVT in case nothing is said.
            # this should only affect including or not pV.
            ensemble = Ensemble.NVT
        out  = cls(bd["name"][0],
                   bd["system"][0],
                   coordinates[bd["coords"][0]],
                   bd["mdps"],
                   ensemble,
                   maxSteps,
                   minFactor)
        # set custom attributes
        out.set_custom_attributes_from_blockdict(bd, systems, coordinates,
                                                 protocols)
        return out

    
    def run_gridpoint_at_dir(self, gridpoint, workdir):
        # labels are 0, 1, 2, ..., NUMBER_OF_MDPS
        labels = [str(x) for x in range(len(self.mdps))]
        # number of steps---either from the last .mdp or from extension
        nsteps = gridpoint.getProtocolLength(self)
        # create initial configuration
        conf = self.conf_fac.construct_configuration(
            os.path.join(workdir, f"{self.name}.gro"))
        # get topology
        topo = gridpoint.getTopologyPath(self.system)
        # perform simulations
        out = simulate_protocol_general(conf.get_path(),
                                        topo,
                                        self.mdps,
                                        nsteps,
                                        labels,
                                        workdir)
        gridpoint.add_protocol_output(self, out)

        
    def prepare_gridpoint_at_dir(self, gridpoint, workdir):
        # When there is reweighting, store path of the topology
        # TODO: I believe this is no longer needed.
        if self.requires_reweight():
            topo = gridpoint.getTopologyPath(self.system),
            gridpoint.add_protocol_output(
                self,
                {'top': topo})


class CustomProtocol(BaseProtocol,
                     CustomizableAttributesMixin):

    @staticmethod
    def _default_calc_initial_len(args_dict):
        return None

    @staticmethod
    def _default_calc_extend(errs_tols, length):
        return None

    @staticmethod
    def _default_get_last_frame(protocol_output):
        raise NotImplementedError("Can't get last frame for custom protocol.")

    def __init__(self, simulator, calc_initial_len=None, calc_extend=None,
                 get_last_frame=None):
        self.properties = []
        # This is initialized later on.
        self.component_properties = {}
        self.surrogate_models = []
        self.simulator = simulator
        if calc_initial_len is None:
            self._calc_initial_len = self._default_calc_initial_len
        else:
            self._calc_initial_len = calc_initial_len
        if calc_extend is None:
            self._calc_extend = self._default_calc_extend
        else:
            self._calc_extend = calc_extend
        if get_last_frame is None:
            self._get_last_frame = self._default_get_last_frame
        else:
            self._get_last_frame = get_last_frame

    def calc_initial_len(self):
        """
        Wrapper for self._calc_initial_len.
        """
        return self._calc_initial_len(self.get_custom_attributes())

    def calc_extend(self, gridpoint, optimizer, protocolsHash):
        """
        Wrapper for self._calc_extend to be more easily used in gridbase.
        """
        errs_tols = gridpoint.get_errs_tols(optimizer, self, protocolsHash)
        length = gridpoint.getProtocolLength(self)
        return self._calc_extend(errs_tols, length,
                                 self.get_custom_attributes())

    def get_last_frame(self, gridpoint):
        protocol_output = gridpoint.protocol_outputs[self.name]
        return ConfigurationFactory.from_file(self._get_last_frame(protocol_output))


    @classmethod
    def from_dict(cls, bd, systems, coordinates, protocols):
        # 'name', 'system', 'coords' and 'type' are required keys
        for key in ['name', 'system', 'coords', 'type']:
            if (key not in bd.keys()):
                raise Exception(f"CustomProtocol must specify a '{key}'.")
        if 'ensemble' in bd.keys():
            ensemble = Ensemble.from_string(bd['ensemble'][0])
        else:
            # adopt NVT in case nothing is said.
            # this should only affect including or not pV.
            ensemble = Ensemble.NVT
        out = CustomProtocolFactory.create(bd['type'][0])
        out.name = bd['name'][0]
        out.ensemble = ensemble
        if ensemble == Ensemble.NPT:
            out._has_pv = True
        else:
            out._has_pv = False
        out.system = bd['system'][0]
        out.type = bd['type'][0]
        out.name = bd['name'][0]
        out.coords = coordinates[bd['coords'][0]]
        out.set_custom_attributes_from_blockdict(bd, systems, coordinates,
                                                 protocols)
        return out

    def run_gridpoint_at_dir(self, gridpoint, workdir):
        initial_length = self.calc_initial_len()
        length = gridpoint.getProtocolLength(self)
        conf = self.coords.construct_configuration(
            os.path.join(workdir, "start.gro"))
        out = self.simulator(
            length,
            gridpoint.getTopologyPath(self.system),
            conf.get_path(),
            initial_length != length,
            self.get_custom_attributes(),
            workdir)
        gridpoint.add_protocol_output(self, out)

    def prepare_gridpoint_at_dir(self, gridpoint, workdir):
        # TODO: If one wants to rw, it would probably be easier to
        # expand the simulator to also encompass cases where a
        # simulate flag is set to False, and then use this here like
        #
        # if self.requires_reweight():
        #     self.simulator(..., simulate=False)
        raise NotImplementedError("Reweighting not implemented for "
                                  "custom protocols.")


class CustomProtocolFactory:

    # this is a dict 'type' : simulator
    ptable = {}

    @classmethod
    def add_custom_protocol(cls, type_name, simulator,
                            calc_initial_len=None, calc_extend=None,
                            get_last_frame=None):
        cls.ptable[type_name] = (simulator, calc_initial_len, calc_extend,
                                 get_last_frame)

    @classmethod
    def create(cls, type_name):
        return CustomProtocol(cls.ptable[type_name][0],
                              calc_initial_len=cls.ptable[type_name][1],
                              calc_extend=cls.ptable[type_name][2],
                              get_last_frame=cls.ptable[type_name][3])


def create_protocols(bd, protocols, coordinates, systems, grid):
    _type = bd['type'][0]
    if _type == 'gmx':
        return GmxProtocol.from_dict(bd, systems, coordinates, protocols)
    elif _type == 'gmx_alchemical':
        out = GmxAlchemicalProtocol.from_dict(bd, systems, coordinates,
                                              protocols, grid)
        # TODO think of a nicer way to set parent
        for p in out.core_protocols:
            p.parent = out
        return out
    else:
        # Custom
        try:
            return CustomProtocol.from_dict(bd, systems, coordinates,
                                            protocols)
        except KeyError:
            raise NotImplementedError(f"Can't create protocol of type {_type}.")
