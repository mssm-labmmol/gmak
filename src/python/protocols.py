from abc import ABC, abstractmethod
import runcmd
from simulate import *
from shutil import copyfile
from surrogate_model import *
import warnings
import re
import copy
import os
from logger import *
from mdputils import *
import atomic_properties

class BaseProtocol:
    """Contains few methods where implementation is common to all protocols."""
    def get_filtering_properties(self):
        raw_properties = self.get_properties()
        # filter those that are timeseries
        standard_properties = []
        for p in raw_properties:
            if p == 'polcorr':
                mu = self.dipole
                alpha = self.polar
                ap = atomic_properties.create_atomic_property(p, mu, alpha)
            elif p == 'dgsolv':
                temp = self.get_temperature()
                ap = atomic_properties.create_atomic_property(p, temp)
            else:
                ap = atomic_properties.create_atomic_property(p)
            if ap.is_timeseries:
                standard_properties.append(p)
        if ('potential' not in standard_properties) and (self.get_mbar_model() is not None):
            standard_properties.append('potential')
        if self.has_pv():
            if ('pV' not in standard_properties) and (self.get_mbar_model() is not None):
                standard_properties.append('pV')
        return standard_properties

    def get_nonfiltering_properties(self):
        raw_properties = self.get_properties()
        filter_properties = self.get_filtering_properties()
        return [p for p in raw_properties if p not in filter_properties]

    def get_properties (self):
        return self.properties

    def get_temperature (self):
        temp = 0.0
        fp = open(self.mdps[-1], "r")
        for line in fp:
            if ((re.match(r"^ref_t.*", line)) or (re.match(r"^ref-t.*", line))):
                temp = float(line.split()[2])
        fp.close()
        return temp

    def add_surrogate_model (self, surrogate_model_string, atomic_property, bool_legacy):
        """Add a (surrogate_model, property) tuple to the protocol."""
        model = init_surrogate_model_from_string(surrogate_model_string, bool_legacy)
        try:
            self.surrogate_models.append((model, atomic_property))
        except AttributeError:
            # this probably means the list has not been initialized
            self.surrogate_models = []
            self.surrogate_models.append((model, atomic_property))

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

    def get_avg_err_estimate_of_property (self, prop, kind):
        """Returns a tuple (EA_k, dEA_k) for atomic property 'prop' estimated via 'kind' surrogate model for all states."""
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


class LiquidProtocol(BaseProtocol):

    type = 'liquid'
    _has_pv = True

    def __init__ (self, name, molecule, nmols, coords, mdps, box, properties):
        self.name = name
        self.molecule = molecule
        self.mdps = mdps
        self.coords = coords
        self.box_size = box
        self.nmols = nmols
        self.properties = properties
        self.surrogate_models = []

    @classmethod
    def from_dict(cls, bd):
        box_len = float(bd['length'][0])
        box     = [box_len, box_len, box_len]
        return cls(bd['name'][0], bd['molecule'][0], int(bd['nmols'][0]),
                   bd['coords'][0], bd['mdps'], box, [])

    def run_gridpoint_at_dir (self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        nsteps = gridpoint.getProtocolSteps(self)
        out_liquid = simulate_protocol_liquid (self.coords,
                                               self.nmols,
                                               self.box_size,
                                               gridpoint.getTopologyPath(self.molecule),
                                               self.mdps, nsteps, labels, workdir)
        gridpoint.add_protocol_output (self, out_liquid)
        return

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        if self.requires_reweight():
            labels = [str(x) for x in range(len(self.mdps))]
            out_liquid = dummy_protocol_liquid (self.coords, self.nmols, \
                    self.box_size, gridpoint.getTopologyPath(self.molecule), self.mdps, labels, workdir)
            gridpoint.add_protocol_output (self, out_liquid)
            return


class GasProtocol(BaseProtocol):

    type = 'gas'
    _has_pv = False

    def __init__ (self, name, molecule, mdps, coords, dipole, polar, properties):
        self.name = name
        self.molecule = molecule
        self.mdps = mdps
        self.coords = coords
        self.dipole = dipole
        self.polar = polar
        self.properties = properties
        self.nmols = 1
        self.surrogate_models = []

    def read_from_stream (self, stream):
        for line in stream:
            if line[0] == '#':
                continue
            if line.rstrip() == '$end':
                return
            if line.split()[0] == 'name':
                self.name = line.split()[1]
            if line.split()[0] == 'molecule':
                self.molecule = line.split()[1]
            if line.split()[0] == 'mdps':
                self.mdps = line.split()[1:]
            if line.split()[0] == 'coords':
                self.coords = line.split()[1]
            if line.split()[0] == 'gasdipole':
                self.dipole = float(line.split()[1])
            if line.split()[0] == 'polarizability':
                self.polar = float(line.split()[1])

    @classmethod
    def from_dict(cls, bd):
        return cls(bd['name'][0], bd['molecule'][0], bd['mdps'], bd['coords'][0],
                   float(bd['gasdipole'][0]), float(bd['polarizability'][0]), [])

    def run_gridpoint_at_dir (self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        nsteps = gridpoint.getProtocolSteps(self)
        out_gas = simulate_protocol_gas (self.coords, gridpoint.getTopologyPath(self.molecule),\
                                         self.mdps, nsteps, labels, workdir)
        gridpoint.add_protocol_output (self, out_gas)
        return

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        if self.requires_reweight():
            labels = [str(x) for x in range(len(self.mdps))]
            out_gas = dummy_protocol_gas (self.coords, gridpoint.getTopologyPath(self.molecule),\
                    self.mdps, labels, workdir)
            gridpoint.add_protocol_output (self, out_gas)
            return

    def set_other_corrections (self, corr):
        self.other_corrections = corr

    def get_avg_err_estimate_of_property (self, prop, kind):
        """Override when polcorr is not necessary."""
        if (prop == 'polcorr') and (self.polar < 0):
            EA, dEA = super().get_avg_err_estimate_of_property('potential', kind)
            EA = copy.deepcopy(EA)
            dEA = copy.deepcopy(dEA)
            for i in range(len(EA)):
                EA[i] = 0.0
                dEA[i] = 0.0
            return EA, dEA
        else:
            EA, dEA = super().get_avg_err_estimate_of_property(prop, kind)
            return EA, dEA


class SlabProtocol(BaseProtocol):

    type = "slab"
    _has_pv = False

    def __init__ (self, name, molecule, mdps, follow, properties, nprocs=-1):
        self.name = name
        self.molecule = molecule
        self.mdps = mdps
        self.properties = properties
        self.nprocs = nprocs
        self.follow = None
        self.surrogate_models = []

    @classmethod
    def from_dict(cls, bd):
        if 'nprocs' in bd.keys():
            return cls(bd['name'][0], bd['molecule'][0], bd['mdps'],
                       bd['follow'], [], int(bd['nprocs']))
        else:
            return cls(bd['name'][0], bd['molecule'][0], bd['mdps'],
                       bd['follow'], [], int(bd['nprocs']))


    def set_follow (self, pro_name):
        self.follow = pro_name

    def run_gridpoint_at_dir (self, gridpoint, workdir):
        labels   =  [str(x) for x in range(len(self.mdps))]
        conf     =  gridpoint.protocol_outputs[self.follow]['gro']
        top      =  gridpoint.protocol_outputs[self.follow]['top']
        liq_tpr  =  gridpoint.protocol_outputs[self.follow]['tpr']
        nsteps   =  gridpoint.getProtocolSteps(self)
        out_slab = simulate_protocol_slab (conf, top, liq_tpr, self.mdps, nsteps, labels, workdir, self.nprocs)
        gridpoint.add_protocol_output (self, out_slab)
        return

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        if self.requires_reweight():
            labels   =  [str(x) for x in range(len(self.mdps))]
            conf     =  gridpoint.protocol_outputs[self.follow]['gro']
            top      =  gridpoint.protocol_outputs[self.follow]['top']
            out_slab = dummy_protocol_slab (conf, top, self.mdps, labels,\
                    workdir)
            gridpoint.add_protocol_output (self, out_slab)
            return


class CoreSolvationProtocol(BaseProtocol):

    type = "_coresolv"
    _has_pv = True

    def __init__ (self, name, moleculeSolvent, moleculeSolute,
                  nmolsSolvent, nmolsSolute, coordsSolvent, coordsSolute, mdps, 
                  properties, followsObj):
        self.name = name
        self.mdps = mdps
        self.molecules = {'solvent': moleculeSolvent, 'solute': moleculeSolute}
        self.coords = {'solvent': coordsSolvent, 'solute': coordsSolute}
        self.nmols = {'solvent': nmolsSolvent, 'solute': nmolsSolute}
        self.properties = properties
        self.follows = followsObj
        self.surrogate_models = []
        self.parent = None
        self.output_files = None

    def run_gridpoint_at_dir (self, gridpoint, workdir, simulate=True):
        labels = [str(x) for x in range(len(self.mdps))]
        nsteps = gridpoint.getProtocolSteps(self.parent)
        if self.follows is not None and simulate:
            runcmd.run("mkdir -p {}".format(workdir))
            copyfile(self.follows.output_files['gro'],
                     workdir + '/solvation.gro')
        out = simulate_protocol_solvation(
            self.coords, self.nmols,
            {'solute':
             gridpoint.getTopologyPath(self.molecules['solute']),
             'solvent':
             gridpoint.getTopologyPath(self.molecules['solvent'])},
            self.mdps, nsteps, labels,
            workdir, simulate,
            self.follows is not None)
        self.output_files = out
        return out

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        if self.requires_reweight():
            return self.run_gridpoint_at_dir(gridpoint, workdir, simulate=False)

    @classmethod
    def from_dict(cls, bd):
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
                follows = None
            else:
                follows = outputProtocols[i - 1]
            # create SolvationProtocol for each state and append it to the
            # list of states returned by this factory
            outputProtocols.append(
                cls(bd['name'][0] + '-{}'.format(i),
                    bd['solvent'][0],
                    bd['solute'][0],
                    int(bd['solvent'][1]),
                    int(bd['solute'][1]),
                    bd['solvent'][2],
                    bd['solute'][2],
                    mdps,
                    [],
                    follows))
        return outputProtocols


class SolvationProtocol(BaseProtocol):

    type = "solv"
    _has_pv = True

    def __init__(self, name, core_protocols, mdps, properties):
        self.name = name
        self.core_protocols = core_protocols
        self.mdps = mdps
        self.properties = properties
        self.surrogate_models = []

    @classmethod
    def from_dict(cls, bd):
        # first create the core protocols
        core_protocols = CoreSolvationProtocol.from_dict(bd)
        # then create the object
        return cls(bd['name'][0],
                   core_protocols,
                   bd['mdps'],
                   [])

    # overriden
    def run_gridpoint_at_dir(self, gridpoint, workdir, simulate=True):
        pre_out = []
        for i, p in enumerate(self.core_protocols):
            # replace workdir name
            core_workdir = gridpoint.makeSimudir(
                gridpoint.baseGrid.makeProtocolSimudir(p))
            pre_out.append(p.run_gridpoint_at_dir(gridpoint,
                                              core_workdir,
                                              simulate=simulate))
        # pre_out is a list of dicts, all of them with the same keys.
        # we have to convert it into the corresponding dict of lists.
        keys = pre_out[0].keys()
        out  = {}
        for k in keys:
            out[k] = []
            for o in pre_out:
                out[k].append(o[k])
        gridpoint.add_protocol_output(self, out)

    # overriden
    def prepare_gridpoint_at_dir(self, gridpoint, workdir):
        if self.requires_reweight():
            pre_out = []
            for i, p in enumerate(self.core_protocols):
                # replace workdir name
                core_workdir = gridpoint.makeSimudir(
                    gridpoint.baseGrid.makeProtocolSimudir(p))
                pre_out.append(p.prepare_gridpoint_at_dir(gridpoint, core_workdir))
            # pre_out is a list of dicts, all of them with the same keys.
            # we have to convert it into the corresponding dict of lists.
            keys = pre_out[0].keys()
            out  = {}
            for k in keys:
                out[k] = []
                for o in pre_out:
                    out[k].append(o[k])
            gridpoint.add_protocol_output(self, out)


class GeneralProtocol(BaseProtocol):

    type = "general"
    _has_pv = False

    def __init__(self, name, system, coords, mdps, properties):
        self.name = name
        self.coords = coords
        self.system = system
        self.mdps = mdps
        self.properties = properties
        self.surrogate_models = []

    @classmethod
    def from_dict(cls, bd):
        return cls(bd["name"][0], bd["system"][0], bd["coords"][0], bd["mdps"], [])

    def run_gridpoint_at_dir(self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        nsteps = gridpoint.getProtocolSteps(self)
        out    = simulate_protocol_general(self.coords,
                                           gridpoint.getTopologyPath(self.system),
                                           self.mdps,
                                           nsteps,
                                           labels,
                                           workdir)
        gridpoint.add_protocol_output(self, out)

    def prepare_gridpoint_at_dir(self, gridpoint, workdir):
        if self.requires_reweight():
            labels = [str(x) for x in range(len(self.mdps))]
            nsteps = gridpoint.getProtocolSteps(self)
            out    = dummy_protocol_general(self.coords,
                                            gridpoint.getTopologyPath(self.system),
                                            self.mdps,
                                            labels,
                                            workdir)
            gridpoint.add_protocol_output(self, out)


class CustomProtocol(BaseProtocol):
    _has_pv = False

    def __init__(self, simulator):
        # __init__ does almost nothing because attributes are mostly set by
        # cls.from_dict
        self.properties = []
        self.surrogate_models = []
        self.simulator = simulator

    @classmethod
    def from_dict(cls, bd):
        # 'name', 'system' and 'type' are a required keys
        for key in ['name', 'system', 'type']:
            if (key not in bd.keys()):
                raise Exception(f"CustomProtocol must specify a '{key}'.")
        out = CustomProtocolFactory.create(bd['type'][0])
        for k in bd.keys():
            if len(bd[k]) == 1:
                setattr(out, k, bd[k][0])
            else:
                setattr(out, k, bd[k])
        return out

    def run_gridpoint_at_dir(self, gridpoint, workdir):
        out = self.simulator(gridpoint.getTopologyPath(self.system),
                             vars(self),
                             workdir)
        gridpoint.add_protocol_output(self, out)

    def prepare_gridpoint_at_dir(self, gridpoint, workdir):
        return


class CustomProtocolFactory:

    # this is a dict 'type' : simulator
    ptable = {}

    @classmethod
    def add_custom_protocol(cls, type_name, simulator):
        cls.ptable[type_name] = simulator

    @classmethod
    def create(cls, type_name):
        return CustomProtocol(cls.ptable[type_name])

def add_custom_protocol(type_name, simulator):
    CustomProtocolFactory.add_custom_protocol(type_name, simulator)

def create_protocols(bd):
    _type = bd['type'][0]
    if _type == 'liquid':
        return [LiquidProtocol.from_dict(bd),]
    elif _type == 'gas':
        return [GasProtocol.from_dict(bd),]
    elif _type == 'slab':
        return [SlabProtocol.from_dict(bd),]
    elif _type == 'general':
        return [GeneralProtocol.from_dict(bd),]
    elif _type == 'solv':
        out = SolvationProtocol.from_dict(bd)
        # TODO think of a nicer way to set parent
        for p in out.core_protocols:
            p.parent = out
        return [out,]
    else:
        # Custom
        try:
            return [CustomProtocol.from_dict(bd),]
        except KeyError:
            raise NotImplementedError(f"Can't create protocol of type {_type}.")

