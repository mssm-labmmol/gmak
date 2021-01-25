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

class BaseProtocol:
    """Contains few methods where implementation is common to all protocols."""
    def get_filtering_properties (self):
        standard_properties = self.get_properties()
        if ('potential' not in standard_properties) and (self.get_mbar_model() is not None):
            standard_properties.append('potential')
        if self.has_pv():
            if ('pV' not in standard_properties) and (self.get_mbar_model() is not None):
                standard_properties.append('pV')
        return standard_properties

    def expand(self):
        return [self,]
    
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

class SlabProtocol(BaseProtocol):

    def __init__ (self):
        self.name = ""
        self.type = "slab"
        self.mdps = []
        self.box_size = 0.0
        self.nprocs = -1

    def __init__ (self, name, mdps, factor, properties, nprocs):
        self.name = name
        self.type = "slab"
        self.mdps = mdps
        self.properties = properties
        self.nprocs = int(nprocs)

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
            if line.split()[0] == 'follow':
                self.follow = line.split()[1]
            if line.split()[0] == 'nprocs':
                self.nprocs = int(line.split()[1])

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
        labels   =  [str(x) for x in range(len(self.mdps))]
        conf     =  gridpoint.protocol_outputs[self.follow]['gro']
        top      =  gridpoint.protocol_outputs[self.follow]['top']
        out_slab = dummy_protocol_slab (conf, top, self.mdps, labels,\
                workdir)
        gridpoint.add_protocol_output (self, out_slab)
        return

    def has_pv (self):
        return False

class LiquidProtocol(BaseProtocol):

    def __init__ (self):
        self.name = ""
        self.type = "liquid"
        self.mdps = []
        self.mdps_rw = []
        self.properties = []
        self.coords = ""
        self.box_size = 0.0

    def __init__ (self, name, nmols, coords, mdps, box, properties):
        self.name = name
        self.type = "liquid"
        self.mdps = mdps
        self.coords = coords
        self.box_size = box
        self.nmols = nmols
        self.properties = properties

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
            if line.split()[0] == 'length':
                box_len = float(line.split()[1])
                self.box_size = [box_len, box_len, box_len]
            if line.split()[0] == 'nmols':
                self.nmols = int(line.split()[1])

    def has_pv (self):
        return True

    def run_gridpoint_at_dir (self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        nsteps = gridpoint.getProtocolSteps(self)
        out_liquid = simulate_protocol_liquid (self.coords, self.nmols,\
                                               self.box_size, gridpoint.getTopologyPath(self.molecule), self.mdps, nsteps, labels, workdir)
        gridpoint.add_protocol_output (self, out_liquid)
        return

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        out_liquid = dummy_protocol_liquid (self.coords, self.nmols, \
                self.box_size, gridpoint.getTopologyPath(self.molecule), self.mdps, labels, workdir)
        gridpoint.add_protocol_output (self, out_liquid)
        return
    
class GasProtocol(BaseProtocol):

    def __init__ (self):
        self.name = ""
        self.type = "gas"
        self.mdps = []
        self.coords = ""
        self.box_size = 0.0

    def __init__ (self, name, coords, dipole, polar, mdps, properties):
        self.name = name
        self.type = "gas"
        self.mdps = mdps
        self.coords = coords
        self.nmols = 1
        self.properties = properties
        self.dipole = dipole
        self.polar = polar

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
    
    def has_pv (self):
        # NOTE: Formally, the gas simulations do have a pV term. However, they
        # are the same in all cases (pV = NkT), so they have no influence over
        # the distribution in phase space.
        return False

    def run_gridpoint_at_dir (self, gridpoint, workdir):
        labels = [str(x) for x in range(len(self.mdps))]
        nsteps = gridpoint.getProtocolSteps(self)
        out_gas = simulate_protocol_gas (self.coords, gridpoint.getTopologyPath(self.molecule),\
                                         self.mdps, nsteps, labels, workdir)
        gridpoint.add_protocol_output (self, out_gas)
        return

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
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

class SolvationProtocol(BaseProtocol):

    type = "solv"

    def __init__ (self):
        self.name = ""
        self.mdps = []
        self.properties = []
        self.coords = {'solvent': '', 'solute': ''}
        self.molecules = {'solvent': '', 'solute': ''}
        self.nmols  = {'solvent': 0, 'solute': 0}
        self.box_size = 0.0
        self.follows = None

    def __init__ (self, name, moleculeSolvent, moleculeSolute,
                  nmolsSolvent, nmolsSolute, coordsSolvent, coordsSolute, mdps, 
                  properties, followsObj):
        self.name = name
        self.mdps = mdps
        self.molecules = {'solvent': moleculeSolvent, 'solute': moleculeSolute}
        self.coords = {'solvent': coordsSolvent, 'solute': coordsSolute}
        self.nmols = {'solvent': nmolsSolvent, 'solute': nmolsSolute}
        self.properties = properties
        self.follows= followsObj

    def has_pv (self):
        return True

    def run_gridpoint_at_dir (self, gridpoint, workdir, simulate=True):
        labels = [str(x) for x in range(len(self.mdps))]
        nsteps = gridpoint.getProtocolSteps(self)
        if self.follows is not None and simulate:
            runcmd.run("mkdir -p {}".format(workdir))
            copyfile(gridpoint.protocol_outputs[self.follows.name]['gro'],
                     workdir + '/solvation.gro')
        out = simulate_protocol_solvation(self.coords, self.nmols,
                                          {'solute':
                                           gridpoint.getTopologyPath(self.molecules['solute']),
                                           'solvent':
                                           gridpoint.getTopologyPath(self.molecules['solvent'])},
                                          self.mdps, nsteps, labels,
                                          workdir, simulate,
                                          self.follows is not None)
        gridpoint.add_protocol_output (self, out)

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        self.run_gridpoint_at_dir(gridpoint, workdir, simulate=False)


class DummyProtocol(BaseProtocol):

    def __init__ (self, name, protocol):
        self.type = 'dgsolv'
        self.name = name
        self.properties = []
        self.mdps = protocol.mdps
        self.coords = protocol.coords
        self.molecules = protocol.molecules
        self.nmols  = {'solvent': 0, 'solute': 0}
        self.box_size = 0.0
        self.follows = None
        self.inner_protocols = []

    def point_to(self, protocols):
        self.inner_protocols = protocols

    def expand(self):
        return self.inner_protocols

    def has_pv (self):
        return False

    def run_gridpoint_at_dir (self, gridpoint, workdir):
        return

    def prepare_gridpoint_at_dir (self, gridpoint, workdir):
        return

class SolvationProtocolFactory(ABC):
    @abstractmethod
    def createSolvationProtocolFromDict(self, args):
        pass

class SolvationFreeEnergyFactory(SolvationProtocolFactory):
    def __init__ (self):
        return

    def createSolvationProtocolFromDict(self, args):
        # read number of lambda points from some mdp file
        imu = mdpUtils()
        imu.parse_file(args['mdps'][-1])
        nlambdas = imu.get_nlambdas()
        if (int(imu.get_option('calc-lambda-neighbors')[0]) != -1):
            warnings.warn('calc-lambda-neighbors will be set to -1')
        outputProtocols = []
        for i in range(nlambdas):
            mdps = []
            for m in args['mdps']:
                prefix = os.path.abspath(m)[:-4]
                fn = prefix + '_{}.mdp'.format(i)
                mdps.append(fn)
                imu.clean()
                imu.parse_file(m)
                imu.set_lambda_state(i)
                imu.set_option('calc-lambda-neighbors', ['-1'])
                imu.write_to_file(fn)
            # create
            if (i == 0):
                follows = None
            else:
                follows = outputProtocols[i - 1]
            outputProtocols.append(SolvationProtocol(args['name'][0] +
                                                     '-{}'.format(i),
                                                     args['solvent'][0],
                                                     args['solute'][0],
                                                     int(args['solvent'][1]),
                                                     int(args['solute'][1]),
                                                     args['solvent'][2],
                                                     args['solute'][2],
                                                     mdps,
                                                     [], follows))
        return outputProtocols
