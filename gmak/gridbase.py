import gmak.runcmd as runcmd
import sys
import re
import os
import shlex
import numpy           as     np
from gmak.property        import DGAlchemicalAnalysis
from gmak.mdputils        import *
from gmak.reweightbase    import *
from gmak.traj_ana        import *
from gmak.traj_filter     import *
from gmak.reweight        import *
import multiprocessing
from gmak.property        import *
from gmak.grid_ana        import *
import warnings
from gmak.surrogate_model import *
from gmak.cartesiangrid   import *
from gmak.gridshifter     import *
import copy
import pickle
import gmak.logger as logger
from gmak.state           import *
import gmak.component_properties as component_properties

def replaceMacros(fn, macros):
    #print("Replacing macros in %s:" % fn)
    fp = open(fn, 'r')
    data = fp.read()
    for macro in macros:
        token = macro.token
        value = macro.value
        #print("\t%s -> %s" % (token, value))
        data = data.replace(token, value)
    fp.close()
    fp = open(fn, 'w')
    fp.write(data)
    fp.close()

def reweightItpChanger (inputTopology, index, outputTopology):
    fp = open(inputTopology, 'r')
    fp_out = open(outputTopology, 'w')
    for line in fp:
        if (line[:8] == '#include'):
            pathToFile = line[line.find('\"')+1:line.rfind('\"')]
            newPathToFile = pathToFile[:pathToFile.find('.itp')] + ("_%d.itp" % (index))
            newLine = "#include \"%s\"\n" % newPathToFile
            fp_out.write(newLine)
        else:
            fp_out.write(line)
    fp_out.close()
    fp.close()

# This function must have global scope for Pool.map.
def core_parallel_clean_reweight_i_at_j (aux_input_list):
    i             =  aux_input_list[0]
    gp            =  aux_input_list[1]
    j             =  aux_input_list[2]
    gp_other      =  aux_input_list[3]
    protocol      =  aux_input_list[4]
    workdir       =  aux_input_list[5]
    workdir_mbar  =  aux_input_list[6]
    properties    =  aux_input_list[7]

    gp.reweight_with_protocol_at_other (protocol, gp_other,
            "%s/%d_at_%d/" % (workdir, i, j))
    # now get properties and clean
    gi = gp 
    gi_id = i
    gj_id = j 
    # get potential (and) pV
    out_file = workdir_mbar + \
            "/reweighted_properties/potential_%d_%d.xvg" % \
            (gi.id, gj_id)
    xtc = gi.rw_outputs[gj_id][protocol.name]['xtc']
    edr = gi.rw_outputs[gj_id][protocol.name]['edr']
    tpr = gi.rw_outputs[gj_id][protocol.name]['tpr']
    trr = gi.rw_outputs[gj_id][protocol.name]['trr']
    gro = gi.rw_outputs[gj_id][protocol.name]['gro']
    # If file already exists, no need to create it again.
    if not (os.path.isfile(out_file)):
        # get
        obtain_property (xtc, edr, gro, tpr, 'potential', out_file)
    if (protocol.has_pv()):
        out_file = workdir_mbar + \
            "/reweighted_properties/pV_%d_%d.xvg" % \
            (gi.id, gj_id)
        # If file already exists, no need to create it again, just set the variable.
        if not (os.path.isfile(out_file)):
            # copy from original trajectory
            original_file = gi.retrieve_component_property_from_protocol('pV', protocol)            
            runcmd.run("cp {} {}".format(original_file, out_file))
    # Delete trajectory files (the ones that are big).
    # We don't need to delete the xtc file because it corresponds to the original trajectory.
    for ff in [trr,edr]:
        if (os.path.isfile(ff)):
            os.remove (ff)

class GridPoint:
    """Class responsible for attributes regarding a single grid point.

    members:

        baseGrid: A reference to the ParameterGrid.

        id: A number identifying this gridpoint. This is necessary to
            use gridpoints as keys in dictionaries.

        is_sample: Flag which determines if this grid point is in
            samples list. Note that we assume all protocols for a grid
            involve the same samples.

        is_simulated: List of protocol names where the point has been
            simulated.

        protocol_outputs: Dictionary containing final output files of
            the protocols executed.

        topology_outputs: Dictionary containing the topology output
            objects for each system

        component_properties 2-d dictionary containing output files for
            the properties extracted,
            e.g. component_properties['liquid']['potential'] is the file
            containing potential energies for protocol liquid.

        rw_outputs: A dictionary where keys are other gridpoints and
            values are dictionaries indexing the output files of the
            reweighting. For instance, rw_outputs[gp][pr]['xtc']
            accesses the xtc file of the reweighting on grid point
            'gp' simulated under protocol 'pr'. This may be empty, if
            no reweighting is explicitly performed.

        estimated_properties: dictionary of estimated_properties
            e.g. estimated_propeties['gamma'] should be a
            PropertyBase::Gamma

    """
    def __init__ (self, baseGrid, idx):
        self.baseGrid = baseGrid
        self.id = idx
        self.is_sample = False
        self.is_simulated = []
        self.protocol_outputs = {}
        self.topology_outputs = {}
        self.component_properties = {}
        # self.rw_outputs = {}
        # self.rw_outputs[self.id] = {}
        self.estimated_properties = {}
        self.protocol_lenspecs = {}

    def resetWithNewId(self, newId):
        self.id = newId
        self.estimated_properties = {}

    def wasSimulatedWithProtocol(self, protocol):
        return (protocol.name in self.is_simulated)

    def setProtocolAsSimulated(self, protocol):
        if (protocol.name not in self.is_simulated):
            self.is_simulated.append(protocol.name)

    def unsetProtocolAsSimulated(self, protocol):
        if (protocol.name in self.is_simulated):
            self.is_simulated.remove(protocol.name)

    def getTopologyPath(self, molecule):
        return self.topology_outputs[molecule]

    def get_property_estimate(self, prop):
        return self.estimated_properties[prop].value

    def get_property_err(self, prop):
        return self.estimated_properties[prop].err

    def prepare_with_protocol_at_dir (self, protocol, workdir):
        if (protocol.requires_reweight()):
            protocol.prepare_gridpoint_at_dir (self, workdir)

    def simulate_with_protocol_at_dir (self, protocol, workdir):
        protocol.run_gridpoint_at_dir (self, workdir)

    def reweight_with_protocol_at_other (self, protocol, gp_other, workdir):
        # from reweight.py
        raise NotImplementedError("Reweighting not implemented.")
        #if (protocol.type == 'slab'):
        #    out_rw = reweight (self.protocol_outputs[protocol.name]['trr'],\
        #            self.protocol_outputs[protocol.name]['gro'],\
        #            gp_other.protocol_outputs[protocol.name]['top'],\
        #            protocol.mdps[-1], workdir)
        #else:
        #    out_rw = reweight (self.protocol_outputs[protocol.name]['xtc'],\
        #            self.protocol_outputs[protocol.name]['gro'],\
        #            gp_other.protocol_outputs[protocol.name]['top'],\
        #            protocol.mdps[-1], workdir)
        #self.add_reweight_output (gp_other, protocol, out_rw)

    def set_as_sample (self):
        self.is_sample = True

    def is_sample (self):
        return self.is_sample

    # obj is a PropertyBase::X object (X=Density,Gamma,dHvap)
    def add_property_estimate(self, prop_id, prop_name, obj):
        if (prop_name != obj.name):
            raise ValueError ("ERROR: expected {} but got {}.\n".format(prop_name, obj.name))
        self.estimated_properties[prop_id] = obj

    # This is to ensure that any method that receives a grid point has access
    # to the file names of the simulations performed.
    #
    # protocol_outputs is a dictionary for the output (which are also
    # dictionaries) of a simulation, where keys are the protocol names.
    #
    def add_protocol_output(self, protocol, output):
        self.protocol_outputs[protocol.name] = output

    def add_reweight_output (self, gp_other, protocol, output):
        try:
            self.rw_outputs[gp_other.id][protocol.name] = output
        # key error means something has not been initialized
        except KeyError:
            self.rw_outputs[gp_other.id] = {}
            self.rw_outputs[gp_other.id][protocol.name] = output

    def get_component_property_from_protocol(self, name, protocol, output):
        output = os.path.abspath(output)
        ap = protocol.get_component_properties()[name]
        topo = self.topology_outputs[protocol.system]
        ap.calc(self.protocol_outputs[protocol.name], output, topology=topo)
        self.add_component_property_output(name, protocol, output)

    def retrieve_component_property_from_protocol(self, propname, protocol):
        return self.component_properties[protocol.name][propname]

    def add_component_property_output_from_name(self, prop, protocolName, output):
        try:
            self.component_properties[protocolName][prop] = output
        except KeyError:
            # this means component_properties[protocol-name] has not been initialized
            self.component_properties[protocolName] = {}
            self.component_properties[protocolName][prop] = output

    def add_component_property_output (self, prop, protocol, output):
        self.add_component_property_output_from_name(prop, protocol.name, output)


    def makeSimudir(self, protocolSimudir, makeDir=True):
        dirname = os.path.abspath(f"{protocolSimudir}/{self.id}")
        if makeDir:
            try:
                os.mkdir(dirname)
            except FileExistsError:
                pass
        return dirname


    def get_errs_tols(self, optimizer, protocol, protocolsHash):
        """Returns a dictionary <property_name:str, {'tol': float, 'err': float}>.

        'tol' is the tolerance for the property specified in the input file.
        'err' is the actual error obtained in the estimation.
        """
        out_dict   = {}
        properties = optimizer.getProperties()
        tolerances = optimizer.getTolerances()
        for prop in properties:
            if protocol.name in protocolsHash[prop]:
                out_dict[prop] = {}
                out_dict[prop]['tol'] = tolerances[prop]
                out_dict[prop]['err'] = self.get_property_err(prop)
        return out_dict

    def initProtocolLengths(self, protocols):
        for prot in protocols:
            # Only init if the value is not set, otherwise you may end
            # up overwriting and extended simulation that was "kept"
            # between grid shifts.
            if prot.name not in self.protocol_lenspecs.keys():
                self.protocol_lenspecs[prot.name] = prot.calc_initial_len()

    def getProtocolLength(self, protocol):
        return self.protocol_lenspecs[protocol.name]

    def getParameterValues(self):
        return self.baseGrid.parSpaceGen.getParameterValues(self.id)

class ParameterGrid:

    _mainString = 'main'

    def __init__ (self, parSpaceGen, topologyBundles, reweighter, shifter, workdir):
        """
        Parameters:
          ParameterSpaceGenerator   parSpaceGen
          dict<TopologyBundle>      topologyBundles
          ReweighterInterface       reweighter
          GridShifter               shifter
          string                    workdir
        """
        self.parSpaceGen     = parSpaceGen
        self.topologyBundles = topologyBundles
        self.grid_points     = []
        self.fixed_points    = []
        self.indexGrid       = CartesianGrid(self.get_size())
        self.reweighter      = reweighter
        self.shifter         = shifter
        self.workdir         = workdir

        # Initialize GridPoints
        for i in range(self.get_linear_size()):
            self.grid_points.append( GridPoint(self, i) ) 

    def merge(self, other):
        self.shifter.merge(other.shifter)

    @staticmethod
    def load_from_binary(fn):
        fp = open(fn, 'rb')
        obj = pickle.load(fp)
        fp.close()
        return obj

    @staticmethod
    def createParameterGrid(parSpaceGen, topologyBundles, samples,
                            reweighterType, reweighterFactory,
                            shifterFactory, shifterArgs, workdir,
                            keep_initial_samples=True, validateFlag=False):
        from gmak.gridshifter import EmptyGridShifter
        parameterGrid        = ParameterGrid(parSpaceGen, topologyBundles,
                                             EmptyReweighter(),
                                             EmptyGridShifter(), workdir)
        parameterGrid.init = True
        parameterGrid.reweighter = reweighterFactory.create(reweighterType, parameterGrid)
        if (validateFlag):
            if shifterArgs is not None:
                shifterArgs['maxshifts'] = [0,] # must be a list to mimim reading from input
        parameterGrid.shifter    = shifterFactory(parameterGrid, shifterArgs)
        parameterGrid.set_samples(samples)
        if (keep_initial_samples):
            parameterGrid.set_fixed_points(samples)
        return parameterGrid

    @staticmethod
    def createParameterGridFromStream(stream, parSpaceGen, topologyBundles,
                                      reweighterFactory, shifterFactory,
                                      shifterArgs, workdir, validateFlag):
        samples       = []
        keep_initial_samples = True
        reweighterType = "standard"
        # Assumes last line read was '$grid'.
        for line in stream:
            if line[0] == '#':
                continue
            if (re.match(r"^\$end.*",line)):
                break
            if (line.split()[0] == 'samples'):
                samples = [int(x) for x in line.split()[1:]]
            if (line.split()[0] == 'reweight'):
                reweighterType = line.split()[1]
            if (line.split()[0] == 'fixsamples'):
                if (line.split()[1] == 'yes'):
                    keep_initial_samples = True
                elif (line.split()[1] == 'no'):
                    keep_initial_samples = False
                else:
                    raise ValueError("fixsamples can only be 'yes' or 'no'")
        grid = ParameterGrid.createParameterGrid(parSpaceGen, topologyBundles,
                                                 samples,
                                                 reweighterType,
                                                 reweighterFactory,
                                                 shifterFactory, shifterArgs,
                                                 workdir, keep_initial_samples,
                                                 validateFlag)
        return grid

    def initProtocolLengths(self, protocols):
        for gp in self.grid_points:
            gp.initProtocolLengths(protocols)

    def setGridpoints(self, gridpoints):
        self.grid_points = gridpoints

    def getCartesianGrid(self):
        return self.indexGrid

    def getParameterNames(self):
        return self.parSpaceGen.getParameterNames()

    def getParameterValues(self):
        return self.parSpaceGen.getAllParameterValues()

    def get_molecules(self):
        return list(self.topologyBundles.keys())

    def writeTopologies(self):
        grid = self.shifter.get_current_number_of_shifts()
        for state, gridpoint in enumerate(self.grid_points):
            if gridpoint.is_sample:
                params = self.parSpaceGen.getStateParameters(state)
                for systemName, topoBundle in self.topologyBundles.items():
                    workdir = os.path.join(self.workdir, topoBundle.name)
                    # ensure that workdir exists; if it does not, create it
                    if not os.path.isdir(workdir):
                        os.mkdir(workdir)
                    topoOut = topoBundle.write_topology(
                        workdir, grid, state, params)
                    # save topology output to gridpoint
                    self.grid_points[state].topology_outputs[systemName] = topoOut


    def writeParameters(self):
        self.parSpaceGen.writeParameters(self.makePrefixOfParameters())

    def get_dim(self):
        return self.parSpaceGen.getDimension(self._mainString)

    def get_size(self):
        return self.parSpaceGen.getSizes(self._mainString)

    def get_linear_size(self):
        return self.parSpaceGen.getNumberOfStates()

    def __getitem__(self, i):
        return self.grid_points[i]

    def get_square_size(self):
        allSizes  = self.get_size()
        firstSize = allSizes[0]
        for s in allSizes:
            if s != firstSize:
                raise ValueError("Can't get square size of non-square grid.")
        return int(firstSize)

    def get_number_of_samples (self):
        num = 0
        for x in self.grid_points:
            if x.is_sample:
                num += 1
        return num

    def get_samples (self):
        sp = []
        for x in self.grid_points:
            if x.is_sample:
                sp.append(x)
        return sp

    def get_samples_id (self):
        sp = []
        for x in self.grid_points:
            if x.is_sample:
                sp.append(x.id)
        return sp

    def save_samples_to_file (self, fn):
        sp = self.get_samples_id()
        np.savetxt(fn, sp, fmt='%d')
        return

    def read_samples_from_file(self, fn):
        sp = np.loadtxt(fn, dtype=int)
        self.set_samples(sp)
        return

    def set_samples (self, samples_list):
        for x in samples_list:
            self.grid_points[x].set_as_sample()

    def set_fixed_points(self, points_list):
        self.fixed_points = copy.deepcopy(points_list)

    def add_sample (self, new_sample):
        self.grid_points[new_sample].set_as_sample()

    def linear2tuple(self, linpos):
        """Returns tuple position for linear position."""
        return self.indexGrid.linear2tuple(linpos)

    def tuple2linear(self, pos):
        """Returns linear position for tuple position."""
        return self.indexGrid.tuple2linear(pos)

    def add_corners(self):
        """Add corners as sampling points. Useful when interpolation is used, to
        guarantee that all points of the grid are estimated."""
        for p in self.indexGrid.getCornersAsLinear():
            self.add_sample(p)

    def add_fixed_points(self):
        """Add fixed points as sampling points."""
        for p in self.fixed_points:
            self.add_sample(p)

    def simulate_with_protocol_at_dir (self, protocol, workdir):
        for i,gp in enumerate(self.grid_points):
            # run 
            if gp.is_sample:
                if not (gp.wasSimulatedWithProtocol(protocol)):
                    logger.globalLogger.putMessage("MESSAGE: GridPoint {} will be simulated or extended up to {} steps.".format(gp.id, gp.getProtocolLength(protocol)), dated=True)
                    logger.globalLogger.indent()
                    gp.simulate_with_protocol_at_dir(protocol, gp.makeSimudir(workdir))
                    gp.setProtocolAsSimulated(protocol)
                    logger.globalLogger.unindent()
                else:
                    logger.globalLogger.putMessage("MESSAGE: GridPoint {} has already been simulated with protocol {}!".format(gp.id, protocol.name))
            # only prepare
            else:
                logger.globalLogger.putMessage('MESSAGE: GridPoint {} is not a sample.'.format(gp.id))
                gp.prepare_with_protocol_at_dir(protocol, gp.makeSimudir(workdir, makeDir=False))


    def compute_and_filter_properties_with_protocol_at_dir(self, protocol, workdir):
        skips_dict = {}
        for i,gp in enumerate(self.grid_points):
            if gp.is_sample:
                filtering_properties = protocol.get_filtering_properties()
                for prop in protocol.get_properties():
                    oprop = os.path.join(gp.makeSimudir(workdir), f"{prop}.xvg")
                    gp.get_component_property_from_protocol(
                        prop,
                        protocol,
                        oprop)
                    if prop in filtering_properties:
                        ofilter = os.path.join(gp.makeSimudir(workdir),
                                               f"filtered_{prop}.xvg")
                        skips_dict[prop] = filter_datafile(oprop, ofilter)
                        # update path
                        gp.component_properties[protocol.name][prop] = ofilter
        return skips_dict


    def reweight(self, protocol, workdir):
        self.reweighter.run(protocol, workdir)

    def retrieveReweightProperty(self, prop):
        return self.reweighter.getPropertyMatrix(prop)

    def retrieveReweightProperties(self):
        return self.reweighter.getFullPropertyMatrix()

    def retrieveReweightNumberOfConfigurations(self):
        return self.reweighter.getConfigurationMatrix()

    def add_property_driver(self, prop_driver):
        try:
            self.prop_drivers[prop_driver.name] = prop_driver
        except AttributeError:
            self.prop_drivers = dict()
            self.prop_drivers[prop_driver.name] = prop_driver

    def compute_final_properties(self):
        for prop_driver in self.prop_drivers.values():
            for gridpoint in self.grid_points:
                prop_driver.compute(gridpoint)

    def setNewCenterForParameters(self, i):
        self.parSpaceGen.setNewCenter(i)

    def setNewOriginForParameters(self, i):
        self.parSpaceGen.setNewOrigin(i)

    # Performs shifting operations.
    def shift (self, optimizer):
        optimizer.reset()
        return self.shifter.shift(optimizer)

    def save_property_values_to_file (self, prop, optimizer):
        # make data
        data = [x.estimated_properties[prop].value for x in self.grid_points]
        data = np.array(data)
        filename = "{}/{}_EA_k.dat".format(self.makeStepPropertiesdir(optimizer), prop)
        np.savetxt(filename, data)

    def save_property_err_to_file (self, prop, optimizer):
        # make data
        data = [x.estimated_properties[prop].err for x in self.grid_points]
        data = np.array(data)
        filename = "{}/{}_dEA_k.dat".format(self.makeStepPropertiesdir(optimizer), prop)
        np.savetxt(filename, data)

    def save_property_diff_to_file (self, prop, ref, optimizer):
        # make data
        data = [x.estimated_properties[prop].value - ref for x in self.grid_points]
        data = np.array(data)
        filename = "{}/{}_diff.dat".format(self.makeStepPropertiesdir(optimizer), prop)
        np.savetxt(filename, data)

    def plot_property_to_file (self, prop, optimizer):
        if (self.get_dim() > 2):
            warnings.warn("Can only plot 1-D or 2-D grids.")
            return
        # make data
        data = [x.estimated_properties[prop].value for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.get_size())
        cbox_label = ""
        cbox_limits_colors = ()
        cbox_limits = ()
        title = self.grid_points[0].estimated_properties[prop].get_label()
        filename = "{}/{}_EA_k.pdf".format(self.makeStepPropertiesdir(optimizer), prop)
        # plot
        # from grid_ana ...
        if (self.get_dim() == 2):
            plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())
        elif (self.get_dim() == 1):
            plot_1d_to_file(filename, title, self.xlabel, data, self.get_samples_id())
        else:
            return

    def plot_property_err_to_file (self, prop, optimizer):
        if (self.get_dim() > 2):
            warnings.warn("Can only plot 1-D or 2-D grids.")
            return
        # make data
        data = [x.estimated_properties[prop].err for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.get_size())
        cbox_label = ""
        cbox_limits_colors = ('white', 'blue')
        cbox_limits = (0.0, np.max(data))
        title = self.grid_points[0].estimated_properties[prop].get_label_err()
        filename = "{}/{}_dEA_k.pdf".format(self.makeStepPropertiesdir(optimizer), prop)
        # plot
        # from grid_ana ...
        if (self.get_dim() == 2):
            plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())
        elif (self.get_dim() == 1):
            plot_1d_to_file(filename, title, self.xlabel, data, self.get_samples_id())
        else:
            return

    def plot_property_diff_to_file (self, prop, propname, ref, optimizer):
        if (self.get_dim() > 2):
            warnings.warn("Can only plot 1-D or 2-D grids.")
            return
        # make data
        data = [x.estimated_properties[prop].value - ref for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.get_size())
        cbox_label = ""
        # hard-coded preferences
        cbox_limits_colors = ('red', 'white', 'blue')
        # hard-coded preferences
        limits = {\
                'density': 20,\
                'dhvap': 4,\
                  'gamma': 10, \
                  'ced': 10 * 1e-4, \
                  'gced': 10}
        cbox_limits = (-limits[propname], limits[propname])
        title = "$\\Delta^\\mathrm{exp}$" + self.grid_points[0].estimated_properties[prop].get_label()
        filename = "{}/{}_diff.pdf".format(self.makeStepPropertiesdir(optimizer), prop)
        # plot
        # from grid_ana ...
        if (self.get_dim() == 2):
            plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())
        elif (self.get_dim() == 1):
            plot_1d_to_file(filename, title, self.xlabel, data, self.get_samples_id())
        else:
            return

    def clean_reweight_residues (self, protocol):
        for gi in self.grid_points:
            if gi.is_sample:
                for gj_id in gi.rw_outputs:
                    gi_id = self.get_samples().index(gi)
                    xtc = gi.rw_outputs[gj_id][protocol.name]['xtc']
                    trr = gi.rw_outputs[gj_id][protocol.name]['trr']
                    edr = gi.rw_outputs[gj_id][protocol.name]['edr']
                    tpr = gi.rw_outputs[gj_id][protocol.name]['tpr']
                    gro = gi.rw_outputs[gj_id][protocol.name]['gro']
                    for res_file in [xtc,trr,edr,tpr]:
                        if (os.path.isfile(res_file)):
                            os.remove(res_file)

    def makeDir(self, dirname):
        system("mkdir -p {}".format(os.path.abspath(dirname)))
        return os.path.abspath(dirname)

    def resetWorkdir(self, workdir):
        self.workdir = workdir

    def makeCurrentWorkdir(self):
        shifts = self.shifter.get_current_number_of_shifts()
        outdir = "{}/grid_{}".format(self.workdir, shifts)
        return self.makeDir(outdir)

    def makeProtocolWorkdir(self, protocol):
        return self.makeDir("{}/{}".format(self.makeCurrentWorkdir(), protocol.name))

    def makeProtocolSimudir(self, protocol):
        return self.makeDir("{}/simu".format(self.makeProtocolWorkdir(protocol), protocol.name))

    def makeProtocolReweightdirs(self, protocol):
        rwdir = self.makeDir("{}/rw".format(self.makeProtocolWorkdir(protocol), protocol.name))
        mbardir = self.makeDir("{}/mbar".format(self.makeProtocolWorkdir(protocol), protocol.name))
        return rwdir, mbardir

    def makeProtocolEstimatedir(self, protocol, method):
        return self.makeDir("{}/{}/estimated_properties".format(self.makeProtocolWorkdir(protocol), method))

    def makeStepPropertiesdir(self, optimizer):
        return self.makeDir("{}/step_{}".format(self.makeCurrentWorkdir(), optimizer.getCurrentIteration()))

    def makePathOfPropertyEstimates(self, protocol, method, prop):
        avg = "{}/{}_EA_k.dat".format(self.makeProtocolEstimatedir(protocol, method), prop)
        err = "{}/{}_dEA_k.dat".format(self.makeProtocolEstimatedir(protocol, method), prop)
        return avg, err

    def makePathOfSamples(self, optimizer):
        return "{}/samples_{}.dat".format(self.makeCurrentWorkdir(), optimizer.getCurrentIteration())

    def makePrefixOfParameters(self):
        return "{}/parameters".format(self.makeCurrentWorkdir())

    def makePathOfAbsoluteIndexes(self):
        return "{}/grid_indexes.dat".format(self.makeCurrentWorkdir())

    def save_to_binary(self, optimizer):
        # save grid
        fn = self.makeStepPropertiesdir(optimizer) + "/grid.bin"
        fp = open(fn, 'wb')
        pickle.dump(self, fp, pickle.HIGHEST_PROTOCOL)
        fp.close()
        # save optimizer
        fn = self.makeStepPropertiesdir(optimizer) + "/optimizer.bin"
        fp = open(fn, 'wb')
        pickle.dump(optimizer, fp, pickle.HIGHEST_PROTOCOL)
        fp.close()

    def make_grid_for_protocols (self, protocols, optimizer):
        # save samples to file
        self.save_samples_to_file(self.makePathOfSamples(optimizer))

        # simulate all
        for protocol in protocols:
            workdir  = self.makeProtocolWorkdir(protocol)
            simu_dir = self.makeProtocolSimudir(protocol)
            logger.globalLogger.putMessage('BEGIN PROTOCOL {}'.format(protocol.name), dated=True)
            logger.globalLogger.indent()
            # simulate sampling points, filter trajectories and properties
            # gridpoints have access to the correct topology paths
            self.simulate_with_protocol_at_dir (protocol, simu_dir)
            logger.globalLogger.unindent()
            logger.globalLogger.putMessage('END PROTOCOL {}'.format(protocol.name), dated=True)

        # calculate properties
        for protocol in protocols:
            simu_dir = self.makeProtocolSimudir(protocol)
            # filter properties needed
            skips_dict = self.compute_and_filter_properties_with_protocol_at_dir(
                protocol,
                simu_dir)

        # also get parameter values
        X_ki = self.getParameterValues()

        # reweight
        for protocol in protocols:
            if (protocol.requires_reweight()):
                # PLACEHOLDER FOR FILTERING THE TRAJECTORIES
                raise NotImplementedError(f"Reweight is not supported right now.")
                logger.globalLogger.putMessage('BEGIN REWEIGHT PROTOCOL {}'.format(protocol.name), dated=True)
                logger.globalLogger.indent()
                rw_dir, mbar_dir = self.makeProtocolReweightdirs(protocol)
                self.reweight(protocol, rw_dir)
                A_pkn = self.retrieveReweightProperties()
                u_kn  = self.retrieveReweightProperty('potential')
                if (protocol.has_pv()):
                    pv_kn = self.retrieveReweightProperty('pV')
                    u_kn  = (u_kn + pv_kn) / (0.0083144626 * protocol.get_temperature())
                else:
                    u_kn  = u_kn / (0.0083144626 * protocol.get_temperature())
                N_k   = self.retrieveReweightNumberOfConfigurations()
                # retrieve MBAR model (one model for all properties!)
                mbar_model = protocol.get_mbar_model()
                # estimate properties
                estimate_dir = self.makeProtocolEstimatedir(protocol, mbar_model.kind)
                mbar_model.computeExpectations(A_pkn, u_kn, N_k, X_ki)
                for p, rw_prop in enumerate(protocol.get_reweighting_properties()):
                    fn_avg, fn_err = self.makePathOfPropertyEstimates(protocol, mbar_model.kind, rw_prop)
                    mbar_model.writeExpectationsToFile(fn_avg, fn_err, p)
                mbar_model.writeLogToDirectory("%s/details" % mbar_dir)
                logger.globalLogger.unindent()
                logger.globalLogger.putMessage('END REWEIGHT PROTOCOL {}'.format(protocol.name), dated=True)

        # non-reweighted properties
        for protocol in protocols:
            interp_models_props = protocol.get_interp_models_props()
            A_psn = []
            for p, (model, prop) in enumerate(interp_models_props):
                A_psn.append([])
                for s, gs in enumerate(self.get_samples()):
                    prop_file = gs.retrieve_component_property_from_protocol (prop, protocol)
                    A_psn[p].append([])
                    A_psn[p][s] = np.loadtxt(prop_file, comments=['#','@'],
                                             usecols=(-1,))
            I_s = self.get_samples_id()
            # estimate properties
            for p, (model, prop) in enumerate(interp_models_props):
                model.computeExpectations([A_psn[p]], I_s, tuple(self.get_size()), X_ki)
                fn_avg, fn_err = self.makePathOfPropertyEstimates(protocol, model.kind, prop)
                model.writeExpectationsToFile(fn_avg, fn_err, 0)

    def run(self, protocols, optimizer, surrogateModelHash, properties,
            protocolsHash, resultsAssembler, plotFlag=False):

        for protocol in protocols:
            if (protocol.requires_corners()):
                self.add_corners()
                break

        self.add_fixed_points()

        # initialize length of simulations
        # this applies only to the samples that have no saved lengths
        self.initProtocolLengths(protocols)

        # create topology files
        self.writeTopologies()

        # create parameters file
        self.writeParameters()

        # create absolute-index file
        self.writeAbsoluteIndexes()

        self.make_grid_for_protocols(protocols, optimizer)

        self.compute_final_properties()

        for prop in optimizer.getProperties():
            referenceValue = optimizer.referenceValues[prop]
            self.save_property_values_to_file (prop, optimizer)
            self.save_property_err_to_file (prop, optimizer)
            self.save_property_diff_to_file (prop, referenceValue, optimizer)
            #
            if (plotFlag):
                self.plot_property_to_file (prop, optimizer)
                self.plot_property_err_to_file (prop, optimizer)
                self.plot_property_diff_to_file (prop, properties[prop], referenceValue, optimizer)
        #
        optimizer.fillWithScores (self)
        optimizer.printToFile (self, self.makeStepPropertiesdir(optimizer) + "/optimizer_data.dat")
        optimizer.printToFile (self, self.makeStepPropertiesdir(optimizer) + "/full_data.dat", sorted=False)
        if (plotFlag):
            optimizer.plotToPdf (self, self.makeStepPropertiesdir(optimizer) + "/optimizer_score.pdf")

        # convert protocolsHash (values = list of protocol names) into
        # protocolsHashByObject (values = list of references to protocols)
        protocolsHashByObject = {}
        for prop in protocolsHash.keys():
            protocolsHashByObject[prop] = []
            for name in protocolsHash[prop]:
                for prot in protocols:
                    if (prot.name == name):
                        protocolsHashByObject[prop].append(prot)

        # set up protocol extensions --- mark gridpoints as unsampled
        # if necessary
        if self.setExtendedProtocolLengths(protocols,
                                           optimizer,
                                           protocolsHash):
            # no new samples, but also don't proceed to shifting
            nextSample = []
        else:
            # if all protocols are converged
            if len(self.fixed_points) == 0:
                # if don't fix samples
                nextSample = optimizer.determineNextSample(self,
                                                           surrogateModelHash,
                                                           protocolsHashByObject)
            else:
                # proceed to shifting
                nextSample = -1

        if (nextSample == -1):
            # update current results
            globalState.update_record_book()
            # update results assembler
            for gs in self.get_samples_id():
                for prop in properties:
                    est = self[gs].get_property_estimate(prop)
                    err = self[gs].get_property_err(prop)
                    pars = self.parSpaceGen.getParameterValues(gs)
                    resultsAssembler.addData(pars, prop, est, err)

            if not self.shift(optimizer):
                globalState.saveToFile() # Save state to file if
                                         # needed for further
                                         # analysis.
                return
        else:
            for sample in nextSample:
                self.add_sample(sample)
        # Recursion
        globalState.saveToFile()
        self.run(protocols,
                 optimizer,
                 surrogateModelHash,
                 properties,
                 protocolsHash,
                 resultsAssembler,
                 plotFlag)

    def setExtendedProtocolLengths(self, protocols, optimizer, protocolsHash):
        """
        Set extended length for protocols.

        Returns True if some protocol, at some gridpoint, needs to be
        extended, and False if everything has converged. The new
        lengths are set if necessary.
        """
        output = False
        for gp in self.grid_points:
            if gp.is_sample:
                for prot in protocols:
                    old_length = gp.protocol_lenspecs[prot.name]
                    new_length = prot.calc_extend(gp, optimizer, protocolsHash)
                    if new_length is not None:
                        logger.globalLogger.putMessage('MESSAGE: GridPoint {} @ Protocol {} :'
                                                ' Steps : {}->{}'.
                                                format(gp.id, prot.name,
                                                       old_length, new_length))
                        gp.protocol_lenspecs[prot.name] = new_length
                        gp.unsetProtocolAsSimulated(prot)
                        output = True
                    else:
                        logger.globalLogger.putMessage('MESSAGE: GridPoint {} @ Protocol {} :'
                                                ' Steps have reached machine precision '
                                                'or there is no need to extend the '
                                                'simulations.'
                                                .format(gp.id, prot.name))
        if output:
            logger.globalLogger.putMessage('MESSAGE: Estimates are not converged,'
                                    ' so some simulations will be extended.')
        else:
            logger.globalLogger.putMessage('MESSAGE: Estimates are converged and '
                                    'simulations will not be extended.')
        return output

    def writeAbsoluteIndexes(self):
        fn = self.makePathOfAbsoluteIndexes()
        X = self.getAbsoluteIndexes()
        fp = open(fn, "w")
        # iterate over tuples
        for idx in X:
            # write to file
            for xi in idx:
                fp.write("%8d" % xi)
            fp.write("\n")
        fp.close()

    def getAbsoluteIndexes(self):
        idxs = []
        trans = self.shifter.index_transform
        # iterate over tuples
        for x in CartesianGridIterator(self.indexGrid):
            # transform to absolute index
            idx = trans.transform(x)
            # append
            idxs.append(idx)
        return idxs
