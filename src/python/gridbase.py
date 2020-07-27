import  re
import  os
import  shlex
import  numpy           as     np
from    reweightbase    import *
from    traj_ana        import *
from    traj_filter     import *
from    reweight        import *
import  multiprocessing
from    property        import * 
from    grid_ana        import *
import  warnings
from    surrogate_model import *
from    cartesiangrid   import *
from    gridshifter     import *
import  copy

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
        print ("Getting potential from %s and storing in %s." % \
        (edr, out_file))
        # get
        obtain_property (xtc, edr, gro, tpr, 'potential', out_file)
    if (protocol.has_pv()):
        out_file = workdir_mbar + \
            "/reweighted_properties/pV_%d_%d.xvg" % \
            (gi.id, gj_id)
        # If file already exists, no need to create it again, just set the variable.
        if not (os.path.isfile(out_file)):
            # copy from original trajectory
            original_file = gi.retrieve_atomic_property_from_protocol('pV', protocol)            
            print("Copying pV: {} -> {}".format(original_file, out_file))
            os.system("cp {} {}".format(original_file, out_file))
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

        atomic_properties: 2-d dictionary containing output files for
            the properties extracted,
            e.g. atomic_properties['liquid']['potential'] is the file
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
        self.atomic_properties = {}
        # self.rw_outputs = {}
        # also initialize subdictionaries - is this needed? <<todo-1>>
        # self.rw_outputs[self.id] = {}
        self.estimated_properties = {}

    def resetWithNewId(self, newId):
        self.id = newId
        self.estimated_properties = {}

    def wasSimulatedWithProtocol(self, protocol):
        return (protocol.name in self.is_simulated)

    def setProtocolAsSimulated(self, protocol):
        if (protocol.name not in self.is_simulated):
            self.is_simulated.append(protocol.name)
    
    def getTopologyPath(self, molecule):
        return self.baseGrid.topologyBundles[molecule].getPathsForStatepath(self.id)

    def get_property_estimate(self, prop):
        return self.estimated_properties[prop].value

    def get_property_err(self, prop):
        return self.estimated_properties[prop].err

    def prepare_with_protocol_at_dir (self, protocol, workdir):
        protocol.prepare_gridpoint_at_dir (self, workdir)

    def simulate_with_protocol_at_dir (self, protocol, workdir):
        if (self.wasSimulatedWithProtocol(protocol)):
            print("MESSAGE: Gridpoint {} has already been simulated with protocol {}!".format(self.id, protocol.name))
        else:
            protocol.run_gridpoint_at_dir (self, workdir)
            self.setProtocolAsSimulated(protocol)

    def reweight_with_protocol_at_other (self, protocol, gp_other, workdir):
        # from reweight.py
        if (protocol.type == 'slab'):
            out_rw = reweight (self.protocol_outputs[protocol.name]['trr'],\
                    self.protocol_outputs[protocol.name]['gro'],\
                    gp_other.protocol_outputs[protocol.name]['top'],\
                    protocol.mdps[-1], workdir)
        else:
            out_rw = reweight (self.protocol_outputs[protocol.name]['xtc'],\
                    self.protocol_outputs[protocol.name]['gro'],\
                    gp_other.protocol_outputs[protocol.name]['top'],\
                    protocol.mdps[-1], workdir)
        self.add_reweight_output (gp_other, protocol, out_rw)

    def set_as_sample (self):
        self.is_sample = True

    def is_sample (self):
        return self.is_sample

    # obj is a PropertyBase::X object (X=Density,Gamma,dHvap)
    def add_property_estimate (self, prop_id, prop_name, obj):
        if (prop_name != obj.name):
            raise ValueError ("ERROR: expected {} but got {}.\n".format(prop_name, obj.name))
        self.estimated_properties[prop_id] = obj

    # This is to ensure that any method that receives a grid point has access
    # to the file names of the simulations performed.
    #
    # protocol_outputs is a dictionary for the output (which are also
    # dictionaries) of a simulation, where keys are the protocol names.
    #
    def add_protocol_output (self, protocol, output):
        self.protocol_outputs[protocol.name] = output

    def add_reweight_output (self, gp_other, protocol, output):
        try:
            self.rw_outputs[gp_other.id][protocol.name] = output
        # key error means something has not been initialized
        except KeyError:
            self.rw_outputs[gp_other.id] = {}
            self.rw_outputs[gp_other.id][protocol.name] = output

    def get_atomic_property_from_protocol (self, name, protocol, output):
        # just to make sure absolute paths are used
        output = os.path.abspath(output)
        # from traj_ana.py
        if (name == 'polcorr'):
            obtain_polcorr (self.protocol_outputs[protocol.name]['xtc'],\
                    self.protocol_outputs[protocol.name]['edr'],\
                    self.protocol_outputs[protocol.name]['gro'],\
                    self.protocol_outputs[protocol.name]['tpr'],\
                    protocol.dipole, protocol.polar, output)
        else:
            obtain_property (self.protocol_outputs[protocol.name]['xtc'],\
                    self.protocol_outputs[protocol.name]['edr'],\
                    self.protocol_outputs[protocol.name]['gro'],\
                    self.protocol_outputs[protocol.name]['tpr'],\
                    name, output)
        self.add_atomic_property_output (name, protocol, output)

    def retrieve_atomic_property_from_protocol(self, propname, protocol):
        return self.atomic_properties[protocol.name][propname]

    def get_number_of_configurations_for_protocol(self, protocol):
        if (self.is_sample):
            confs = [0] * len(self.atomic_properties[protocol.name])      
            for i, prop_keys in enumerate(self.atomic_properties[protocol.name]):
                confs[i] = np.loadtxt(self.atomic_properties[protocol.name][prop_keys], comments=['@','#'], usecols=(0,)).shape[0]
                if (i > 0):
                    if confs[i] != confs[i-1]:
                        raise ValueError("Incompatible number of configurations in gridpoint {} for protocol {}: {} vs {}.".format(self.id, protocol.name, confs[i], confs[i-1]))
            return confs[0]
        else:
            return 0

    def add_atomic_property_output (self, prop, protocol, output):
        try:
            self.atomic_properties[protocol.name][prop] = output
        except KeyError:
            # this means atomic_properties[protocol-name] has not been initialized
            self.atomic_properties[protocol.name] = {}
            self.atomic_properties[protocol.name][prop] = output

    def filter_xtc_in_protocol (self, protocol, properties, odir):
        for prop in properties:
            self.get_atomic_property_from_protocol (prop, protocol,\
                    odir + "/" + prop + ".xvg")
        # from traj_filter
        extract_uncorrelated_frames (self.protocol_outputs[protocol.name]['xtc'],\
                self.protocol_outputs[protocol.name]['tpr'], \
                [self.atomic_properties[protocol.name][x] for x in properties],\
                odir + '/filtered_trajectory.xtc',\
                [odir + '/filtered_' + prop + '.xvg' for prop in properties])
                                     
        # update gridpoint trajectory
        self.protocol_outputs[protocol.name]['xtc'] = \
                os.path.abspath(odir + '/filtered_trajectory.xtc')
        # update path of filtered properties
        for x in properties:
            self.atomic_properties[protocol.name][x] = odir + '/filtered_' + x + '.xvg'

    def filter_trr_in_protocol (self, protocol, properties, odir):
        for prop in properties:
            self.get_atomic_property_from_protocol (prop, protocol,\
                    odir + "/" + prop + ".xvg")
        # from traj_filter
        extract_uncorrelated_frames (self.protocol_outputs[protocol.name]['trr'],\
                self.protocol_outputs[protocol.name]['tpr'], \
                [self.atomic_properties[protocol.name][x] for x in properties],\
                odir + '/filtered_trajectory.trr',\
                [odir + '/filtered_' + prop + '.xvg' for prop in properties])        
        # update gridpoint trajectory
        self.protocol_outputs[protocol.name]['trr'] = \
                os.path.abspath(odir + '/filtered_trajectory.trr')
        # update path of filtered properties
        for x in properties:
            self.atomic_properties[protocol.name][x] = odir + '/filtered_' + x + '.xvg'

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
        self.indexGrid       = CartesianGrid(self.get_size())
        self.reweighter      = reweighter
        self.shifter         = shifter
        self.workdir         = workdir
        
        """
        plotting stuff that I do not wish to alter right now.
        yMHG ter jul 14 14:40:02 -03 2020
        """
        self.xlabel     = "$X$"
        self.ylabel     = "$Y$"

        # Initialize GridPoints
        for i in range(self.get_linear_size()):
            self.grid_points.append( GridPoint(self, i) ) 

    @staticmethod
    def createParameterGrid(parSpaceGen, topologyBundles, samples, xlabel, ylabel, reweighterType, reweighterFactory, shifterFactory, shifterArgs, workdir):
        from gridshifter import EmptyGridShifter
        parameterGrid        = ParameterGrid(parSpaceGen, topologyBundles, EmptyReweighter(), EmptyGridShifter(), workdir)
        parameterGrid.xlabel = xlabel
        parameterGrid.ylabel = ylabel
        parameterGrid.reweighter = reweighterFactory.create(reweighterType, parameterGrid)
        parameterGrid.shifter    = shifterFactory(parameterGrid, shifterArgs)
        parameterGrid.set_samples(samples)
        return parameterGrid
        
    @staticmethod
    def createParameterGridFromStream(stream, parSpaceGen, topologyBundles, reweighterFactory, shifterFactory, shifterArgs, workdir):
        samples       = []
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
            if (line.split()[0] == 'labels'):
                splitted = shlex.split(line)
                xlabel = splitted[1]
                ylabel = splitted[2]
        grid = ParameterGrid.createParameterGrid(parSpaceGen, topologyBundles, samples, xlabel, ylabel, reweighterType, reweighterFactory, shifterFactory, shifterArgs, workdir)
        return grid

    def setGridpoints(self, gridpoints):
        self.grid_points = gridpoints
        
    def getCartesianGrid(self):
        return self.indexGrid
        
    def get_molecules(self):
        return list(self.topologyBundles.keys())

    def incrementPrefixOfTopologies(self):
        for topoBundle in self.topologyBundles.values():
            topoBundle.incrementPrefix()

    def writeTopologies(self):
        for topoBundle in self.topologyBundles.values():
            for i in range(self.get_linear_size()):
                self.parSpaceGen.setState(i)
                topoBundle.writeFilesForStatepath(i)

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

    def simulate_with_protocol_at_dir (self, protocol, workdir):
        for i,gp in enumerate(self.grid_points):
            # run 
            if gp.is_sample:
                # simulate
                gp.simulate_with_protocol_at_dir (protocol, workdir + "/" + str(i) + "/")
                # filter
                if (protocol.type == 'slab'):
                    gp.filter_trr_in_protocol (protocol, protocol.get_filtering_properties(), workdir + "/" + str(i) + "/")
                else:
                    gp.filter_xtc_in_protocol (protocol, protocol.get_filtering_properties(), workdir + "/" + str(i) + "/")
            # only prepare
            else:
                gp.prepare_with_protocol_at_dir (protocol, workdir + "/" + str(i) + "/")

    def reweight(self, protocol, workdir):
        self.reweighter.run(protocol, workdir)

    def retrieveReweightProperty(self, prop):
        return self.reweighter.getPropertyMatrix(prop)
        
    def retrieveReweightProperties(self):
        return self.reweighter.getFullPropertyMatrix()
        
    def retrieveReweightNumberOfConfigurations(self):
        return self.reweighter.getConfigurationMatrix()

    # # standard reweight + potential/pV recalculation
    # def parallel_clean_reweight_with_protocol_at_dir (self, protocol, workdir, workdir_mbar):
    #     properties = protocol.get_properties()
    #     ncpus_per_process = 2 
    #     nprocesses = multiprocessing.cpu_count() / ncpus_per_process
    #     pool = multiprocessing.Pool(nprocesses)
    #     pool_arg_list = []
    #     for i,gp in enumerate(self.grid_points):
    #         if gp.is_sample:
    #             for j,gp_other in enumerate(self.grid_points):
    #                 pool_arg_list.append([i,gp,j,gp_other,protocol,workdir,workdir_mbar,properties])
    #     pool.map(core_parallel_clean_reweight_i_at_j, pool_arg_list)
    
    # # fast reweight + potential/pV recalculation
    # def fast_clean_reweight_with_protocol_at_dir (self, protocol, reweightHash, workdir, workdir_mbar):
    #     properties = protocol.get_properties()
    #     for i,gp in enumerate(self.grid_points):
    #         if (gp.is_sample):
    #             # Clean this list.
    #             preReweightOutputs = []
    #             listOfPotentialFiles = []
    #             # Before reweighting, rerun the trajectory progressively
    #             # deactivating the energy components associated to each parameter.
    #             numberOfRuns = 2 * int( reweightHash['npars']) + 1
    #             for j in range(numberOfRuns):
    #                 # Modify the topology of the original simulation to refer to
    #                 # the proper *.itp file. Only 'j' is needed, as the names of
    #                 # the files are hard-coded.
    #                 os.system("mkdir -p %s/decouple_%d_at_%d" % (workdir, i, j))
    #                 topologyName = "%s/decouple_%d_at_%d/topol.top" % (workdir, i, j)
    #                 reweightItpChanger (gp.protocol_outputs[protocol.name]['top'], j,\
    #                         topologyName)
    #                 # Run the reweighting.
    #                 if (protocol.type == 'slab'):
    #                     preReweightOutput = reweight (gp.protocol_outputs[protocol.name]['trr'],\
    #                             gp.protocol_outputs[protocol.name]['gro'],\
    #                             topologyName,\
    #                             protocol.mdps[-1],\
    #                             "%s/decouple_%d_at_%d" % (workdir, i, j))
    #                 else:
    #                     preReweightOutput = reweight (gp.protocol_outputs[protocol.name]['xtc'],\
    #                             gp.protocol_outputs[protocol.name]['gro'],\
    #                             topologyName,\
    #                             protocol.mdps[-1],\
    #                             "%s/decouple_%d_at_%d" % (workdir, i, j))
    #                 # Save.
    #                 xtc = preReweightOutput['xtc']
    #                 edr = preReweightOutput['edr']
    #                 gro = preReweightOutput['gro']
    #                 tpr = preReweightOutput['tpr']
    #                 outPotential =  "%s/decouple_%d_at_%d/potential.xvg" % (workdir, i, j)
    #                 preReweightOutputs.append(preReweightOutput)
    #                 listOfPotentialFiles.append(outPotential)
    #                 # Get the potential energy of the reweighted trajectory.
    #                 obtain_property (xtc, edr, gro, tpr, 'potential', outPotential)
    #             # Run the reweighing program to obtain the reweighted potential
    #             # energies. This will always be necessary for MBAR.
    #             filesArg = " ".join(listOfPotentialFiles)
    #             fileGrid = reweightHash['parameters']
    #             anaRwBinary = reweightHash['program']
    #             command  = "%s -f %s -p %s -r %d -o %s/reweighted_properties/potential_%d_%%d.xvg -time" % (anaRwBinary, filesArg, fileGrid, i, workdir_mbar, i)
    #             print ("ANA_RW: %s" % command )
    #             # Create containing directory.
    #             os.system("mkdir -p %s/reweighted_properties/" % workdir_mbar)
    #             os.system(command)
    #             # Reweight other properties.
    #             # Redefine variables.
    #             # Note that rw = 0 corresponds to reweighting in the original state.
    #             xtc = preReweightOutputs[0]['xtc']
    #             edr = preReweightOutputs[0]['edr']
    #             tpr = preReweightOutputs[0]['tpr']
    #             gro = preReweightOutputs[0]['gro']
    #             gi  = gp
    #             # Calculate pV if necessary -- it is the same for all reweighted states.
    #             # In practice, this means COPYING the (filtered) pV data of the original simulation.
    #             if (protocol.has_pv()):
    #                 for gj in self.grid_points:
    #                     gj_id = gj.id 
    #                     out_file = workdir_mbar + \
    #                         "/reweighted_properties/pV_%d_%d.xvg" % \
    #                         (gi.id, gj_id)
    #                     # If file already exists, no need to create it again, just set the variable.
    #                     if not (os.path.isfile(out_file)):
    #                         # copy from original trajectory
    #                         original_file = gi.retrieve_atomic_property_from_protocol('pV', protocol)            
    #                         print("Copying pV: {} -> {}".format(original_file, out_file))
    #                         os.system("cp {} {}".format(original_file, out_file))

    # def reweight_mechanical_properties_for_protocol (self, protocol, mbar_dir):
    #     """Reweight mechanical properties for_protocol and stores results in rw_dir."""
    #     for prop in protocol.get_reweighting_properties():
    #         # first, ignore potential and pV
    #         if (prop == 'potential') or (prop == 'pV'):
    #             continue
    #         elif (prop == 'density') or (prop == 'polcorr') or (prop == 'volume'):
    #             # same for all states
    #             for gp in self.get_samples():
    #                 id = gp.id
    #                 input_file = gp.retrieve_atomic_property_from_protocol(prop, protocol)
    #                 for gpo in self.grid_points:
    #                     id_j = gpo.id
    #                     out_file = "%s/reweighted_properties/%s_%d_%d.xvg" % (mbar_dir, prop, id, id_j)
    #                     print("Copying density: {} -> {}".format(input_file, out_file))
    #                     os.system("cp %s %s" % (input_file, out_file))
    #         elif (prop == 'gamma'):
    #             # something needs to be done
    #             print ("ERROR: Reweighting of property " + prop + " is not implemented.")
    #             exit()                
    #         else:
    #             print ("ERROR: Reweighting of property " + prop + " is not implemented.")
    #             exit()

    def compute_final_properties (self, prop_id, prop, propProtocols, protocols, kind):
        protocolObjs = []
        for desiredProt in propProtocols:
            for prot in protocols:
                if prot.name == desiredProt:
                    protocolObjs.append(prot)

        if (prop == 'density'):
            estimateValue, estimateErr = protocolObjs[0].get_avg_err_estimate_of_property(prop, kind)
            for gp in self.grid_points:
                estimateObj   = Density (estimateValue[gp.id], estimateErr[gp.id])
                gp.add_property_estimate (prop_id, prop, estimateObj)
        if (prop == 'gamma'):
            estimateValue, estimateErr = protocolObjs[0].get_avg_err_estimate_of_property(prop, kind)
            for gp in self.grid_points:
                estimateObj   = Gamma (estimateValue[gp.id], estimateErr[gp.id])                
                gp.add_property_estimate (prop_id, prop, estimateObj)
        if (prop == 'dhvap'):
            estimateValueLiq, estimateErrLiq = protocolObjs[0].get_avg_err_estimate_of_property('potential', kind)
            nmols            = protocolObjs[0].nmols
            #
            estimateValueGas, estimateErrGas = protocolObjs[1].get_avg_err_estimate_of_property('potential', kind)                
            #
            estimateValuePol, estimateErrPol = protocolObjs[1].get_avg_err_estimate_of_property('polcorr', kind)                
            corr             = protocolObjs[1].other_corrections
            temp             = protocolObjs[1].get_temperature()
            #
            for gp in self.grid_points:
                estimateObj   = dHvap (estimateValueLiq[gp.id], estimateValueGas[gp.id], estimateValuePol[gp.id],
                                       estimateErrLiq[gp.id], estimateErrGas[gp.id], estimateErrPol[gp.id], corr, nmols, temp)
                gp.add_property_estimate (prop_id, prop, estimateObj)
        if (prop == 'ced'):
            estimateValueLiq, estimateErrLiq = protocolObjs[0].get_avg_err_estimate_of_property('potential', kind)
            nmols            = protocolObjs[0].nmols
            #
            estimateValueGas, estimateErrGas = protocolObjs[1].get_avg_err_estimate_of_property('potential', kind)                
            #
            estimateValuePol, estimateErrPol = protocolObjs[1].get_avg_err_estimate_of_property('polcorr', kind)                
            corr             = protocolObjs[1].other_corrections
            temp             = protocolObjs[1].get_temperature()
            #
            estimateValueVol, estimateErrVol = protocolObjs[0].get_avg_err_estimate_of_property('volume', kind)
            #
            for gp in self.grid_points:
                estimateObj   = CohesiveEnergyDensity (estimateValueLiq[gp.id], estimateValueGas[gp.id],
                                                       estimateValueVol[gp.id], estimateValuePol[gp.id],
                                                       estimateErrLiq[gp.id], estimateErrGas[gp.id],
                                                       estimateErrVol[gp.id], estimateErrPol[gp.id], corr, nmols)
                gp.add_property_estimate (prop_id, prop, estimateObj)
        if (prop == 'gced'):
            estimateValueLiq, estimateErrLiq = protocolObjs[0].get_avg_err_estimate_of_property('potential', kind)
            nmols            = protocolObjs[0].nmols
            #
            estimateValueGas, estimateErrGas = protocolObjs[1].get_avg_err_estimate_of_property('potential', kind)                
            #
            estimateValuePol, estimateErrPol = protocolObjs[1].get_avg_err_estimate_of_property('polcorr', kind)                
            corr             = protocolObjs[1].other_corrections
            temp             = protocolObjs[1].get_temperature()
            #
            estimateValueVol, estimateErrVol = protocolObjs[0].get_avg_err_estimate_of_property('volume', kind)
            #
            for gp in self.grid_points:
                estimateObj   = GammaViaCohesiveEnergyDensity (estimateValueLiq[gp.id], estimateValueGas[gp.id],
                                                       estimateValueVol[gp.id], estimateValuePol[gp.id],
                                                       estimateErrLiq[gp.id], estimateErrGas[gp.id],
                                                       estimateErrVol[gp.id], estimateErrPol[gp.id], corr, nmols)
                gp.add_property_estimate (prop_id, prop, estimateObj)

    def setNewCenterForParameters(self, i):
        self.parSpaceGen.setNewCenter(i)

    # Performs shifting operations.
    def shift (self, optimizer):
        optimizer.reset()
        return self.shifter.shift(optimizer)

    def read_property_values_err_from_file (self, prop, propType, filename_value, filename_err):
        values = np.loadtxt(filename_value)
        err = np.loadtxt(filename_err)
        for i,x in enumerate(self.grid_points):
            x.estimated_properties[prop] = init_property_from_string(propType, values[i], err[i])

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

    def makeCurrentWorkdir(self):
        shifts = self.shifter.getCurrentNumberOfShifts()
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

    def make_grid_for_protocol (self, protocol, optimizer):
        workdir  = self.makeProtocolWorkdir(protocol)
        simu_dir = self.makeProtocolSimudir(protocol)

        # create topology files
        self.writeTopologies()

        # save samples to file
        self.save_samples_to_file(self.makePathOfSamples(optimizer))

        # simulate sampling points, filter trajectories and properties
        # gridpoints have access to the correct topology paths 
        self.simulate_with_protocol_at_dir (protocol, simu_dir)

        # if any surrogate model requires reweight
        if (protocol.requires_reweight()):
            rw_dir, mbar_dir = self.makeProtocolReweightdirs(protocol)
            # BEGIN_OF_OLD_CODE
            # # standard reweight -- using 'rerun'
            # if (reweightHash['type'] == 'standard'):
            #     self.parallel_clean_reweight_with_protocol_at_dir (protocol, rw_dir, mbar_dir)
            # # fast reweight -- linear-basis function approach
            # elif (reweightHash['type'] == 'fast'):
            #     self.fast_clean_reweight_with_protocol_at_dir (protocol, reweightHash, rw_dir, mbar_dir)
            # else:
            #     raise ValueError("Reweight type must be 'standard' or 'fast'.")
            # # reweight mechanical properties besides 'potential' and 'pV'
            # self.reweight_mechanical_properties_for_protocol(protocol, mbar_dir)
            # END_OF_OLD_CODE
            self.reweight(protocol, rw_dir)
            # BEGIN_OF_OLD_CODE
            # # mount matrix with reweighted properties for all configurations
            # A_pkn = []
            # for p, rw_prop in enumerate(protocol.get_reweighting_properties()):
            #     A_pkn.append([])
            #     for k in range(self.get_linear_size()):
            #         A_pkn[p].append([])
            #         for s in self.get_samples_id():
            #             prop_file = "%s/reweighted_properties/%s_%d_%d.xvg" % (mbar_dir, rw_prop, s, k)
            #             A_pkn[p][k] = np.append(A_pkn[p][k], np.loadtxt(prop_file, usecols=(1,)))
            # A_pkn = np.array(A_pkn)
            # # mount u and pV matrixes
            # u_kn = []
            # for k in range(self.get_linear_size()):
            #     u_kn.append([])
            #     for s in self.get_samples_id():
            #         prop_file = "%s/reweighted_properties/potential_%d_%d.xvg" % (mbar_dir, s, k)
            #         u_kn[k] = np.append(u_kn[k], np.loadtxt(prop_file, usecols=(1,)))
            # u_kn = np.array(u_kn)
            # # sum pV if protocol requires it 
            # if (protocol.has_pv()):
            #     pv_kn = []
            #     for k in range(self.get_linear_size()):
            #         pv_kn.append([])
            #         for s in self.get_samples_id():
            #             prop_file = "%s/reweighted_properties/potential_%d_%d.xvg" % (mbar_dir, s, k)
            #             pv_kn[k] = np.append(pv_kn[k], np.loadtxt(prop_file, usecols=(1,)))
            #     pv_kn = np.array(pv_kn)
            #     u_kn += pv_kn
            # # mount N_k
            # N_k = []
            # for k, gp in enumerate(self.grid_points):
            #     N_k.append(gp.get_number_of_configurations_for_protocol(protocol))
            # N_k = np.array(N_k)
            # END_OF_OLD_CODE
            A_pkn = self.retrieveReweightProperties()
            u_kn  = self.retrieveReweightProperty('potential')
            if (protocol.has_pv()):
                pv_kn = self.retrieveReweightProperty('pV')
                u_kn  = (u_kn + pv_kn) / (0.83144626 * protocol.get_temperature())
            else:
                u_kn  = u_kn / (0.83144626 * protocol.get_temperature())                
            N_k   = self.retrieveReweightNumberOfConfigurations()
            # retrieve MBAR model (one model for all properties!)
            mbar_model = protocol.get_mbar_model()
            # estimate properties
            estimate_dir = self.makeProtocolEstimatedir(protocol, mbar_model.kind)
            mbar_model.computeExpectations(A_pkn, u_kn, N_k)
            for p, rw_prop in enumerate(protocol.get_reweighting_properties()):
                fn_avg, fn_err = self.makePathOfPropertyEstimates(protocol, mbar_model.kind, rw_prop)
                mbar_model.writeExpectationsToFile(fn_avg, fn_err, p)
            mbar_model.writeLogToDirectory("%s/details" % mbar_dir)

        # non-reweighted properties
        interp_models_props = protocol.get_interp_models_props()
        A_psn = []
        for p, (model, prop) in enumerate(interp_models_props):
            A_psn.append([])
            for s, gs in enumerate(self.get_samples()):
                prop_file = gs.retrieve_atomic_property_from_protocol (prop, protocol)                
                A_psn[p].append([])
                A_psn[p][s] = np.loadtxt(prop_file, usecols=(1,))
        I_s = self.get_samples_id()
        # estimate properties
        for p, (model, prop) in enumerate(interp_models_props):
            for kind in [model.kind, 'nearest']:
                model.computeExpectations([A_psn[p]], I_s, tuple(self.get_size()))
                fn_avg, fn_err = self.makePathOfPropertyEstimates(protocol, kind, prop)
                model.writeExpectationsToFile(fn_avg, fn_err, 0) # note that it is always property 0!

    def run(self, protocols, optimizer, surrogateModelHash, properties, protocolsHash, plotFlag=False):

        for protocol in protocols:
            if (protocol.requires_corners()):
                self.add_corners()
                break

        for protocol in protocols:
            self.make_grid_for_protocol(protocol, optimizer)
        
        for prop in properties:
            referenceValue = optimizer.referenceValues[prop]
            kind = surrogateModelHash[prop]
            self.compute_final_properties(prop, properties[prop], protocolsHash[prop], protocols, kind)
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
        if (plotFlag):
            optimizer.plotToPdf (self, self.makeStepPropertiesdir(optimizer) + "/optimizer_score.pdf")
        nextSample = optimizer.determineNextSample (self, surrogateModelHash)
        print ("Next sample is %d"  % nextSample)
        if (nextSample == -1):
            if not self.shift(optimizer):
                return
        else:
            self.add_sample(nextSample)
        # Recursion
        self.run(protocols, optimizer, surrogateModelHash, properties, protocolsHash, plotFlag)
            
    # type-hinted header is commented because it is not supported in old Python versions
    #def create_refined_subgrid(self, factors_list: list, model_str: str, propid2type: dict):            
    def create_refined_subgrid(self, factors_list, model_str, propid2type):
        # first check everything is compatible
        if (len(factors_list) != self.dim):
            raise ValueError
        for a_i in factors_list:
            if not(isinstance(a_i, int)) or not(a_i >= 1):
                raise ValueError
        # initialize model
        model = init_surrogate_model_from_string(model_str, False)
        if (model.kind == 'mbar'):
            raise ValueError("Can't use mbar for generating subgrid.")
        # copy this grid and alter data
        subgrid = ParameterGrid()
        subgrid.xlabel = self.xlabel
        subgrid.ylabel = self.ylabel
        subgrid.dim = self.dim 
        subgrid.size = [int(a_i * d_i - a_i + 1) for a_i, d_i in zip(factors_list, self.size)]
        subgrid.linear_size = int(np.prod(subgrid.size))
        # the list of grid_points is a bit tricky to be created
        # what needs to be done is the following:
        #
        #     a - If idx does not correspond to an old sample, create a new empty grid point with id = idx.
        #     b - If idx corresponds to an old estimated point, copy the corresponding grid point from self.grid_points[idx]
        #
        # by definition of refinement, an idx corresponds to an old
        # estimated point iff each coordinate x_i of linear2tuple(idx) is
        # such that x_i % a_i == 0, where a_i \in factors_list.
        subgrid.grid_points = [None] * subgrid.linear_size
        I_e = []
        for gp in self.grid_points:
            old_coordinates = list(self.linear2tuple(gp.id))
            new_coordinates = [int(x_i * a_i) for x_i, a_i in zip(old_coordinates, factors_list)]
            new_linear_position = subgrid.tuple2linear(tuple(new_coordinates))
            subgrid.grid_points[new_linear_position] = copy.deepcopy(gp)
            # change id
            subgrid.grid_points[new_linear_position].id = new_linear_position
            I_e.append(new_linear_position)
        # those that do not come from old estimated points become new grid points
        for idx, gp in enumerate(subgrid.grid_points):
            if (gp is None):
                subgrid.grid_points[idx] = GridPoint('', idx)
        # re-build property matrix for estimation
        list_of_property_names = list(subgrid.grid_points[I_e[0]].estimated_properties.keys())
        A_pe = [None] * len(list_of_property_names)
        dA_pe = [None] * len(list_of_property_names)        
        for p, prop_id in enumerate(list_of_property_names):
            A_pe[p] = [None] * self.linear_size # self.linear_size corresponds to the number of previously estimated points
            dA_pe[p] = [None] * self.linear_size # self.linear_size corresponds to the number of previously estimated points            
            for e in range(self.linear_size):
                A_pe[p][e] = subgrid.grid_points[I_e[e]].get_property_estimate(prop_id)
                dA_pe[p][e] = subgrid.grid_points[I_e[e]].get_property_err(prop_id)
        A_pe = np.array(A_pe)
        dA_pe = np.array(dA_pe)
        (EA_pk, dEA_pk) = model.computeExpectationsFromAvgStd(A_pe, dA_pe, I_e, tuple(subgrid.size))
        for p, prop_id in enumerate(list_of_property_names):
            for gp in subgrid.grid_points:
                propertyObj = init_property_from_string(propid2type[prop_id], EA_pk[p, gp.id], dEA_pk[p, gp.id])
                gp.add_property_estimate(prop_id, propid2type[prop_id], propertyObj)

        # return new instance
        return subgrid
