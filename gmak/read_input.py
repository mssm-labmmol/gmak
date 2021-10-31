#from RunFromInput import ParameterGrid
import gmak.runcmd as runcmd
from gmak.reweightbase import ReweighterFactory
from gmak.gridoptimizer import gridOptimizer
from gmak.gridshifter import *
from gmak.parameters import *
from gmak.variations import *
from gmak.gridbase import *
from gmak.topology import *
import gmak.configurations
import gmak.atomic_properties as atomic_properties
from gmak.protocols import create_protocols
import os

# --------------------------------------------------------------------
# Custom errors for input.
# --------------------------------------------------------------------

class InputError(ValueError):
    pass


class InputMissingRequiredOptionError(InputError):
    pass

# --------------------------------------------------------------------
# Blockdict functions
# --------------------------------------------------------------------

def convert(value):
    types = [int, float, str]
    for t in types:
        try:
            return t(value)
        except ValueError:
            pass


def dict2types(bd):
    """Converts a block dict to the appropriate types."""
    return {k: [convert(x) for x in v] for k,v in bd.items()}
    

def block2dict(stream, endString, comment):
    out_dict = {}
    for line in stream:
        if line[0] == comment:
            continue
        if line.rstrip() == endString:
            break
        splittedLine = line.split()
        out_dict[splittedLine[0]] = splittedLine[1:]
    return dict2types(out_dict)

# --------------------------------------------------------------------
# Block-reading functions
# --------------------------------------------------------------------

def verify_required(require, blockdict, blockname):
    for required in require:
        if required not in blockdict:
            raise InputMissingRequiredOptionError(
                f"Option \"{required}\" not found in block ${blockname}.")

def get_optional(blockdict, optional_args):
    n = len(optional_args)
    out = []
    for arg in optional_args:
        if arg in blockdict.keys():
            if len(blockdict[arg]) > 1:
                raise InputError("get_optional only works for 1 value.")
            out.append(blockdict[arg][0])
        else:
            out.append(None)
    return tuple(out)

def read_coordinates(blockdict, protocols, grid):
    """
    Creates the Configuration object for a system based on the
    contents of the input file given in the `blockdict` dictionary.

    The list of protocols and the `grid` object are needed for
    follow-type protocols.
    """
    # Required options for any configuration. The other options may be
    # required by specific configuration types or be optional.
    require = [
        "name",
        "type",
    ]
    verify_required(require, blockdict, "coordinates")
    _type = blockdict['type'][0]
    _name = blockdict['name'][0]
    
    if _type == "any":
        # Required options: coords.
        verify_required(["coords"], blockdict, "coordinates/any")
        fn = blockdict['coords'][0]
        conf = gmak.configurations.ConfigurationFactory.from_file(fn)
        out = _name, gmak.configurations.FullCoordsConfigurationFactory(conf)
    elif _type == "gmx_liquid":
        # Required options: coords, nmols, box.
        verify_required(["coords", "nmols", "box"], blockdict, "coordinates/gmx_liquid")
        fn = blockdict['coords'][0]
        nmols = blockdict['nmols'][0]
        box_list = blockdict['box']
        box = gmak.configurations.BoxFactory.from_string_list(box_list)
        conf = gmak.configurations.ConfigurationFactory.from_file(fn)
        out = _name, gmak.configurations.GmxLiquidConfigurationFactory(
            conf, nmols, box)
    elif _type == "gmx_gas":
        # Required options: coords.
        verify_required(["coords"], blockdict, "coordinates/gmx_gas")
        fn = blockdict['coords'][0]
        box = gmak.configurations.RectangularBox(0.0, 0.0, 0.0)
        conf = gmak.configurations.ConfigurationFactory.from_file(fn)
        out = _name, gmak.configurations.GmxGasConfigurationFactory(conf)
    elif _type == "gmx_slab":
        # Required options: coords, nmols, box.
        # Optional: axis(default=z), factor(default=5.0)
        verify_required(["coords", "nmols", "box"], blockdict, "coordinates/gmx_slab")
        axis, factor = get_optional(blockdict, ['axis', 'factor']) 
        fn = blockdict['coords'][0]
        conf = gmak.configurations.ConfigurationFactory.from_file(fn)
        nmols = blockdict['nmols'][0]
        box_list = blockdict['box']
        box = gmak.configurations.BoxFactory.from_string_list(box_list)
        out =  _name, gmak.configurations.GmxSlabConfigurationFactory(
            conf, nmols, box, axis=axis, factor=factor)
    elif _type == "gmx_solvation":
        # Required options: coords (2), nmols (2), box.
        verify_required(["coords", "nmols", "box"], blockdict, "coordinates/gmx_solvation")
        fn_solute  = blockdict['coords'][0]
        fn_solvent = blockdict['coords'][1]
        conf_solute = gmak.configurations.ConfigurationFactory.from_file(fn_solute)
        conf_solvent = gmak.configurations.ConfigurationFactory.from_file(fn_solvent)
        nmols_solute = blockdict['nmols'][0]
        nmols_solvent = blockdict['nmols'][1]
        box_list = blockdict['box']
        box = gmak.configurations.BoxFactory.from_string_list(box_list)
        out =  _name, gmak.configurations.GmxSolvationConfigurationFactory(
            conf_solute, nmols_solute,
            conf_solvent, nmols_solvent,
            box)
    elif _type == "gmx_slab_follow":
        # Required options: follow.
        verify_required(['follow'], blockdict, "coordinates/gmx_slab_follow")
        axis, factor = get_optional(blockdict, ['axis', 'factor']) 
        protocol_name = blockdict['follow'][0]
        # Find protocol with given name.
        for prot in protocols:
            if prot.name == protocol_name:
                # Early exit.
                out = _name, gmak.configurations.GmxSlabFollowExtendConfigurationFactory(prot, grid, axis, factor)
                break
        # Protocol not found.
        raise InputError(f"Followed protocol \"{protocol_name}\" was not found.")
    elif _type == "follow":
        # Required options: follow.
        verify_required(['follow'], blockdict, "coordinates/follow")
        protocol_name = blockdict['follow'][0]
        # Find protocol with given name.
        for prot in protocols:
            if prot.name == protocol_name:
                # Early exit.
                out = _name, gmak.configurations.FollowProtocolConfigurationFactory(prot, grid)
                break
        # Protocol not found.
        raise InputError(f"Followed protocol \"{protocol_name}\" was not found.")
    else:
        raise InputError(f"Unknown configuration type: \"{_type}\".")

    # custom attributes
    out[1].set_custom_attributes_from_blockdict(blockdict)
    return out


class ParameterIO:
    def __init__(self, stream, parRefFactory, domainSpaceFactory):
        self.stream             = stream
        self.using              = None
        self.names              = []
        self.types              = []
        self.parRefs            = []
        self.parRefFactory      = parRefFactory
        self.domainSpaceFactory = domainSpaceFactory
        self.parSpaceGen        = ParameterSpaceGenerator()

    def getParameterSpaceGenerator(self):
        # This is called after all '$variation' blocks are read.
        return self.parSpaceGen

    def _addName(self, options):
        self.names.append( options[0] )
        return True

    def _updateUsing(self, options):
        usingString = options[0]
        self.using = self.parSpaceGen.getDomainSpace(usingString)
        return True

    def _addParameterReferences(self, options):
        _theseParameters = []
        for refstring in options:
            splitted_refstring = refstring.split('/')
            if len(splitted_refstring) == 2:
                [atom, parameter] = splitted_refstring
            elif len(splitted_refstring) == 1:
                atom, = splitted_refstring
                parameter = None
            else:
                raise ValueError("Can't parse refstring %s" % refstring)
            newReference = self.parRefFactory.create(atom, parameter)
            _theseParameters.append(newReference)
        # Add list of parameter references to list of list of parameter references.
        self.parRefs.append(_theseParameters)
        return True

    def _addTypeAndDelegateRest(self, options):
        self.types.append( options[0] )
        # Create domain space and add to ParameterSpaceGenerator.
        ds = self.domainSpaceFactory.readFromTypeAndStream(self.types[-1], self.stream, self.using)
        self.parSpaceGen.addMember(self.names[-1], ds, self.parRefs[-1])
        return False

    def read(self):
        # Dictionary of behaviors based on identifiers.
        funct_dict   = {
            'name': self._addName,
            'pars': self._addParameterReferences,
            'type': self._addTypeAndDelegateRest,
            'using': self._updateUsing,
        }
        for line in self.stream:
            if line[0] == '#':
                continue
            # Assumes the given stream has just read a '$variations' line.
            splittedLine = line.split()
            identifier   = splittedLine[0]
            options      = splittedLine[1:]
            # Execute behavior - funct_dict returns False after '$end' is reached in behavior.
            if not funct_dict[identifier](options):
                break
        return

# class TestParameterIO:
#     def __init__(self):
#         from io import StringIO
#         dsfac = DomainSpaceFactory

#         # for forcefield
#         txt = "standard\nCH3 1.0 2.0 3.0 4.0\nCH2 0.1 0.2 0.3 0.4\nCH1 0 2.0 1.0 2.0\n$end"
#         strio = StringIO(txt)
#         nbff = NonbondedForcefieldFactory.createFromStream(strio)
#         parreffac = ParameterReferenceFactory(nbff, EmptyBondedForcefield(), [EmptyTopology()])

#         # for variations
#         txt = "$variation\n"
#         txt += "name        main\n"
#         txt += "pars        CH3/c6           CH3/c12\n"
#         txt += "type        cartesian\n"
#         txt += "start       9.3471353e-03    2.6266240e-05\n"
#         txt += "step        2.336783825e-05  6.566559999999999e-08\n"
#         txt += "size        33               33\n"
#         txt += "$end\n"
#         txt += "\n"
#         txt += "$variation\n"
#         txt += "name     ties\n"
#         txt += "pars     CH2/c6  CH2/c12\n"
#         txt += "using    main\n"
#         txt += "type     scale\n"
#         txt += "factors  0.1      20\n"
#         txt += "$end\n"
#         strio = StringIO(txt)
#         self.pario = ParameterIO(strio, parreffac, dsfac)

#     def run(self):
#         # set at first variation
#         self.pario.stream.readline()
#         # do the thing
#         self.pario.read()
#         # set at next variation
#         self.pario.stream.readline()
#         self.pario.stream.readline()
#         # do the thing
#         self.pario.read()
#         # end
#         psgen = self.pario.getParameterSpaceGenerator()
#         # debug
#         psgen.debugPrint()

def initialize_from_input (input_file, bool_legacy, validateFlag=False):
    output_grid = None
    output_gridshifter = None
    output_protocols = []
    output_properties = {}
    output_protocolsHash = {}
    output_workdir = ""
    output_optimizer = None
    output_surrogateModel = {}

    # Dictionary of ConfigurationFactory objects where each key is the
    # 'name' option, i.e. a system/molecule name chosen by the user.
    output_coordinates = {}

    output_subgrid = {'method': 'cubic', 'factors': [1, 1]}
    fp = open(input_file, "r")

    # NEW
    outputNonbondedForcefield = EmptyNonbondedForcefield()
    outputBondedForcefield    = EmptyBondedForcefield()
    parameterIO               = None
    moleculesDict             = {}
    outputTopoBundles         = {}
    # END_NEW

    for line in fp:
        # skip blank and comment lines
        if line[0] == '#':
            continue
        if len(line.split()) < 1:
            continue

        # Read block-by-block.

        # $coordinates
        if (line.split()[0] == "$coordinates"):
            blockdict = block2dict(fp, "$end", "#")
            name, conf_factory = read_coordinates(blockdict,
                                                  output_protocols,
                                                  output_grid)
            output_coordinates[name] = conf_factory

        # grid
        if (line.split()[0] == "workdir"):
            output_workdir = os.path.abspath(line.split()[1])
        if (line.rstrip() == "$atomtypes"):
            outputNonbondedForcefield = NonbondedForcefieldFactory.createFromStream(fp)
        if (line.rstrip() == "$variation"):
            if parameterIO is None:
                parameterIO = ParameterIO(fp, ParameterReferenceFactory(outputNonbondedForcefield, outputBondedForcefield, [EmptyTopology()]), DomainSpaceFactory)
                parameterIO.read()
            else:
                parameterIO.read()
        # This is specific to GROMACS stuff.
        if (line.rstrip() == '$molecules') or (line.rstrip() == '$systems'):
            blockDict = block2dict(fp, '$end', '#')
            keys = blockDict.keys()
            for i, name in enumerate(blockDict['names']):
                prefix = "{}/{}/{}_0".format(output_workdir, name, name)
                molecularBlockDict = {}
                for k in keys:
                    molecularBlockDict[k] = blockDict[k][i]
                outputTopoBundles[name] = TopologyBundleFactory.createBundle(
                    'gromacs',
                    molecularBlockDict,
                    prefix,
                    (outputNonbondedForcefield, outputBondedForcefield))
                runcmd.run("mkdir -p {}/{}".format(output_workdir, name))
        # This is for custom topology.
        if (line.rstrip() == '$system'):
            blockDict = block2dict(fp, '$end', '#')
            name = blockDict['name'][0]
            stype = blockDict['type'][0]
            prefix = "{}/{}/{}_0".format(output_workdir, name, name)
            parRefList = parameterIO.getParameterSpaceGenerator().getParameterReferences()
            outputTopoBundles[name] = TopologyBundleFactory.createBundle(
                stype,
                blockDict,
                prefix,
                parRefList)
        if (line.rstrip() == "$grid"):
            # This also signifies that no more $variation blocks exist
            # beyond this point.
            outputParameterSpaceGen = parameterIO.getParameterSpaceGenerator()
            output_grid = ParameterGrid.createParameterGridFromStream(
                fp,
                outputParameterSpaceGen,
                outputTopoBundles,
                ReweighterFactory,
                create_gridshifter,
                output_gridshifter,
                output_workdir,
                validateFlag)
        if (line.rstrip() == "$protocol"):
            blockDict = block2dict(fp, '$end', '#')
            new_protocols = create_protocols(blockDict, output_protocols)
            output_protocols += new_protocols
        if (line.rstrip() == '$compute'):
            for line in fp:
                if line[0] == '#':
                    continue
                if line.rstrip() == '$end':
                    break
                propId   = line.split()[0]
                if (validateFlag):
                    surrModel = 'empty'
                else:
                    surrModel = line.split()[1]
                propRead = line.split()[2]
                nameRead = line.split()[3]

                propIds = [propId]
                surrModels = [surrModel]

                for propId, surrModel in zip(propIds, surrModels):
                    output_surrogateModel[propId] = surrModel
                    output_properties[propId] = propRead
                    if (propRead == 'density'):
                        output_protocolsHash[propId] = [nameRead]
                        # find protocol with name given
                        protocols = filter (lambda x: x.name == nameRead, output_protocols)
                        # actually no list -- this is only one protocol
                        for protocol in protocols:
                            protocol.add_surrogate_model(surrModel, 'density', bool_legacy)
                            protocol.properties.append('density')
                            # always add potential for safety
                            #if 'potential' not in protocol.properties:
                            #    protocol.properties.append('potential')
                    elif (propRead == 'dhvap'):
                        nameLiq = line.split()[3]
                        nameGas = line.split()[4]
                        corr    = float(line.split()[5])

                        if (nameGas == 'none'):
                            output_protocolsHash[propId] = [nameLiq, nameGas, corr]
                        else:
                            output_protocolsHash[propId] = [nameLiq, nameGas]
                        # find protocol with name given - Liq
                        protocols = filter (lambda x: x.name == nameLiq, output_protocols)
                        for protocol in protocols:
                            protocol.add_surrogate_model(surrModel, 'potential', bool_legacy) 
                            if 'potential' not in protocol.properties:
                                protocol.properties.append('potential')

                        # find protocol with name given - Gas
                        protocols = filter (lambda x: x.name == nameGas, output_protocols)
                        for protocol in protocols:
                            protocol.add_surrogate_model(surrModel, 'potential', bool_legacy)
                            if 'potential' not in protocol.properties:
                                protocol.properties.append('potential')
                            if protocol.polar >= 0:
                                protocol.add_surrogate_model(surrModel, 'polcorr', bool_legacy)                         
                                protocol.properties.append('polcorr')
                            protocol.set_other_corrections(corr)
                    elif (propRead == 'ced') or (propRead == 'gced'):
                        nameLiq = line.split()[3]
                        nameGas = line.split()[4]
                        corr    = float(line.split()[5])
                        output_protocolsHash[propId] = [nameLiq, nameGas]
                        # find protocol with name given - Liq
                        protocols = filter (lambda x: x.name == nameLiq, output_protocols)
                        for protocol in protocols:
                            protocol.add_surrogate_model(surrModel, 'potential', bool_legacy)
                            protocol.add_surrogate_model(surrModel, 'volume', bool_legacy)                         
                            if 'potential' not in protocol.properties:
                                protocol.properties.append('potential')
                            if 'volume' not in protocol.properties:
                                protocol.properties.append('volume')
                        # find protocol with name given - Gas
                        protocols = filter (lambda x: x.name == nameGas, output_protocols)
                        for protocol in protocols:
                            protocol.add_surrogate_model(surrModel, 'potential', bool_legacy)
                            if 'potential' not in protocol.properties:
                                protocol.properties.append('potential')
                            if (protocol.polar >= 0):
                                protocol.add_surrogate_model(surrModel, 'polcorr', bool_legacy)                                                         
                                protocol.properties.append('polcorr')
                            protocol.set_other_corrections(corr)
                    elif (propRead == 'gamma'):
                        # find protocol with name given
                        output_protocolsHash[propId] = [nameRead]
                        protocols = filter (lambda x: x.name == nameRead, output_protocols)
                        for protocol in protocols:
                            protocol.add_surrogate_model(surrModel, 'gamma', bool_legacy)           
                            protocol.properties.append('gamma')
                            # potential for safety
                            #if 'potential' not in protocol.properties:
                            #    protocol.properties.append('potential')
                    elif (propRead == 'dgsolv'):
                        output_protocolsHash[propId] = [nameRead]
                        # find solvation protocol
                        protocols = filter(lambda x: x.name == nameRead, output_protocols)
                        # add surrogate model and atomic property
                        for protocol in protocols:
                            protocol.add_surrogate_model(surrModel, propRead, bool_legacy)
                            protocol.properties.append(propRead)
                    else:
                        # Perhaps it is a custom property, try to create it
                        # just to check.
                        try:
                            atomic_property = atomic_properties.create_atomic_property(propRead)
                        except atomic_properties.PropertyNotInitialized:
                            # If it's not, then it's nothing.
                            print("ERROR: Property \"%s\" is not supported.\n"
                                  % typeRead)
                            exit()
                        output_protocolsHash[propId] = [nameRead]
                        protocols = filter(lambda x: x.name == nameRead,
                                           output_protocols)
                        for protocol in protocols:
                            protocol.add_surrogate_model(surrModel,
                                                         propRead,
                                                         bool_legacy)
                            protocol.properties.append(propRead)
        if (line.rstrip() == '$optimize'):
            blockDict = block2dict(fp, '$end', '#')
            output_optimizer = gridOptimizer.from_dict(blockDict, validateFlag)
        if (line.rstrip() == '$parameters'):
            output_paramLoop.readFromStream (fp)
        if (line.rstrip() == '$subgrid'):
            for line in fp:
                if line[0] == '#':
                    continue
                if line.rstrip() == '$end':
                    output_subgrid['parspacegen'] = outputParameterSpaceGen
                    break
                option = line.split()[0]
                if (option == 'method'):
                    output_subgrid[option] = line.split()[1]
                if (option == 'factors'):
                    output_subgrid[option] = [int(x) for x in line.split()[1:]]
        if (line.rstrip() == '$gridshift'):
            output_gridshifter = block2dict(fp, '$end', '#')
    fp.close()

    return (output_workdir, output_grid, output_protocols, output_properties, \
            output_protocolsHash, output_optimizer, output_surrogateModel, output_subgrid)
