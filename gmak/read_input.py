import gmak.runcmd as runcmd
from gmak.reweightbase import ReweighterFactory
from gmak.gridoptimizer import gridOptimizer
from gmak.gridshifter import *
from gmak.parameter_space import *
from gmak.variations import *
from gmak.gridbase import *
from gmak.topology import *
from gmak.systems import *
import gmak.property
import gmak.configurations
import gmak.config
import gmak.component_properties as component_properties
from gmak.protocols import create_protocols
import gmak.interaction_parameter
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


def read_coordinates(blockdict, systems, coordinates, protocols, grid):
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
        if (prot.name != protocol_name):
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
        if (prot.name != protocol_name):
            raise InputError(f"Followed protocol \"{protocol_name}\" was not found.")
    else:
        raise InputError(f"Unknown configuration type: \"{_type}\".")

    # custom attributes
    out[1].set_custom_attributes_from_blockdict(blockdict, systems,
                                                coordinates, protocols)
    return out


def read_compute(blockdict: dict,
                 protocols: dict,
                 coordinates: dict,
                 systems: dict,
                 output_protocolsHash: dict,
                 bool_legacy: bool,
                 validate_flag: bool):
    base_options = [
        'name',
        'type',
        'surrogate_model',
        'protocols',
    ]
    verify_required(base_options, blockdict, "compute")
    name = blockdict['name'][0]
    type = blockdict['type'][0]
    if validate_flag:
        sms = 'empty'
    else:
        sms  = blockdict['surrogate_model'][0]
    prots = blockdict['protocols']
    try:
        component_props = blockdict['components']
    except KeyError:
        component_props = None
    propDriver = gmak.property.PropertyDriverFactory.create(
        name,
        type,
        sms,
        prots,
        [protocols[p] if p != "none" else None for p in prots],
        list(systems.values()),
        component_props)
    propDriver.set_custom_attributes_from_blockdict(
        blockdict, systems, coordinates, protocols)
    output_protocolsHash[name] = prots
    aps = propDriver.create_component_properties()
    for ap, prot in zip(aps, prots):
        if ap is None:
            # this means skipping something, e.g. gas-phase potential
            # energy, polcorr.
            continue
        else:
            # fill protocol
            prot_obj = protocols[prot]
            prot_obj.add_property(ap.name)
            prot_obj.add_component_property(ap)
            prot_obj.add_surrogate_model(sms, ap.name, bool_legacy)
    return propDriver


class ParameterOptionIterator:

    def __init__(self, iterable):
        self._iter = iterable.__iter__()
        self.ret = self._iter.__next__()

    def __iter__(self):
        return self

    def __next__(self):
        kwargs = {}
        if self.ret is None:
            raise StopIteration
        ret = self.ret
        # the range is the number of possible options
        while(True):
            try:
                option = self._iter.__next__()
                self.ret = option
            except StopIteration:
                self.ret = None
                return ret, kwargs
            if option == 'exclude':
                try:
                    value = self._iter.__next__()
                    kwargs['exclude_pairs'] = value
                except StopIteration:
                    raise ValueError
            elif option == 'include':
                try:
                    value = self._iter.__next__()
                    kwargs['include_pairs'] = value
                except StopIteration:
                    raise ValueError
            else:
                return ret, kwargs


class ParameterIO:
    def __init__(self, domainSpaceFactory):
        self.using              = None
        self.names              = []
        self.types              = []
        self.parRefs            = []
        self.domainSpaceFactory = domainSpaceFactory
        self.parSpaceGen        = ParameterSpaceGenerator()

    def getParameterSpaceGenerator(self):
        # This is called after all '$variation' blocks are read.
        return self.parSpaceGen

    def _addName(self, options):
        self.names.append(options['name'][0])

    def _updateUsing(self, options):
        try:
            usingString = options['using'][0]
            self.using = self.parSpaceGen.getDomainSpace(usingString)
        except KeyError:
            pass

    def _addParameterReferences(self, options):
        _theseParameters = []
        for refstring, kwargs in ParameterOptionIterator(options['pars']):
            param = gmak.interaction_parameter.InteractionParameter.from_string(
                refstring, **kwargs)
            _theseParameters.append(param)
        # Add list of parameter references to list of list of parameter references.
        self.parRefs.append(_theseParameters)

    def _addTypeAndDelegateRest(self, options):
        self.types.append(options['type'][0])
        # Create domain space and add to ParameterSpaceGenerator.
        ds = self.domainSpaceFactory.readFromTypeAndOptions(self.types[-1], options, self.using)
        self.parSpaceGen.addMember(self.names[-1], ds, self.parRefs[-1])

    def read(self, options):
        self._addName(options)
        self._updateUsing(options)
        self._addParameterReferences(options)
        self._addTypeAndDelegateRest(options)


def initialize_from_input (input_file, bool_legacy, validateFlag=False):
    output_grid = None
    output_gridshifter = None
    output_protocols = []
    # dict form of output_protocols
    output_protocols_dict = {}
    output_properties = {}
    output_protocolsHash = {}
    output_workdir = ""
    output_optimizer = None
    output_surrogateModel = {}

    # Dictionary of ConfigurationFactory objects where each key is the
    # 'name' option, i.e. a system/molecule name chosen by the user.
    output_coordinates = {}

    fp = open(input_file, "r")

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
                                                  outputTopoBundles,
                                                  output_coordinates,
                                                  output_protocols,
                                                  output_grid)
            output_coordinates[name] = conf_factory

        # grid
        if (line.split()[0] == "workdir"):
            output_workdir = os.path.abspath(line.split()[1])
        if (line.rstrip() == "$variation"):
            blockDict = block2dict(fp, '$end', '#')
            if parameterIO is None:
                parameterIO = ParameterIO(DomainSpaceFactory)
            parameterIO.read(blockDict)
        # This is for custom topology.
        if (line.rstrip() == '$system'):
            blockDict = block2dict(fp, '$end', '#')
            name = blockDict['name'][0]
            stype = blockDict['type'][0]
            outputTopoBundles[name] = create_system(stype, name, blockDict, outputTopoBundles,
                                                    output_coordinates, output_protocols,
                                                    exclude=None)
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
            new_protocol = create_protocols(blockDict, output_protocols_dict,
                                            output_coordinates,
                                            outputTopoBundles, output_grid)
            output_protocols.append(new_protocol)
            output_protocols_dict[new_protocol.name] = new_protocol
        if (line.rstrip() == '$compute'):
            blockDict = block2dict(fp, '$end', '#')
            propDriver = read_compute(blockDict,
                                      output_protocols_dict,
                                      output_coordinates,
                                      outputTopoBundles,
                                      output_protocolsHash,
                                      bool_legacy,
                                      validateFlag)
            # Add property driver to grid.
            output_grid.add_property_driver(propDriver)
        if (line.rstrip() == '$optimize'):
            blockDict = block2dict(fp, '$end', '#')
            output_optimizer = gridOptimizer.from_dict(blockDict, validateFlag)
        if (line.rstrip() == '$parameters'):
            output_paramLoop.readFromStream (fp)
        if (line.rstrip() == '$gridshift'):
            output_gridshifter = block2dict(fp, '$end', '#')
    fp.close()

    return (output_workdir, output_grid, output_protocols,
            output_properties, output_protocolsHash, output_optimizer,
            output_surrogateModel)
