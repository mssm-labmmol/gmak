#from RunFromInput import ParameterGrid
from reweightbase import ReweighterFactory
from gridoptimizer import gridOptimizer
from gridshifter import *
from parameters import *
from variations import *
from gridbase import *
from topology import * 
from coverage import coverInterface
from protocols import LiquidProtocol, GasProtocol, SlabProtocol, SolvationProtocol, SolvationFreeEnergyFactory, DummyProtocol

import os
import sys 

def block2dict(stream, endString, comment):
    out_dict = {}
    for line in stream:
        if line[0] == comment:
            continue
        if line.rstrip() == endString:
            break
        splittedLine = line.split()
        out_dict[splittedLine[0]] = splittedLine[1:]
    return out_dict

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
        # This knowledge is not the responsibility of this class.
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
            [atom, parameter] = refstring.split('/')
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
        
def initialize_from_input (input_file, bool_legacy):
    output_grid = None
    output_gridshifter = None
    output_protocols = []
    output_properties = {}
    output_protocolsHash = {}
    output_workdir = ""
    output_optimizer = gridOptimizer ()
    output_gaCoverInterface = coverInterface ()
    output_surrogateModel = {}
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
        if line[0] == '#':
            continue
        # grid
        if len(line.split()) < 1:
            continue
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
        if (line.rstrip() == '$molecules'):
            blockDict     = block2dict(fp, '$end', '#')
            for i, name in enumerate(blockDict['names']):
                prefix                     = "{}/{}/{}_0".format(output_workdir, name, name)
                moleculesDict[name]        = {}
                moleculesDict[name]['itp'] = blockDict['itp'][i]
                #moleculesDict[name]['gro'] = blockDict['gro'][i]
                outputTopoBundles[name]    = TopologyBundleFactory.createBundle('gromacs', blockDict['itp'][i], prefix,
                                                outputNonbondedForcefield, outputBondedForcefield)
                os.system("mkdir -p {}/{}".format(output_workdir, name))
                                                                                
        if (line.rstrip() == "$grid"):
            # This also signifies that no more $variation blocks exist
            # beyond this point.
            outputParameterSpaceGen = parameterIO.getParameterSpaceGenerator()
            output_grid = ParameterGrid.createParameterGridFromStream (fp, outputParameterSpaceGen, outputTopoBundles, ReweighterFactory, GridShifter.createFromGridAndDict, output_gridshifter, output_workdir)
        if (line.rstrip() == "$protocol"):
            line = next(fp)
            if (line.split()[0] == 'type'):
                typeRead = line.split()[1]
                if (typeRead == 'liquid'):
                    new_protocol = LiquidProtocol("",0,"",[""],[],[])
                    new_protocol.read_from_stream (fp)
                    output_protocols.append(new_protocol)
                elif (typeRead == 'gas'):
                    new_protocol = GasProtocol("","",0.0,0.0,[],[])
                    new_protocol.read_from_stream (fp)
                    output_protocols.append(new_protocol)
                elif (typeRead == 'slab'):
                    new_protocol = SlabProtocol("",[],5.0,[],-1)
                    new_protocol.read_from_stream (fp)
                    output_protocols.append(new_protocol)
                elif (typeRead == 'dgsolv'):
                    argsdict = block2dict(fp, '$end', '#')
                    new_protocols = SolvationFreeEnergyFactory().createSolvationProtocolFromDict(argsdict)
                    base_protocol = DummyProtocol(argsdict['name'][0], new_protocols[0])
                    base_protocol.point_to(new_protocols)
                    for prot in new_protocols:
                        output_protocols.append(prot)
                    output_protocols.append(base_protocol)
                else:
                    print ("ERROR: Type \"%s\" is not supported.\n" % typeRead)
                    exit()
            else:
                print ("ERROR: First line after a $protocol flag MUST assign its type.\n")
                exit()
        if (line.rstrip() == '$compute'):
            for line in fp:
                if line[0] == '#':
                    continue
                if line.rstrip() == '$end':
                    break
                propId   = line.split()[0]
                surrModel = line.split()[1]
                propRead = line.split()[2]
                nameRead = line.split()[3]

                if (surrModel == 'linear') or (surrModel == 'cubic'):
                    propIds = [propId, propId + '_nearest']
                    surrModels = [surrModel, 'nearest']
                else:
                    # mbar
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
                        # find protocol with name given
                        output_protocolsHash[propId] = [nameRead]
                        base_protocol = filter(lambda x: x.name == nameRead, output_protocols)
                        for protocol in base_protocol:
                            protocol.add_surrogate_model(surrModel, 'dgsolv', bool_legacy)           
                            protocol.properties.append('dgsolv')
                        # also
                        protocols = filter(lambda x: x.name.startswith(nameRead + '-'), output_protocols)
                        for protocol in protocols:
                            protocol.surrogate_models = []
                    else:
                        print ("ERROR: Property \"%s\" is not supported.\n" % typeRead)
                        exit()
        if (line.rstrip() == '$optimize'):
            output_optimizer.readFromStream (fp)
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
        if (line.rstrip() == '$gacover'):
            output_gaCoverInterface.readFromStream (fp)
    fp.close()

    return (output_workdir, output_grid, output_protocols, output_properties, \
            output_protocolsHash, output_optimizer, output_gaCoverInterface, output_surrogateModel, output_subgrid)
