# import base force-field classes
from parameters import *
from variations import * 
from abc import ABC, abstractmethod, abstractproperty

class TopologyInfo:
    def __init__(self, atoms, interPairs, intraPairs, bonds, angles, dihedrals, impropers):
        self.atoms = atoms
        self.interPairs = interPairs
        self.intraPairs = intraPairs
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.impropers = impropers

    def getAtoms(self):
        return [at.asDict() for at in self.atoms]

    def getIntermolecularPairs(self):
        return [pair.asDict() for pair in self.interPairs]

    def getIntramolecularPairs(self):
        return [pair.asDict() for pair in self.intraPairs]

    def getBonds(self):
        return [x.asDict() for x in self.bonds]

    def getAngles(self):
        return [x.asDict() for x in self.angles]

    def getDihedrals(self):
        return [x.asDict() for x in self.dihedrals]

    def getImpropers(self):
        return [x.asDict() for x in self.impropers]

class EmptyTopologyInfo(TopologyInfo):
    def __init__(self):
        super().__init__([], [], [], [], [], [], [])

class Topology:
    def __init__(self, bondedForcefield, nonbondedForcefield, topoInfo):
        self.bff = bondedForcefield
        self.nbff = nonbondedForcefield
        self.topo = topoInfo

    def updateNonbondedForcefield(self, newNbff):
        self.nbff = newNbff

    def updateBondedForcefield(self, newBff):
        self.bff = newBff

    def getBondtypes(self):
        return self.bff.getBondtypes()

    def getAngletypes(self):
        return self.bff.getAngletypes()

    def getDihedraltypes(self):
        return self.bff.getDihedraltypes()

    def getImpropertypes(self):
        return self.bff.getImpropertypes()

    def getAtomtypes(self):
        return self.nbff.getAtomtypes()

    def getPairtypes(self):
        return self.nbff.getPairtypes()

    def getAtoms(self):
        return self.topo.getAtoms()

    def getIntermolecularPairs(self):
        return self.topo.getIntermolecularPairs()

    def getIntramolecularPairs(self):
        return self.topo.getIntramolecularPairs()

class EmptyTopology(Topology):
    def __init__(self):
        super().__init__(EmptyBondedForcefield(), EmptyNonbondedForcefield(), EmptyTopologyInfo())

class AbstractTopologyInput(ABC):

    @abstractmethod
    def getTopology(self, fn):
        pass

# Concretions.
class GromacsDummyTopologyInput(AbstractTopologyInput):
    """
    This concretion does nothing.
    """
    def __init__(self, itp_fn):
        pass
    
    def getTopology(self):
        return EmptyTopology()
        
class AbstractTopologyOutput(ABC):

    @abstractmethod
    def writeToFiles(self, topology):
        pass

class GromacsDummyTopologyOutput(AbstractTopologyOutput):

    def __init__(self, itp_input, itp_fn):
        self.fn = itp_fn
        self.itp_input = itp_input

    def _l2z(self, name):
        if name[0] == 'H':
            return 1
        if name[0] == 'O':
            return 8
        if name[0] == 'C':
            return 6
        if name[0] == 'N':
            return 7

    def _writeBonded(self, topology):
        self.fp.write("; placeholder for bondedtypes\n")
        self.fp.write("\n")

    def _writeNonbonded(self, topology):
        atomtypes = topology.getAtomtypes()
        pairtypes = topology.getPairtypes()
        # Write atomtypes.
        self.fp.write("[ atomtypes ]\n")
        for label, pars in atomtypes:
            parameters = dict(pars)
            self.fp.write("%-5s%4d%6.3f%6.3f%3s%18.7f%18.7f\n" % (label, self._l2z(label), 0.0, 0.0, "A",
                                                             parameters['c6'], parameters['c12']))
        self.fp.write('\n')

        # Write normal pairs.
        self.fp.write("[ nonbond_params ]\n")
        for label, pars in pairtypes:
            parameters = dict(pars)
            if (label[0] != label[1]):            
                self.fp.write("%-6s%-6s%6d%18.7f%18.7f\n" % (label[0], label[1], 1,
                            parameters['c6'], parameters['c12']))
        self.fp.write('\n')

        # Write special pairs.
        self.fp.write("[ pairtypes ]\n")
        for label, pars in pairtypes:
            parameters = dict(pars)
            self.fp.write("%-6s%-6s%6d%18.7f%18.7f\n" % (label[0], label[1], 1,
                            parameters['cs6'], parameters['cs12']))
        self.fp.write('\n')

    def _writeTopoInfo(self, topology):
        _fp = open(self.itp_input, 'r')
        for line in _fp:
            self.fp.write(line)
        _fp.close()
    
    def writeToFiles(self, topology):
        self.fp = open(self.fn, 'w')

        # Write defaults block.
        self.fp.write("[ defaults ]\n; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n  1             1               no              1.0     1.0\n")
        self.fp.write('\n')

        # Write bonded.
        self._writeBonded(topology)

        # Write nonbonded.
        self._writeNonbonded(topology)

        # Write topo info.
        self._writeTopoInfo(topology)
        
        self.fp.close()

# ----------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------

class TopologyTestGromacsNoForcefield:

    def __init__(self, input_itp, output_itp):
        self.input_itp = input_itp
        self.output_itp = output_itp

    def run(self):
        topoinput    = GromacsDummyTopologyInput(self.input_itp)
        outputwriter = GromacsDummyTopologyOutput(self.input_itp, self.output_itp)
        topo         = topoinput.getTopology()
        outputwriter.writeToFiles(topo)

class TopologyTestGromacs:

    def __init__(self, input_itp, output_itp, nonbonded):
        self.input_itp  = input_itp
        self.output_itp = output_itp
        self.nonbonded  = nonbonded

    def run(self):
        topoinput    = GromacsDummyTopologyInput(self.input_itp)
        outputwriter = GromacsDummyTopologyOutput(self.input_itp, self.output_itp)
        topo         = topoinput.getTopology()

        topo.updateNonbondedForcefield(self.nonbonded)
        outputwriter.writeToFiles(topo)

class TopologyTestGromacsWithVariations:
    def __init__(self, input_itp, output_itp_preffix, nonbonded, spacegenerator):
        self.input_itp          = input_itp
        self.output_itp_preffix = output_itp_preffix
        self.nonbonded          = nonbonded
        self.spacegen           = spacegenerator

    def run(self):
        topoinput    = GromacsDummyTopologyInput(self.input_itp)
        topo         = topoinput.getTopology()

        for i in range(self.spacegen.getNumberOfStates()):
            output_itp = "{}_{}.itp".format(self.output_itp_preffix, i+1)
            outputwriter = GromacsDummyTopologyOutput(self.input_itp, output_itp)
            self.spacegen.setState(i)
            topo.updateNonbondedForcefield(self.nonbonded)
            outputwriter.writeToFiles(topo)

class TopologyBundle:

    def __init__(self, topology, topologyOutput, topologyOutputSetter):
        self.topology = topology
        self.topologyOutput = topologyOutput
        self.topologyOutputSetter = topologyOutputSetter # TODO: Setter comment

    def getTopology(self):
        return self.topology

    def writeFilesForStatepath(self, state):
        self.topologyOutputSetter.setState(self.topologyOutput, state)  # TODO: Setter comment
        self.topologyOutput.writeFiles()

    def getPathsForStatepath(self, state):
        self.topologyOutputSetter.setState(self.topologyOutput, state) # TODO: Setter comment
        self.topologyOutput.getFiles()

class TopologyBundleFactory:

    @staticmethod
    def _createBundleGromacs(itpPath, itpOutputPrefix, nonbondedForcefield, bondedForcefield):
        """
        Input object is an itp path.
        Output object is an itp path prefix.
        """
        # Initialize objects.
        _inp = GromacsDummyTopologyInput(itpPath)
        _top = _inp.getTopology()
        _set = GromacsDummyTopologyOutputSetter(itpOutputPrefix) # TODO: Setter comment
        _out = GromacsDummyTopologyOutput(itpPath, '')
        # By default, set output to state zero.
        _set.setState(_out, 0) # TODO: Setter comment
        # Update forcefield elements.
        _top.updateBondedForcefield(bondedForcefield)
        _top.updateNonbondedForcefield(nonbondedForcefield)
        # Create and return bundle.
        return TopologyBundle(_top, _out, _set)

    @staticmethod
    def createBundle(ctrlString, inputObject, outputObject, nonbondedForcefield, bondedForcefield):
        funct_dict = {
            'gromacs': TopologyBundleFactory.createBundleGromacs,
        }
        return funct_dict[ctrlString](inputObject, outputObject, nonbondedForcefield, bondedForcefield)

if __name__ == '__main__':
    from io import StringIO
    
    nonbonded = NonbondedForcefieldFactory.createFromStream(StringIO("standard\nCH3 1.0 2.0 3.0 4.0\nCH2 0.1 0.2 0.3 0.4\nCH1 0 2.0 1.0 2.0\n$end"))
    parameter_reference = NonbondedParameterReference(nonbonded.atomtypes, 'CH1', 'cs6')
    parameter_domain    = DomainSpace([VariationCartesian(1, 5, [0.0], [0.1], [5])])
    space_gen = ParameterSpaceGenerator()
    space_gen.addMember( parameter_domain, [parameter_reference] )
    
    test_1 = TopologyTestGromacsNoForcefield("/home/yan/programs/gridmaker/debug/abstr/et.itp", "/home/yan/programs/gridmaker/debug/abstr/et-output.itp")
    test_2 = TopologyTestGromacs("/home/yan/programs/gridmaker/debug/abstr/et.itp", "/home/yan/programs/gridmaker/debug/abstr/et-output-forcefield.itp", nonbonded)
    test_3 = TopologyTestGromacsWithVariations("/home/yan/programs/gridmaker/debug/abstr/et.itp", "/home/yan/programs/gridmaker/debug/abstr/et-varied", nonbonded, space_gen)
    
    test_1.run()
    test_2.run()
    test_3.run()
