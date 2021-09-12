# import base force-field classes
import gmak.runcmd as runcmd
from gmak.parameters import *
from gmak.variations import *
import os 
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

    @abstractmethod
    def getFiles(self):
        pass

class GromacsDummyTopologyOutput(AbstractTopologyOutput):

    def __init__(self, itp_input, itp_fn):
        self.fn = itp_fn
        self.itp_input = itp_input

    def _l2z(self, name):
        # TODO More atoms???
        if name[0] == 'X':
            # this is for testing purposes only
            return 1
        if name[0] == 'H':
            return 1
        if name[0] == 'O':
            return 8
        if name[0] == 'C':
            return 6
        if name[0] == 'N':
            return 7
        if name == 'AR':
            return 18

    def _alterFile(self, newFile):
        self.fn = newFile

    def _writeBonded(self, fp, topology):
        fp.write("; placeholder for bondedtypes\n")
        fp.write("\n")

    def _writeNonbonded(self, fp, topology):
        atomtypes = topology.getAtomtypes()
        pairtypes = topology.getPairtypes()
        # Write atomtypes.
        fp.write("[ atomtypes ]\n")
        for label, pars in atomtypes:
            parameters = dict(pars)
            fp.write("%-5s%4d%6.3f%6.3f%3s%18.7e%18.7e\n" % (label, self._l2z(label), 0.0, 0.0, "A",
                                                             parameters['c6'], parameters['c12']))
        fp.write('\n')

        # Write normal pairs.
        fp.write("[ nonbond_params ]\n")
        for label, pars in pairtypes:
            parameters = dict(pars)
            if (label[0] != label[1]):
                fp.write("%-6s%-6s%6d%18.7e%18.7e\n" % (label[0], label[1], 1,
                            parameters['c6'], parameters['c12']))
        fp.write('\n')

        # Write special pairs.
        fp.write("[ pairtypes ]\n")
        for label, pars in pairtypes:
            parameters = dict(pars)
            fp.write("%-6s%-6s%6d%18.7e%18.7e\n" % (label[0], label[1], 1,
                            parameters['cs6'], parameters['cs12']))
        fp.write('\n')

    def _writeTopoInfo(self, fp, topology):
        _fp = open(self.itp_input, 'r')
        for line in _fp:
            fp.write(line)
        _fp.close()

    def getFiles(self):
        return self.fn

    def writeToFiles(self, topology):
        fp = open(self.fn, 'w')

        # Guard force field part.
        fp.write("#ifndef FORCEFIELD_INCLUDE\n")
        fp.write("#define FORCEFIELD_INCLUDE\n")
        # Write defaults block.
        fp.write("[ defaults ]\n; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n  1             1               no              1.0     1.0\n")
        fp.write('\n')
        # Write bonded.
        self._writeBonded(fp, topology)
        # Write nonbonded.
        self._writeNonbonded(fp, topology)
        fp.write("#endif\n")

        # Write topo info.
        self._writeTopoInfo(fp, topology)

        fp.close()

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

# ----------------------------------------------------------------------
# Bundle things
# ----------------------------------------------------------------------

# Setter
class AbstractTopologyOutputSetter(ABC):

    @abstractmethod
    def setState(self, abstractTopologyOutput, state):
        pass

    @abstractmethod
    def incrementPrefix(self):
        pass

class GromacsDummyTopologyOutputSetter(AbstractTopologyOutputSetter):

    def __init__(self, itpOutputPrefix, ext):
        self.prefix = itpOutputPrefix
        self.ext = ext

    def incrementPrefix(self):
        indexOfUnderscore = self.prefix.rfind('_')
        number            = int(self.prefix[indexOfUnderscore+1:])
        front             = self.prefix[:indexOfUnderscore]
        number           += 1
        self.prefix       = "{}_{}".format(front, number)

    def setState(self, abstractTopologyOutput, state):
        newFile = "{}_{}.{}".format(self.prefix, state, self.ext)
        newFile = os.path.abspath(newFile)
        abstractTopologyOutput._alterFile(newFile)

# Bundle
class TopologyBundle:

    def __init__(self, topology, topologyOutput, topologyOutputSetter):
        self.topology = topology
        self.topologyOutput = topologyOutput
        self.topologyOutputSetter = topologyOutputSetter 

    def getTopology(self):
        return self.topology

    def incrementPrefix(self):
        self.topologyOutputSetter.incrementPrefix()

    def writeFilesForStatepath(self, state):
        self.topologyOutputSetter.setState(self.topologyOutput, state)  
        self.topologyOutput.writeToFiles(self.topology)

    def getPathsForStatepath(self, state):
        self.topologyOutputSetter.setState(self.topologyOutput, state) 
        return self.topologyOutput.getFiles()

class TopologyBundleFactory:

    @staticmethod
    def _createBundleGromacs(itpPath, itpOutputPrefix,
                             nonbondedForcefield, bondedForcefield,
                             ext):
        """
        Input object is an itp path.
        Output object is an itp path prefix.
        """
        # Initialize objects.
        _inp = GromacsDummyTopologyInput(itpPath)
        _top = _inp.getTopology()
        _set = GromacsDummyTopologyOutputSetter(itpOutputPrefix, ext)
        _out = GromacsDummyTopologyOutput(itpPath, '')
        # By default, set output to state zero.
        _set.setState(_out, 0)
        # Update forcefield elements.
        _top.updateBondedForcefield(bondedForcefield)
        _top.updateNonbondedForcefield(nonbondedForcefield)
        # Create and return bundle.
        return TopologyBundle(_top, _out, _set)

    @staticmethod
    def createBundle(ctrlString, inputObject, outputObject,
                     nonbondedForcefield, bondedForcefield,
                     fullTop):
        funct_dict = {
            'gromacs': TopologyBundleFactory._createBundleGromacs,
        }
        return funct_dict[ctrlString](inputObject, outputObject,
                                      nonbondedForcefield, bondedForcefield,
                                      fullTop)

class TestTopologyBundleGromacs:
    def __init__(self, itpFile, itpPrefix, nonbondedForcefield, bondedForcefield, spacegen):
        self.bundle = TopologyBundleFactory.createBundle('gromacs', itpFile, itpPrefix, nonbondedForcefield, bondedForcefield)
        self.spaceGen = spacegen
        self.nstates = spacegen.getNumberOfStates()

    def run(self):
        for i in range(self.nstates):
            self.spaceGen.setState(i)
            self.bundle.writeFilesForStatepath(i)

