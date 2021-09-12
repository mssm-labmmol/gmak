import  numpy       as      np
import gmak.runcmd as runcmd
from gmak.variations  import  DomainSpace
from    math        import  sqrt
from    copy        import  deepcopy
import  re

def geometric_mean(x, y):
    return np.sqrt(x*y)

# -----------------------------------------------------------------------------
# base classes for parameters
# -----------------------------------------------------------------------------

# Macro parameter
class MacroParameter:
    def __init__(self, token, value=None):
        self.token = str(token)
        if (value is not None):
            self.value = str(value)
        else:
            self.value = "0"

    # this is just for compatibility with other parts of the code
    def dereference(self):
        return self

    def alter(self, value):
        self.value = str(value)

    def get_full_string(self):
        return self.token

# Interface for NonbondedParameter
class NonbondedParameter:
    def value(self): 
        return self._value
    def label(self):
        return self._label
    def alter(self, value):
        self._value = value

# Special for mixing parameters.
#
# This abstraction ensures that these specific types of nb parameters
# are altered when the composing parameters are altered.
class NonbondedParameterMix:
    def __init__(self, label, nb_factory, par_i, par_j, comb_func):
        self.par_i = par_i
        self.par_j = par_j
        self.nbpar = nb_factory.create(comb_func(par_i.value(), par_j.value()), label)
        self.comb_func = comb_func

    def label(self):
        return self.nbpar.label()

    def value(self):
        self.nbpar.alter(self.comb_func( self.par_i.value(), self.par_j.value() ))
        return self.nbpar.value()
    
# Subtypes
class C6(NonbondedParameter):
    def __init__(self, value, label):
        self._value = value
        self._label = label
class CS6(NonbondedParameter):
    def __init__(self, value, label):
        self._value = value
        self._label = label
class C12(NonbondedParameter):
    def __init__(self, value, label):
        self._value = value
        self._label = label
class CS12(NonbondedParameter):
    def __init__(self, value, label):
        self._value = value
        self._label = label

class NonbondedParameterFactory:
    @staticmethod
    def create(value, typeString):
        if (typeString == 'c6'):
            return C6(value, typeString)
        if (typeString == 'cs6'):
            return CS6(value, typeString)
        if (typeString == 'c12'):
            return C12(value, typeString)
        if (typeString == 'cs12'):
            return CS12(value, typeString)

class NonbondedParameters:
    def __init__(self, factory, labels, values):
        self.parameters = {}
        for value, label in zip(values, labels): 
            self.parameters[label] = factory.create(value, label)
    def alter(self, labels, values):
        for value, label in zip(values, labels):
            self.parameters[label].alter(value)
    def getNonbondedParameterReference(self, label):
        return self.parameters[label]
    def getValue(self, label):
        return self.parameters[label].value()
    def debugPrint(self):
        for x in NonbondedParametersIterator(self):
            print(x)

class NonbondedParametersMix(NonbondedParameters):
    def __init__(self, factory, labels, pars_i, pars_j, comb_matrix, comb_funcs):
        self.parameters = {}
        labels = comb_matrix.keys()
        for l in labels:
            nbmix = NonbondedParameterMix(l, factory, pars_i.getNonbondedParameterReference(comb_matrix[l][0]),
                      pars_j.getNonbondedParameterReference(comb_matrix[l][1]), comb_funcs[l])
            
            self.parameters[l] = nbmix
        

class NonbondedParametersIterator:
    def __init__(self, nb):
        self.curr = -1
        self.keys = list(nb.parameters.keys())
        self.nb   = nb
    def __iter__(self):
        return self
    def __next__(self):
        self.curr += 1
        if self.curr >= len(self.nb.parameters):
            raise StopIteration
        return (self.nb.parameters[self.keys[self.curr]].label(), self.nb.parameters[self.keys[self.curr]].value())

# -----------------------------------------------------------------------------
# Atomtype
# -----------------------------------------------------------------------------

class Atomtype:

    def __init__(self, typeString, nonbondedParameters):
        self.type = typeString
        self.nonbondedParameters = nonbondedParameters

    def getLabel(self):
        return self.type

    def getNonbondedParameters(self):
        return [x for x in NonbondedParametersIterator(self.nonbondedParameters)]

    def getNonbondedParameterValue(self, label):
        return self.nonbondedParameters.getValue(label)

    def alterNonbondedParameters(self, labels, values):
        self.nonbondedParameters.alter(labels, values)

    def getNonbondedParameterReference(self, label):
        return self.nonbondedParameters.getNonbondedParameterReference(label)

    def getNonbondedParametersObject(self):
        return self.nonbondedParameters

    def decoupleNonbondedParameter(self, label):
        self.alter([label], [0.0])

    def __eq__(self, other):
        return (self.type == other.type)

    def debugPrint(self):
        print("atom {}".format(self.type))
        self.nonbondedParameters.debugPrint()

# -----------------------------------------------------------------------------
# ParameterReference
# -----------------------------------------------------------------------------

class NonbondedParameterReference(NonbondedParameter):
    """A reference to a particular parameter of an atomtype.
    Useful to directly alter this parameter."""
    def __init__(self, atomtypes, typeString, parString):
        for at in atomtypes:
            if (at.type == typeString):
                self.reference = at.getNonbondedParameterReference(parString)
                self.atomtype  = at
                return
        raise ValueError("Trying to reference an Atomtype/NonbondedParameter which does not exist.")
    
    def dereference(self):
        return self.reference

    def get_atomtype(self):
        try:
            return self.atomtype.type
        except AttributeError: # In old runs, self.atomtype is not
                               # defined while reading the input file,
                               # so an AttributeError is raised at
                               # this point when restarting from an
                               # old .bin file.
            return "XX"

    def get_full_string(self):
        return "{}/{}".format(self.get_atomtype(), self.reference.label())
        

class ParameterReferenceFactory:

    def __init__(self, nonbondedForcefield, bondedForcefield, topologies):
        self.nbff             = nonbondedForcefield
        self.bff              = bondedForcefield
        self.topos            = topologies
        self.bondedStrings    = ['q', 'm']
        
    def create(self, typestring, parameter):
        """
        str typestring : string identifying atomtype (e.g. 'OW') or bondedtype (e.g. 'CH2-OW' for bond)
        str parameter  : string identifying parameter (e.g. 'c6', 'q')
        """
        if ('-' in typestring):
            raise NotImplementedError("Creation of bonded types is not"
                                      " implemented.")
        elif (typestring[0] == '@'):
            return MacroParameter(typestring)
        else:
            if parameter in self.bondedStrings:
                raise NotImplementedError("Creation of bonded types is not implemented.")
            else:
                return NonbondedParameterReference(self.nbff.getAtomtypeObjects(),
                                                   typestring, parameter)
    
# -----------------------------------------------------------------------------
# ParameterSpaceGenerator
# -----------------------------------------------------------------------------

class ParameterSpaceGenerator:
    """
    Contains tuples of (DomainSpace,
    [NonbondedParameterReference/MacroParameter]) indexed in a
    dictionary.  Handles connecting data from DomainSpace with
    parameter values.

    """
    def __init__(self):
        self.members = {}

    def addMember(self, name, domainSpace, nbParRefList):
        """Add a (domainSpace, list<NonbondedParameterReference>) tuple."""
        # check compatibility of number of states
        for ds, pl in self.members.values():
            if (ds.get_linear_size() != domainSpace.get_linear_size()):
                raise ValueError("Trying to add a DomainSpace with non-matching linear size.")
        # check compatibility of dimension
        if (domainSpace.get_dim() != len(nbParRefList)):
            raise ValueError("Trying to associate DomainSpace of dimension {} to {} parameters.".format(
                domainSpace.get_dim(), len(nbParRefList)))
        self.members[name] = (domainSpace, nbParRefList)

    def setState(self, i):
        """Applies DomainSpace values at state 'i' to the parameters referenced by the object."""
        for domainSpace, parList in self.members.values():
            values = domainSpace.get(i)
            for j, parRef in enumerate(parList):
                par = parRef.dereference()
                par.alter( values[j] )

    def writeMainVariationBlock(self, stream):
        name = 'main'
        stream.write("$variation\n")
        # Writing name and parameter strings requires ParameterSpaceGenerator
        stream.write("name {}\n".format(name))
        atomtype_strings = [x.get_full_string() for x in self.members[name][1]]
        stream.write("pars " + " ".join(atomtype_strings) + "\n")
        # Writing the rest is delegated to DomainSpace
        self.members[name][0].write_block_to_stream(stream)
        stream.write("$end\n\n")

    def getParameterNames(self):
        return [x.get_full_string() for x in self.members['main'][1]]

    def getParameterValues(self, state):
        return self.members['main'][0].get(state)

    def getAllParameterValues(self):
        return self.members['main'][0].get_data()

    def getMacros(self):
        macros = []
        for name in self.members.keys():
            for refs in self.members[name][1]:
                if isinstance(refs, MacroParameter):
                    macros.append(refs)
        return macros

    def writeParameters(self, prefix):
        for name in self.members.keys():
            domainSpace, parList = self.members[name]
            outfile = prefix + "_" + name + ".dat"
            domainSpace.write_to_file(outfile)

    def setNewCenter(self, i):
        for name in self.members:
            self.members[name][0].set_new_center(i)

    def getDimension(self, name):
        return self.members[name][0].get_dim()

    def getSizes(self, name):
        return self.members[name][0].get_sizes()

    def getNumberOfStates(self):
        return list(self.members.values())[0][0].get_linear_size()

    def getDomainSpace(self, name):
        return self.members[name][0]

    def debugPrint(self):
        from sys import stderr
        outputStream = stderr
        for name in self.members:
            print("; domain space {}".format(name), file=outputStream)
            self.members[name][0].write_to_stream(outputStream)

    def copy(self):
        newObject = ParameterSpaceGenerator()
        for name in self.members:
            newObject.addMember(name, deepcopy(self.members[name][0]), self.members[name][1])
        return newObject

    def rescale(self, factors):
        for name in self.members:
            self.members[name][0].rescale(factors)

# -----------------------------------------------------------------------------
# NonbondedForcefield
# -----------------------------------------------------------------------------

class BasePairtype:

    def __init__(self, atomtype_i, atomtype_j, nb_factory, comb_matrix, comb_funcs):
        """Parameters:
        atomtype_i: <Atomtype>
        atomtype_j: <Atomtype>
        nb_factory: NonbondedParameterFactory
        comb_matrix: dict of tuples, e.g. {'c6': ('c6', 'c6'), 'c12': ('c12', 'c12')}
        comb_funcs:  dict of functions, e.g. {'c6': lambda x,y : np.sqrt(x*y), 'c12': lambda x,y : np.sqrt(x*y)} 
        """

        self.type_i = atomtype_i
        self.type_j = atomtype_j

        labels  = comb_matrix.keys()
        pars_i  = atomtype_i.getNonbondedParametersObject()
        pars_j  = atomtype_j.getNonbondedParametersObject()

        self.nonbondedParameters = NonbondedParametersMix(nb_factory, labels, pars_i, pars_j, comb_matrix, comb_funcs)

    def getParameter(self, label):
        return self.nonbondedParameters.getValue(label)

    def getParameters(self):
        return [x for x in NonbondedParametersIterator(self.nonbondedParameters)]

    def __eq__(self, other):
        return ((self.type_i == other.type_i) and (self.type_j == other.type_j)) or ((self.type_i == other.type_j) and (self.type_j == other.type_i))

    def getLabel(self):
        return (self.type_i.getLabel(), self.type_j.getLabel())
    
    def debugPrint(self):
        print("pair {}-{}".format(self.type_i.type, self.type_j.type))
        print("interaction parameters")
        self.nonbondedParameters.debugPrint()

class StandardLJPairtype(BasePairtype):
    """
    Facilitates initialization of a BasePairType in case of standard LJ interaction
    with geometric-mean combination rule.
    """
    def __init__(self, atomtype_i, atomtype_j, nb_factory):
        comb_matrix = {'c6': ('c6', 'c6'), 'c12': ('c12', 'c12'), 'cs6': ('cs6', 'cs6'), 'cs12': ('cs12', 'cs12')}
        comb_funcs  = {'c6': geometric_mean, 'c12': geometric_mean, 'cs6': geometric_mean, 'cs12': geometric_mean}
        super().__init__(atomtype_i, atomtype_j, nb_factory, comb_matrix, comb_funcs)

class PairtypeFactory:
    @staticmethod
    def createFromAtomtypes(ffType, atomtypes):
        natoms = len(atomtypes)
        pairtypes = []
        if ffType == 'standard':
            # standard LJ force-field
            for i in range(natoms):
                for j in range(i, natoms):
                    thisPairtype = StandardLJPairtype(atomtypes[i], atomtypes[j], NonbondedParameterFactory)
                    pairtypes.append(thisPairtype)
        else:
            raise ValueError("Forcefield type {} is not supported.".format(ffType))
        
        return pairtypes

class NonbondedForcefield:
    """This is the main API through which the program accesses parameter values and
    atomtype names and writes it to disk. Note that parameter values are not copied
    into this class, so any alteration in the dependencies will reflect in the values
    you can access via this class."""
    
    def __init__(self, atomtypes, pairtypes):
        """Parameters:
        atomtypes:     list of <Atomtype> objects that compose the forcefield
        pairtypes:     list of <Pairtype> objects for 1-4 interactions that compose the forcefield
        """
        self.atomtypes = atomtypes
        self.pairtypes = pairtypes
        
    def getPairtypes(self):
        """Return example: 
        [ (('CH2','CH2'), [('c6', 1.2e-03), ('c12', 1.7e-06)]), 
          (('CH2','CH3'), [('c6', 1.6e-03), ('c12', 1.1e-06)]), 
          ...
        ]
        """        
        return [(x.getLabel(), x.getParameters()) for x in self.pairtypes]

    def getAtomtypes(self):
        """Return example: 
        [ ('CH2', [('c6', 1.2e-03), ('c12', 1.7e-06)]), 
          ('CH3', [('c6', 1.6e-03), ('c12', 1.1e-06)]), 
          ...
        ]
        """        
        return [(x.getLabel(), x.getNonbondedParameters()) for x in self.atomtypes]

    def getAtomtypeObjects(self):
        return self.atomtypes
        
    def debugPrint(self):
        for at in self.atomtypes:
            at.debugPrint()
        for pt in self.pairtypes:
            pt.debugPrint()

class EmptyNonbondedForcefield(NonbondedForcefield):
    def __init__(self):
        super().__init__([], [])

class NonbondedForcefieldFactory:

    @staticmethod
    def createFromStream(fp):
        # assumes a line "$atomtypes" has just been read
        ffType    = 'undefined'
        atomtypes = []
        for line in fp:
            if line[0] == '#':
                continue
            if line.rstrip() == '$end':
                break
            # This if clause is only relevant for the first
            # non-commented line, and sets the forcefield
            # type.
            if (ffType == 'undefined'):
                ffType = line.rstrip()
                continue

            # Create parameters for atomtype.
            if (ffType == 'standard'):
                nb = NonbondedParameters(NonbondedParameterFactory, ['c6', 'c12', 'cs6', 'cs12'], [float(x) for x in line.split()[1:]])
            else:
                raise ValueError("Forcefield type {} is not supported.".format(ffType))

            # Create atomtype.
            at = Atomtype(line.split()[0], nb)
            atomtypes.append(at)
            
        # After break.
        # Create pairtypes (standard and 1-4) based on the atomtypes.
        pts = PairtypeFactory.createFromAtomtypes(ffType, atomtypes)
        return NonbondedForcefield(atomtypes, pts)

# ---------------------------------------------------------------------
# BondedForcefield
# ---------------------------------------------------------------------

class BondedForcefield:
    def getBondtypes(self):
        return self.bondtypes

    def getAngletypes(self):
        return self.angletypes

    def getDihedraltypes(self):
        return self.dihedraltypes

    def getImpropertypes(self):
        return self.impropertypes

class EmptyBondedForcefield(BondedForcefield):
    def __init__(self):
        self.bondtypes     = []
        self.angletypes    = []
        self.dihedraltypes = []
        self.impropertypes = []

