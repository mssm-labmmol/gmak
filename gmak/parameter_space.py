from abc import ABC, abstractmethod
import numpy as np
from gmak.variations import DomainSpace
from gmak.interaction_parameter import InteractionParameter
from math import sqrt
from copy import deepcopy
import re
from typing import List

class ParameterSpaceGenerator:
    """
    Contains tuples of (DomainSpace, List[InteractionParameter]) indexed in a
    dictionary.
    """
    def __init__(self):
        self.members = {}

    def addMember(self,
                  name: str,
                  domainSpace: DomainSpace,
                  params: List[InteractionParameter]):
        """
        Adds a member to the ParameterSpaceGenerator object.
        """
        # check compatibility of number of states
        for ds, pl in self.members.values():
            if (ds.get_linear_size() != domainSpace.get_linear_size()):
                raise ValueError("Trying to add a DomainSpace with non-matching linear size.")
        # check compatibility of dimension
        if (domainSpace.get_dim() != len(params)):
            raise ValueError("Trying to associate DomainSpace of dimension {} to {} parameters.".format(
                domainSpace.get_dim(), len(params)))
        self.members[name] = (domainSpace, params)

    def setState(self, i):
        """
        Applies the variation values at state i to the interaction parameters.
        """
        for domainSpace, parList in self.members.values():
            values = domainSpace.get(i)
            for j, par in enumerate(parList):
                par.value = values[j]

    def getStateParameters(self, i):
        """
        Gets parameter values at state i without applying them to the parameter
        reference -> returns a modified copy of the original parameters.
        """
        output_parameters = []
        for domainSpace, parList in self.members.values():
            values = domainSpace.get(i)
            for j, par in enumerate(parList):
                par_copy = deepcopy(par)
                par_copy.value = values[j]
                output_parameters.append(par_copy)
        return output_parameters

    def getParameterNames(self):
        """
        Returns the list of names of the interaction parameters associated with
        the main variation.
        """
        return [x.name for x in self.members['main'][1]]

    def getParameterReferences(self):
        out = []
        for x in self.members.keys():
            out += [y for y in self.members[x][1]]
        return out

    def getParameterValues(self, state):
        return self.members['main'][0].get(state)

    def getAllParameterValues(self):
        return self.members['main'][0].get_data()

    def writeParameters(self, prefix):
        for name in self.members.keys():
            domainSpace, parList = self.members[name]
            outfile = prefix + "_" + name + ".dat"
            domainSpace.write_to_file(outfile)

    def setNewCenter(self, i):
        for name in self.members:
            self.members[name][0].set_new_center(i)

    def setNewOrigin(self, i):
        for name in self.members:
            self.members[name][0].set_new_origin(i)

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

    def writeMainVariationBlock(self, stream):
        """
        This method looks a bit out of place.
        """
        name = 'main'
        stream.write("$variation\n")
        # Writing name and parameter strings requires ParameterSpaceGenerator
        stream.write("name {}\n".format(name))
        atomtype_strings = [x.name for x in self.members[name][1]]
        stream.write("pars " + " ".join(atomtype_strings) + "\n")
        # Writing the rest is delegated to DomainSpace
        self.members[name][0].write_block_to_stream(stream)
        stream.write("$end\n\n")
