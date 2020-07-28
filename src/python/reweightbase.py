import numpy as np
from abc import ABC, abstractmethod
import reweight
import traj_ana
from shutil import copyfile
import os
from collections import OrderedDict

class GromacsReweighterAdapter:
    @staticmethod
    def getInputData(originGridpoint, destGridpoint, protocol):
        properties = protocol.get_reweighting_properties()
        outputDict = {}
        for prop in properties:
            outputDict[prop] = originGridpoint.retrieve_atomic_property_from_protocol(prop, protocol)
        outputDict['gro'] = originGridpoint.protocol_outputs[protocol.name]['gro']
        outputDict['top'] = destGridpoint.protocol_outputs[protocol.name]['top'] # This is really destGridpoint!
        outputDict['mdp'] = protocol.mdps[-1]
        if (protocol.type == 'slab'):
            outputDict['xtc'] = originGridpoint.protocol_outputs[protocol.name]['trr']
        else:
            outputDict['xtc'] = originGridpoint.protocol_outputs[protocol.name]['xtc']
        return outputDict
                
class ReweighterInterface(ABC):

    @abstractmethod
    def runSimulations(self, protocol, workdir): pass

    @abstractmethod
    def calculateProperties(self, protocol, workdir): pass

    def getConfigurationMatrix(self):
        return np.array(self._configurationMatrix)

    def getPropertyMatrix(self, property):
        return np.array(self._propertyMatrix[property])

    def getFullPropertyMatrix(self):
        return np.array(list(self._propertyMatrix.values()))

    def run(self, protocol, workdir):
        self.runSimulations(protocol, workdir)
        self.calculateProperties(protocol, workdir)
        #self._cleanFiles()

class EmptyReweighter(ReweighterInterface):
    def __init__(self):
        return
    def runSimulations(self, protocol, workdir):
        return
    def calculateProperties(self, protocol, workdir):
        return
    def getConfigurationMatrix(self):
        raise NotImplementedError("Trying to extract data from EmptyReweighter.")
    def getPropertyMatrix(self, property):
        raise NotImplementedError("Trying to extract data from EmptyReweighter.")
    def getFullPropertyMatrix(self):
        raise NotImplementedError("Trying to extract data from EmptyReweighter.")
    def run(self, protocol, workdir):
        return

class StandardReweighterInterface(ReweighterInterface):

    def __init__(self, parameterGrid, rerunFunction, analyzeFunction, reweightAdapter):
        # State that can not be directly altered by mehods. 
       self.parameterGrid = parameterGrid
       self.rerun         = rerunFunction
       self.analyze       = analyzeFunction
       self.adapter       = reweightAdapter
       # State that can be altered by methods.
       self._outputFiles         = [[None for j in range(parameterGrid.get_linear_size())]
                                    for i in range(parameterGrid.get_linear_size())]
       self._configurationMatrix = [0 for i in range(parameterGrid.get_linear_size())]

    def _initProperties(self, protocol):
        self.properties = protocol.get_reweighting_properties()
        self._propertyMatrix = OrderedDict()
        for prop in self.properties:
            self._propertyMatrix[prop] = [[] for i in range(self.parameterGrid.get_linear_size())]
        self._configurationMatrix = [0 for i in range(self.parameterGrid.get_linear_size())]

    def _cleanFiles(self):
        # TODO: This function needs to be fixed before using.
        #       'xj' is a dict, not a string
        for xi in self._outputFiles:
            for xj in xi:
                if (xj is not None):
                    os.remove(xj)

    def runSimulations(self, protocol, workdir):
        for stateOrigin in self.parameterGrid.get_samples_id():
            for stateDest in range(self.parameterGrid.get_linear_size()):
                pointOrigin = self.parameterGrid[stateOrigin]
                pointDest   = self.parameterGrid[stateDest]
                _input      = self.adapter.getInputData(pointOrigin, pointDest, protocol)
                _dir        = "{}/{}_{}/".format(workdir, stateOrigin, stateDest)
                _output     = self.rerun(_input, _dir)
                self._outputFiles[stateOrigin][stateDest] = _output

    def _calculatePropertiesError(self, inputData, rwData, prop, out):
        raise NotImplemented("You are trying to reweight a property that has no strategy to be reweighted.")

    def _calculatePropertiesCopy(self, inputData, rwData, prop, out):
        copyfile(inputData[prop], out)
        return list(np.loadtxt(out, comments=['@','#'], usecols=(1,)))
                
    def _calculatePropertiesAnalyze(self, inputData, rwData, prop, out):
        return list(self.analyze(rwData, prop, out))
    
    def calculateProperties(self, protocol, workdir):
        func_dict = {
            'potential' : self._calculatePropertiesAnalyze,
            'density'   : self._calculatePropertiesCopy,
            'pV'        : self._calculatePropertiesCopy,
            'gamma'     : self._calculatePropertiesError,
            'volume'    : self._calculatePropertiesCopy,
            'polcorr'   : self._calculatePropertiesCopy}
        
        self._initProperties(protocol)
        for prop in self.properties:
            for stateOrigin in self.parameterGrid.get_samples_id():
                for stateDest in range(self.parameterGrid.get_linear_size()):
                    pointOrigin  = self.parameterGrid[stateOrigin]
                    pointDest    = self.parameterGrid[stateDest]
                    _input       = self.adapter.getInputData(pointOrigin, pointDest, protocol)
                    _rwdata      = self._outputFiles[stateOrigin][stateDest]
                    _out         = "{}/{}_{}/{}.xvg".format(workdir, stateOrigin, stateDest, prop)
                    propertyData = func_dict[prop](_input, _rwdata, prop, _out)
                    self._configurationMatrix[stateOrigin] = len(propertyData)
                    self._propertyMatrix[prop][stateDest] += propertyData

class GromacsStandardReweighter(StandardReweighterInterface):
    def __init__(self, parameterGrid):
        super().__init__(parameterGrid, reweight.reweightWrapper, traj_ana.analyzeWrapper, GromacsReweighterAdapter)

class ReweighterFactory:
    @staticmethod
    def create(typeString, parameterGrid):
        if (typeString == 'fast'):
            raise NotImplemented("No fast method implemented.")
        elif (typeString == 'standard'):
            return GromacsStandardReweighter(parameterGrid)
