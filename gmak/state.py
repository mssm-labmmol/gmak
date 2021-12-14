import numpy as np
import  pickle
from os import system
from    datetime  import  datetime
import  os

class State:

    def __init__(self,
                 base_workdir=None,
                 grid=None,
                 protocols=None,
                 properties=None,
                 protocolsHash=None,
                 optimizer=None,
                 surrogateModelHash=None):
        self.base_workdir = base_workdir
        self.grid = grid
        self.protocols = protocols
        self.properties = properties
        self.protocolsHash = protocolsHash
        self.optimizer = optimizer
        self.surrogateModelHash = surrogateModelHash

    def merge(self, other):
        self.grid.merge(other.grid)
        for pi, pj in zip(self.protocols, other.protocols):
            pi.merge(pj)
        self.optimizer.merge(other.optimizer)

    def makeCurrentWorkdir(self):
        shifts = 0
        outdir = "{}/grid_{}".format(self.base_workdir, shifts)
        while(os.path.exists(outdir)):
            shifts += 1
            outdir = "{}/grid_{}".format(self.base_workdir, shifts)
        shifts -= 1
        outdir = "{}/grid_{}".format(self.base_workdir, shifts)
        return self.makeDir(outdir)

    def resetWorkdir(self, workdir):
        self.base_workdir = workdir

    def setFromInput(self, inputTuple):
        self.__init__(*inputTuple)

    def makeStepPropertiesdir(self, optimizer):
        return self.makeDir("{}/step_{}".format(self.makeCurrentWorkdir(), optimizer.getCurrentIteration()))

    def readCurrentSamples(self):
        gridPath = self.makeCurrentWorkdir()
        step = self.optimizer.getCurrentIteration()
        return np.loadtxt(f"{gridPath}/samples_{step}.dat", dtype=int)

    def makeDir(self, dirname):
        system("mkdir -p {}".format(os.path.abspath(dirname)))
        return os.path.abspath(dirname)

    def setFromBinary(self, binFile):
        newObject = pickle.load(binFile)
        self.__dict__.update(newObject.__dict__)

    def get_samples_id (self):
        self.readCurrentSamples()

    def getInitializationState(self):
       return (self.base_workdir,
               self.grid,
               self.protocols,
               self.properties,
               self.protocolsHash,
               self.optimizer,
               self.surrogateModelHash)

    def saveToFile(self):
        #print(vars(self))
        _fn = self.base_workdir + '/state_' + str(os.getpid()) + '.bin'
        _fp = open(_fn, 'wb')
        pickle.dump(self, _fp, pickle.HIGHEST_PROTOCOL)
        _fp.close()

globalState = State()
