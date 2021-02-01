import  pickle
from    datetime  import  datetime

class State:

    def __init__(self,
                 base_workdir=None,
                 grid=None,
                 protocols=None,
                 properties=None,
                 protocolsHash=None,
                 optimizer=None,
                 gaCoverInterface=None,
                 surrogateModelHash=None,
                 subgridHash=None):
        self.base_workdir = base_workdir 
        self.grid = grid 
        self.protocols = protocols 
        self.properties = properties 
        self.protocolsHash = protocolsHash 
        self.optimizer = optimizer 
        self.gaCoverInterface = gaCoverInterface 
        self.surrogateModelHash = surrogateModelHash 
        self.subgridHash = subgridHash 

    def setFromInput(self, inputTuple):
        self.__init__(*inputTuple)

    def setFromBinary(self, binFile):
        newObject = pickle.load(binFile)
        self.__dict__.update(newObject.__dict__)

    def getInitializationState(self):
       return (self.base_workdir,
               self.grid,
               self.protocols,
               self.properties,
               self.protocolsHash,
               self.optimizer,
               self.gaCoverInterface,
               self.surrogateModelHash,
               self.subgridHash)

    def saveToFile(self):
        print(vars(self))
        _fn = self.base_workdir + '/state.' + datetime.now().strftime('%Y-%m-%d_%H:%M:%S') + '.bin'
        _fp = open(_fn, 'wb')
        pickle.dump(self, _fp, pickle.HIGHEST_PROTOCOL)
        _fp.close()

globalState = State()
