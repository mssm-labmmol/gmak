import numpy as np
import pickle
import os
from datetime import datetime
from gmak.cartesiangrid import CartesianGridIterator


class RecordNameNotFound(Exception):
    pass


class RecordMissingName(Exception):
    pass


class RecordBook:
    def __init__(self, names, tags=None, complete=True):
        self.names = names
        self.records = []
        self._current_record = {}
        self.complete = complete
        self.tags = tags

    def register(self):
        if self.complete:
            for name in self.names:
                if name not in self._current_record.keys():
                    raise RecordMissingName(str(name))
        self.records.append(self._current_record)
        self._current_record = {}

    def set_record(self, name, value):
        if name not in self.names:
            raise RecordNameNotFound(str(name))
        self._current_record[name] = value

    def to_dataframe(self, tag=None):
        import pandas as pd
        df = pd.DataFrame.from_records(self.records)
        names = []
        if tag is not None:
            for name, taglist in self.tags.items():
                for itag in taglist:
                    if itag == tag:
                        if name not in names:
                            names.append(name)
            df = df[names]
        try:
            df.columns = pd.MultiIndex.from_tuples(df.columns)
        except:
            pass
        return df

    @classmethod
    def init_for_state(cls, state):
        names = ["grid", "gridpoint", "aidx"]
        tags = {
            "grid": ["idx"],
            "gridpoint": ["idx"],
            "aidx": ["aidx"],
        }
        # parameters
        for i,_ in enumerate(state.grid.parSpaceGen.getParameterNames()):
            # names are ignored and we use ("X", "1"), ("X", "2"), etc.
            name = ("X", str(i+1))
            names.append(name)
            tags[name] = ["X"]
        # properties
        for prop in state.optimizer.getProperties():
            names.append((prop, "mu"))
            names.append((prop, "sigma"))
            names.append((prop, "diff"))
            tags[(prop, "mu")] = ["property"]
            tags[(prop, "sigma")] = ["property"]
            tags[(prop, "diff")] = ["property"]
        # score
        names.append(("score", "mu"))
        tags[("score", "mu")] = ["score"]
        return cls(names, tags)

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
        self.record_book = RecordBook.init_for_state(self)

    def makeStepPropertiesdir(self, optimizer):
        return self.makeDir("{}/step_{}".format(self.makeCurrentWorkdir(), optimizer.getCurrentIteration()))

    def readCurrentSamples(self):
        gridPath = self.makeCurrentWorkdir()
        step = self.optimizer.getCurrentIteration()
        return np.loadtxt(f"{gridPath}/samples_{step}.dat", dtype=int)

    def makeDir(self, dirname):
        os.system("mkdir -p {}".format(os.path.abspath(dirname)))
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

    def update_record_book(self):
        # fill data
        gridshift = self.grid.shifter.get_current_number_of_shifts()
        for linear_idx, tuple_idx in enumerate(CartesianGridIterator(self.grid.indexGrid)):
            self.record_book.set_record("grid", gridshift)
            self.record_book.set_record("gridpoint", linear_idx)
            self.record_book.set_record("aidx", self.grid.shifter.index_transform.transform(tuple_idx))
            for i, x in enumerate(self.grid.parSpaceGen.getParameterValues(linear_idx)):
                self.record_book.set_record(("X", str(i+1)), x)
            for prop in self.optimizer.getProperties():
                self.record_book.set_record((prop, "mu"), self.grid[linear_idx].get_property_estimate(prop))
                self.record_book.set_record((prop, "sigma"), self.grid[linear_idx].get_property_err(prop))
                self.record_book.set_record((prop, "diff"), self.grid[linear_idx].get_property_estimate(prop) - self.optimizer.getReferenceValue(prop))
            self.record_book.set_record(("score", "mu"), self.optimizer.getScore(linear_idx))
            self.record_book.register()

globalState = State()
