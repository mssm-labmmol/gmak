# TODO import ParameterGrid
import numpy as np
import copy

class gridOptimizer:
    """
    Class responsible for building a score from the property estimates
    and determining the next promising sampling point.
    """
    def __init__ (self, maxSteps=5, percentCutoff=.25):
        self.maxSteps = maxSteps
        self.nsteps   = 0
        self.percentCutoff = percentCutoff
        # Dictionaries indexed by property id's
        self.referenceValues = {}
        self.referenceTolerances = {}
        self.referenceWeights = {}
        # ```stateScores``` is a list of tuples (gridpoint_id, score)
        self.stateScores = []

    def reset(self):
        self.nsteps = 0

    def getCurrentIteration(self):
        return self.nsteps

    def readFromStream (self, stream):
        for line in stream:
            if (line[0] == '#'):
                continue
            elif (line.rstrip() == '$end'):
                return
            elif (line.split()[0] == 'maxsteps'):
                self.maxSteps = int(line.split()[1])
            elif (line.split()[0] == 'ncut'):
                self.percentCutoff = float(line.split()[1])
            # In all other cases, the line has syntax
            # PROPID    REFERENCEVALUE     WEIGHT      TOLERANCE
            else:
                splittedLine = line.split()
                propertyName = splittedLine[0]
                self.referenceValues[propertyName] = float(splittedLine[1])
                if (self.referenceValues[propertyName] == 0):
                    raise ValueError("Reference values for properties cannot be zero!")
                self.referenceWeights[propertyName] = float(splittedLine[2])
                self.referenceTolerances[propertyName] = float(splittedLine[3])

    def pushNearest(self, smHash):
        names = list(self.referenceValues.keys())
        for name in names:
            if (smHash[name] != 'mbar'):
                self.referenceValues[name + '_nearest'] = self.referenceValues[name]
                self.referenceWeights[name + '_nearest'] = 0.00
                self.referenceTolerances[name + '_nearest'] = 1.00e+23

    def rankScores (self):
        self.stateScores.sort(key=lambda x: x[1])

    def fillWithScores (self, grid):
        # clean stateScores first
        self.stateScores = []
        for gP in grid.grid_points:
            gPScore = 0.0
            weightSum = 0.0
            for prop_id in self.referenceTolerances.keys():
                gPValue = gP.get_property_estimate(prop_id)
                weightSum += self.referenceWeights[prop_id]
                sumValue = self.referenceWeights[prop_id] * ((gPValue - self.referenceValues[prop_id]) / self.referenceValues[prop_id]) ** 2
                gPScore += sumValue
            gPScore = np.sqrt(gPScore / weightSum)
            self.stateScores.append((gP.id, gPScore))
        self.rankScores()

    def printToFile (self, grid, filename, sorted=True):
        fp = open(filename, "w")
        properties = self.referenceTolerances.keys()
        fp.write("# %3s" % "id")
        for prop in properties:
            fp.write("%16s" % prop)
            fp.write("%8s" % "err")
        fp.write("%16s\n" % "score")
        if sorted:
            localStateScore = self.stateScores
        else:
            localStateScore = copy.deepcopy(self.stateScores)
            localStateScore.sort(key=lambda x: x[0])
        for x in localStateScore:
            fp.write("%5d" % x[0])
            for prop in properties:
                propValue = grid[x[0]].get_property_estimate(prop)
                propErr   = grid[x[0]].get_property_err(prop)
                fp.write("%16.4f%8.4f" % (propValue, propErr))
            fp.write("%16.6e\n" % x[1])
        fp.close()

    def plotToPdf (self, grid, pdf_filename):
        from grid_ana import plot_grid_to_file, plot_1d_to_file
        import copy
        
        title = "Score"
        cbox_label = ""
        cbox_limits = ()
        cbox_limits_colors = ()
        # re-sort state scores by index
        data = copy.deepcopy(self.stateScores)
        data.sort(key=lambda x: x[0])
        data = [x[1] for x in data]
        data = np.array(data)
        if (grid.get_dim() == 1):
            plot_1d_to_file (pdf_filename, title, grid.xlabel, data, grid.get_samples_id())
        elif (grid.get_dim() == 2):
            data = data.reshape (grid.get_size()[0], grid.get_size()[1])
            plot_grid_to_file (pdf_filename, title, grid.xlabel, grid.ylabel, cbox_label, cbox_limits, cbox_limits_colors, data, grid.get_samples_id())
        else:
            warnings.warn("Can only plot scores for 1-D or 2-D grid.")
            return

    def determineNextSample (self, grid, smHash):
        if (self.nsteps >= self.maxSteps):
            return -1
        nTested = 0
        properties = self.referenceTolerances.keys()
        for x in self.stateScores:
            # ignore sampled
            if (x[0] in grid.get_samples_id()):
                nTested += 1
                continue
            if (nTested > self.percentCutoff * grid.get_linear_size()):
                return -1
            for prop in properties:
                # skip properties without weight - this includes those ending with '_nearest'
                if (self.referenceWeights[prop] == 0.0):
                    continue
                if (smHash[prop] == 'mbar'):
                    propErr = grid[x[0]].get_property_err(prop)
                else:
                    propErr = np.abs(grid[x[0]].get_property_estimate(prop) - grid[x[0]].get_property_estimate(prop + '_nearest'))
                if (propErr > self.referenceTolerances[prop]):
                    self.nsteps += 1
                    return x[0]
            nTested += 1
