# TODO import ParameterGrid
import numpy as np

class gridOptimizer:
    """
    Class responsible for building a score from the property estimates
    and determining the next promising sampling point.
    """
    def __init__ (self, maxSteps=5, percentCutoff=.25):
        self.maxSteps = maxSteps
        self.percentCutoff = percentCutoff
        # Dictionaries indexed by property id's
        self.referenceValues = {}
        self.referenceTolerances = {}
        self.referenceWeights = {}
        # ```stateScores``` is a list of tuples (gridpoint_id, score)
        self.stateScores = []

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
            # In all other cases, the line should satisfy the syntax
            # PROPID    REFERENCEVALUE     WEIGHT      TOLERANCE
            else:
                splittedLine = line.split()
                propertyName = splittedLine[0]
                self.referenceValues[propertyName] = float(splittedLine[1])
                self.referenceWeights[propertyName] = float(splittedLine[2])
                self.referenceTolerances[propertyName] = float(splittedLine[3])

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
                gPScore += self.referenceWeights[prop_id] * ((gPValue - self.referenceValues[prop_id]) / self.referenceValues[prop_id]) ** 2
            gPScore = np.sqrt(gPScore / weightSum)
            self.stateScores.append((gP.id, gPScore))
        self.rankScores()

    def printToFile (self, grid, filename):
        fp = open(filename, "w")
        properties = self.referenceTolerances.keys()
        fp.write("# %3s" % "id")
        for prop in properties:
            fp.write("%12s" % prop)
            fp.write("%12s" % "err")
        fp.write("%16s\n" % "score")
        for x in self.stateScores:
            fp.write("%5d" % x[0])
            for prop in properties:
                propValue = grid[x[0]].get_property_estimate(prop)
                propErr   = grid[x[0]].get_property_err(prop)
                fp.write("%12.4f%12.4f" % (propValue, propErr))
            fp.write("%16.6e\n" % x[1])
        fp.close()

    def plotToPdf (self, grid, pdf_filename):
        from grid_ana import plot_grid_to_file
        import copy
        
        title = "Score"
        cbox_label = ""
        cbox_limits = ()
        cbox_limits_colors = ()
        if (grid.dim == 2):
            # re-sort state scores by index
            data = copy.deepcopy(self.stateScores)
            data.sort(key=lambda x: x[0])
            data = [x[1] for x in data]
            data = np.array(data)
            data = data.reshape (grid.size[0], grid.size[1])
            plot_grid_to_file (pdf_filename, title, grid.xlabel, grid.ylabel, cbox_label, cbox_limits, cbox_limits_colors, data, grid.get_samples_id())
        else:
            warnings.warn("Can only plot scores for 2-D grid.")
            return

    def determineNextSample (self, grid):
        nTested = 0
        properties = self.referenceTolerances.keys()
        for x in self.stateScores:
            # ignore sampled
            if (x[0] in grid.get_samples_id()):
                nTested += 1
                continue
            if (nTested > self.percentCutoff * grid.linear_size):
                return -1
            for prop in properties:
                propErr = grid[x[0]].get_property_err(prop)
                if (propErr > self.referenceTolerances[prop]):
                    return x[0]
            nTested += 1
