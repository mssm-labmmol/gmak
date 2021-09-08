# TODO import ParameterGrid
import runcmd
import numpy as np
import copy
import scipy.stats as st
from   logger import *

def computeNumericError(propertyEstimates, propertyErrors, propertyReferences, propertyWeis, confidenceLevel=.95):
    scoreMax = 0.0
    scoreMin = 0.0
    numberOfProperties = len(propertyEstimates)
    K = (confidenceLevel) **(1.0/numberOfProperties)
    for p in range(numberOfProperties):
        Z = st.norm.ppf(1 - 0.5 * (1 - K), scale=propertyErrors[p])
        if (propertyEstimates[p] >= propertyReferences[p]):
            devMax = propertyEstimates[p] + Z - propertyReferences[p]
            devMin = propertyEstimates[p] - Z - propertyReferences[p]
            if (devMin < 0):
                devMin = 0
        else:
            devMax = propertyReferences[p] - propertyEstimates[p] + Z
            devMin = propertyReferences[p] - propertyEstimates[p] - Z
            if (devMin < 0):
                devMin = 0
        scoreMax += propertyWeis[p] * (devMax/propertyReferences[p]) ** 2
        scoreMin += propertyWeis[p] * (devMin/propertyReferences[p]) ** 2
    scoreMax = np.sqrt(scoreMax/np.sum(propertyWeis))
    scoreMin = np.sqrt(scoreMin/np.sum(propertyWeis))
    return (scoreMin, scoreMax)


class gridOptimizer:
    """
    Class responsible for building a score from the property estimates
    and determining the next promising sampling point.
    """
    def __init__ (self, maxSteps=0, percentCutoff=1.0):
        self.maxSteps = maxSteps
        self.nsteps   = 0
        self.percentCutoff = percentCutoff
        self.confidenceLevels = [0.68, 0.80, 0.90, 0.95]
        # Dictionaries indexed by property id's
        self.referenceValues = {}
        self.referenceTolerances = {}
        self.referenceWeights = {}
        # ```stateScores``` is a list of tuples (gridpoint_id, score)
        self.stateScores = []

    def getProperties(self):
        return list(self.referenceValues.keys())

    def getTolerances(self):
        return self.referenceTolerances

    def reset(self):
        self.nsteps = 0

    def getCurrentIteration(self):
        return self.nsteps

    def readFromStream (self, stream, validateFlag):
        for line in stream:
            if (line[0] == '#'):
                continue
            elif (line.rstrip() == '$end'):
                return
            elif (line.split()[0] == 'maxsteps'):
                if (validateFlag):
                    self.maxSteps = 0
                else:
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

    def rankScores (self):
        self.stateScores.sort(key=lambda x: x[1])
        self.stateScoreIntervals[:] = [self.stateScoreIntervals[x[0]] for x in self.stateScores]

    def getRankedBest(self, n):
        try:
            return [self.stateScores[i][0] for i in range(n)]
        except IndexError:
            return self.stateScores[:][0] 

    def fillWithScores (self, grid):
        # Clean stateScores first
        self.stateScores = []
        # self.stateScoreIntervals[i][j] is j-th confidence interval for i-th gridpoint.
        self.stateScoreIntervals = [[-1, [(0.0, 0.0) for j in range(len(self.confidenceLevels))]] for i in range(grid.get_linear_size())]
        for gP in grid.grid_points:
            propertyEstimates = []
            propertyErrs = []
            propertyReferences = []
            propertyWeis = []
            gPScore = 0.0
            weightSum = 0.0
            for prop_id in self.referenceTolerances.keys():
                # ignore properties with no weight
                if self.referenceWeights[prop_id] == 0.0:
                    continue
                gPValue = gP.get_property_estimate(prop_id)
                propertyEstimates.append(gPValue)
                propertyErrs.append(gP.get_property_err(prop_id))
                propertyReferences.append(self.referenceValues[prop_id])
                propertyWeis.append(self.referenceWeights[prop_id])
                weightSum += self.referenceWeights[prop_id]
                sumValue = self.referenceWeights[prop_id] * ((gPValue - self.referenceValues[prop_id]) / self.referenceValues[prop_id]) ** 2
                gPScore += sumValue
            gPScore = np.sqrt(gPScore / weightSum)
            self.stateScores.append((gP.id, gPScore))
            self.stateScoreIntervals[gP.id][0] = gP.id
            for ci, cl in enumerate(self.confidenceLevels):
                self.stateScoreIntervals[gP.id][1][ci] = computeNumericError(propertyEstimates, propertyErrs, propertyReferences, propertyWeis, confidenceLevel=cl)
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
            localIntervals = self.stateScoreIntervals
        else:
            localStateScore = copy.deepcopy(self.stateScores)
            localStateScore.sort(key=lambda x: x[0])
            localIntervals = copy.deepcopy(self.stateScoreIntervals)
            localIntervals.sort(key=lambda x: x[0])
        for x in localStateScore:
            fp.write("%5d" % x[0])
            for prop in properties:
                propValue = grid[x[0]].get_property_estimate(prop)
                propErr   = grid[x[0]].get_property_err(prop)
                fp.write("%16.4f%8.4f" % (propValue, propErr))
            fp.write("%16.6e\n" % x[1])
        fp.close()
        fp = open(filename + '.cis', "w")
        fp.write("# Scores with Confidence Intervals\n")
        fp.write("#%15s%16s" % ("id", "central"))
        for cl in self.confidenceLevels:
            for spec in ["_min", "_max"]:
                fp.write("%16s" % (str(int(100*cl)) + spec))
        fp.write("\n")
        for x, s in zip(localIntervals, localStateScore):
            fp.write("%16d" % x[0])
            fp.write("%16.6e" % s[1])
            for y in x[1]:
                fp.write("%16.6e%16.6e" % (y[0], y[1]))
            fp.write("\n")
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

    def determineNextSample(self, grid, smHash, protHash):
        globalLogger.putMessage('BEGIN OPTIMIZER')
        globalLogger.indent()
        properties = self.referenceTolerances.keys()
        nTested = 0
        outSamples = []
        properties = self.referenceTolerances.keys()
        if (self.nsteps >= self.maxSteps):
            globalLogger.putMessage('MESSAGE: Estimates are all converged.')
            return -1
        for x in self.stateScores:
            # ignore sampled
            if (x[0] in grid.get_samples_id()):
                nTested += 1
                continue
            if (nTested > self.percentCutoff * grid.get_linear_size()):
                globalLogger.putMessage('MESSAGE: All {} best GridPoints are OK!'
                                        .format(int(self.percentCutoff * grid.get_linear_size())))
                globalLogger.unindent()
                globalLogger.putMessage('END OPTIMIZER')
                return -1
            for prop in properties:
                # skip properties without weight
                if (self.referenceWeights[prop] == 0.0):
                    continue
                if (smHash[prop] == 'mbar'):
                    propErr = grid[x[0]].get_property_err(prop)
                    if (propErr > self.referenceTolerances[prop]):
                        self.nsteps += 1
                        outSamples.append(x[0])
                        globalLogger.putMessage('MESSAGE: GridPoint {} will be simulated next.'.format(x[0]))
                        globalLogger.unindent()
                        globalLogger.putMessage('END OPTIMIZER')
                        return [x[0]]
                else: # i.e., interpolated properties
                    propErr = grid[x[0]].get_property_err(prop)
                    if (propErr > self.referenceTolerances[prop]):
                        self.nsteps += 1
                        outSamples.append(x[0])
                        globalLogger.putMessage('MESSAGE: GridPoint {} will be simulated next.'.format(x[0]))
                        globalLogger.unindent()
                        globalLogger.putMessage('END OPTIMIZER')
                        return [x[0]]
            self.nsteps += 1
            nTested += 1
