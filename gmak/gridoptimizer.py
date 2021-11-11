# TODO import ParameterGrid
import gmak.runcmd as runcmd
import numpy as np
import copy
import scipy.stats as st
import gmak.logger as logger

def wrmsd_calc_score_err(propertyEstimates,
                         propertyErrors,
                         propertyWeis,
                         propertyReferences,
                         confidenceLevel=.95):
    keys = list(propertyEstimates.keys())
    numberOfProperties = len(keys)
    ests = [propertyEstimates[p] for p in keys]
    refs = [propertyReferences[p] for p in keys]
    weis = [propertyWeis[p] for p in keys]
    errs = [propertyErrors[p] for p in keys]

    scoreMax = 0.0
    scoreMin = 0.0
    K = (confidenceLevel) **(1.0/numberOfProperties)
    for p in range(numberOfProperties):
        Z = st.norm.ppf(1 - 0.5 * (1 - K), scale=errs[p])
        if (ests[p] >= refs[p]):
            devMax = ests[p] + Z - refs[p]
            devMin = ests[p] - Z - refs[p]
            if (devMin < 0):
                devMin = 0
        else:
            devMax = refs[p] - ests[p] + Z
            devMin = refs[p] - ests[p] - Z
            if (devMin < 0):
                devMin = 0
        scoreMax += weis[p] * (devMax) ** 2
        scoreMin += weis[p] * (devMin) ** 2
    scoreMax = np.sqrt(scoreMax/np.sum(weis))
    scoreMin = np.sqrt(scoreMin/np.sum(weis))
    return (scoreMin, scoreMax)


def wrmsd_calc_score(propertyEstimates,
                     propertyErrs,
                     propertyWeis,
                     propertyReferences):
    keys = list(propertyEstimates.keys())
    ests = [propertyEstimates[p] for p in keys]
    refs = [propertyReferences[p] for p in keys]
    weis = [propertyWeis[p] for p in keys]
    score  = (np.array(ests) - np.array(refs)) ** 2
    score *= np.array(weis)
    score  = np.sqrt(np.sum(score)/np.sum(weis))
    return score


class ScoreNotInitialized(Exception):
    pass

class ScoreFactory:

    ptable = {}

    @classmethod
    def add_custom_score(cls, name, score, score_err):
        cls.ptable[name] = (score, score_err)

    @classmethod
    def create(cls, name):
        try:
            return cls.ptable[name]
        except KeyError:
            raise ScoreNotInitialized(f"Could not initialize score function {name}.")



class gridOptimizer:
    """
    Class responsible for building a score from the property estimates
    and determining the next promising sampling point.
    """
    def __init__ (self,
                  maxSteps,
                  percentCutoff,
                  refs,
                  tols,
                  weis,
                  score,
                  score_err,
                  confidence_levels):
        self.maxSteps = maxSteps
        self.nsteps   = 0
        self.percentCutoff = percentCutoff
        # Dictionaries indexed by property id's
        self.referenceValues = refs
        self.referenceTolerances = tols
        self.referenceWeights = weis
        # ```stateScores``` is a list of tuples (gridpoint_id, score)
        self.stateScores = []
        self.score = score
        self.score_err = score_err
        self.confidenceLevels = confidence_levels

    @classmethod
    def from_dict(cls, bd, validateFlag):
        reserved_keys = ['function', 'maxsteps', 'ncut']
        if 'function' not in bd.keys():
            _score = wrmsd_calc_score
            _score_err = wrmsd_calc_score_err
            _conf_ls = [0.68, 0.80, 0.90, 0.95]
        else:
            # then choose from custom types
            _score, _score_err = ScoreFactory.create(bd['function'][0])
            _conf_ls = None
        # mount dicts
        _refs = {}
        _tols = {}
        _weis = {}
        for key in bd.keys():
            if key not in reserved_keys:
                prop_id = key
                _refs[prop_id] = float(bd[key][0])
                _weis[prop_id] = float(bd[key][1])
                _tols[prop_id] = float(bd[key][2])

        return cls(maxSteps=bd['maxsteps'][0],
                   percentCutoff=bd['ncut'][0],
                   refs=_refs,
                   tols=_tols,
                   weis=_weis,
                   score=_score,
                   score_err=_score_err,
                   confidence_levels=_conf_ls)

    def getProperties(self):
        return list(self.referenceValues.keys())

    def getTolerances(self):
        return self.referenceTolerances

    def reset(self):
        self.nsteps = 0

    def getCurrentIteration(self):
        return self.nsteps

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
        # self.stateScoreIntervals[i][j] is j-th confidence interval for i-th
        # gridpoint.
        if self.confidenceLevels is not None:
            self.stateScoreIntervals = [
                [-1, [(0.0, 0.0) for j in range(len(self.confidenceLevels))]]
                for i in range(grid.get_linear_size())]
        else:
            self.stateScoreIntervals = [[-1, (0.0, 0.0)]
                                        for i in
                                        range(grid.get_linear_size())]

        properties = self.getProperties()
        for gP in grid.grid_points:
            ests = {p: gP.get_property_estimate(p) for p in properties}
            errs = {p: gP.get_property_err(p) for p in properties}
            # score
            gPScore = self.score(ests,
                                 errs,
                                 self.referenceWeights,
                                 self.referenceValues)
            self.stateScores.append((gP.id, gPScore))
            # score errors
            self.stateScoreIntervals[gP.id][0] = gP.id
            if self.score_err is not None:
                if self.confidenceLevels is None:
                        self.stateScoreIntervals[gP.id][1] = self.score_err(
                            ests,
                            errs,
                            self.referenceWeights,
                            self.referenceValues)
                else:
                    for ci, cl in enumerate(self.confidenceLevels):
                        self.stateScoreIntervals[gP.id][1][ci] = self.score_err(
                            ests,
                            errs,
                            self.referenceWeights,
                            self.referenceValues,
                            confidenceLevel=cl)
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
        if self.score_err is not None:
            fp = open(filename + '.cis', "w")
            fp.write("# Scores with Uncertainties\n")
            fp.write("#%15s%16s" % ("id", "central"))
            if self.confidenceLevels is not None:
                for cl in self.confidenceLevels:
                    for spec in ["_min", "_max"]:
                        fp.write("%16s" % (str(int(100*cl)) + spec))
            else:
                for spec in ["min", "max"]:
                    fp.write("%16s" % spec)
            fp.write("\n")
            for x, s in zip(localIntervals, localStateScore):
                fp.write("%16d" % x[0])
                fp.write("%16.6e" % s[1])
                if self.confidenceLevels is not None:
                    for y in x[1]:
                        fp.write("%16.6e%16.6e" % (y[0], y[1]))
                else:
                    fp.write("%16.6e%16.6e" % (x[1][0], x[1][1]))
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
        logger.globalLogger.putMessage('BEGIN OPTIMIZER')
        logger.globalLogger.indent()
        properties = self.referenceTolerances.keys()
        nTested = 0
        outSamples = []
        properties = self.referenceTolerances.keys()
        if (self.nsteps >= self.maxSteps):
            logger.globalLogger.putMessage('MESSAGE: Estimates are all converged.')
            return -1
        for x in self.stateScores:
            # ignore sampled
            if (x[0] in grid.get_samples_id()):
                nTested += 1
                continue
            if (nTested > self.percentCutoff * grid.get_linear_size()):
                logger.globalLogger.putMessage('MESSAGE: All {} best GridPoints are OK!'
                                        .format(int(self.percentCutoff * grid.get_linear_size())))
                logger.globalLogger.unindent()
                logger.globalLogger.putMessage('END OPTIMIZER')
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
                        logger.globalLogger.putMessage('MESSAGE: GridPoint {} will be simulated next.'.format(x[0]))
                        logger.globalLogger.unindent()
                        logger.globalLogger.putMessage('END OPTIMIZER')
                        return [x[0]]
                else: # i.e., interpolated properties
                    propErr = grid[x[0]].get_property_err(prop)
                    if (propErr > self.referenceTolerances[prop]):
                        self.nsteps += 1
                        outSamples.append(x[0])
                        logger.globalLogger.putMessage('MESSAGE: GridPoint {} will be simulated next.'.format(x[0]))
                        logger.globalLogger.unindent()
                        logger.globalLogger.putMessage('END OPTIMIZER')
                        return [x[0]]
            self.nsteps += 1
            nTested += 1


def add_custom_score(type_name, calc_score, calc_score_err=None):
    """
    Adds a custom score function to the program. In the input file, it can be
    referenced with the type ``type_name``.

    :param type_name: Name of the type of the custom function
    :type type_name: str
    :param calc_score: The score function (see :py:func:`~gmak.custom_scores.calc_score`)
    :type calc_score: callable
    :param calc_score_err: (optional) The score uncertainty function (see
        :py:func:`~gmak.custom_scores.calc_score_err`). Defaults to
        :py:obj:`None`, which means that the uncertainties are not calculated.
    :type calc_score_err: callable
    """
    ScoreFactory.add_custom_score(type_name, calc_score, calc_score_err)
