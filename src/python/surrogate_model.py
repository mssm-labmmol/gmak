"""This file defines the base class SurrogateModel and the derived
classes for each type of surrogate model: MBAR, the three
interpolation methods nearest neighbor, linear and (bi)cubic spline
interpolation, and Gaussian Process Regression.
"""

# TODO: Rewrite it so that there are two types of SurrogateModel:
# those based on Reweight (MBAR) and those based on Avgs and Stds (all
# interpolations + GPR).
#
# TODO: The AvgStd stuff can be in the base class.

from os import system
import pymbar
import numpy as np
import sys
from scipy.interpolate import griddata
# Import scikit-learn things
from sklearn import gaussian_process
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel, Matern, DotProduct, ExpSineSquared, RationalQuadratic

def init_surrogate_model_from_string (string, bool_legacy):
    if(string == 'mbar'):
        if (bool_legacy):
            return MBARLegacy()
        else:
            return MBAR()
    elif(string == 'gpr'):
        return GaussianProcessRegressionInterpolation()
    elif(string == 'empty'):
        return EmptySurrogateModel()
    else:
        return Interpolation(string)

class SurrogateModel:

    def requiresCorners(self):
        """
        Is it necessary to sample the grid corners?
        """
        return self.corners

    def requiresReweight(self):
        """
        Is it necessary to reweight potential energy and mechanical properties?
        """
        return self.reweight

    def writeExpectationsToFile(self, fn_avg, fn_err, which_property):
        np.savetxt(fn_avg, self.EA_pk[which_property,:])
        np.savetxt(fn_err, self.dEA_pk[which_property,:])

    def writeLogToDirectory(self, dir_path):
        return

class EmptySurrogateModel(SurrogateModel):

    corners = False
    reweight = False
    kind = 'empty'

    def __init__(self):
        return

    def _computeAvgStd(self, A_psn):
        numberProps = len(A_psn)
        numberStates = len(A_psn[0])
        A_ps = []
        dA_ps = []
        for i in range(numberProps):
            A_s = []
            dA_s = []
            for s in range(numberStates):
                numberOfConfs = len(A_psn[i][s])
                if (numberOfConfs < 2):
                    raise Exception
                elif (numberOfConfs == 2):
                    # In this case, first line is average and second line is error.
                    A_s.append(A_psn[i][s][0])
                    dA_s.append(A_psn[i][s][1])
                else:
                    A_s.append(np.mean(A_psn[i][s]))
                    dA_s.append(np.std(A_psn[i][s], ddof=1) / np.sqrt(len(A_psn[i][s])))
            A_ps.append(A_s)
            dA_ps.append(dA_s)
        A_ps = np.array(A_ps)
        dA_ps = np.array(dA_ps)
        return (A_ps, dA_ps)

    def computeExpectationsFromAvgStd(self, A_ps, dA_ps, I_s, gridShape):
        # make grid from shape
        gridDomain = np.meshgrid(*[np.arange(s) for s in gridShape], indexing='ij')
        # determine grid indices for sampled states
        sampleIndices = []
        linearIdx = 0
        for idx in np.ndindex(gridShape):
            if linearIdx in I_s:
                sampleIndices.append(idx)
            linearIdx += 1
        # interpolate property data
        A_pk = []
        dA_pk = []
        numberProps = A_ps.shape[0]
        for i in range(numberProps):
            A_k  = np.zeros(tuple(gridShape))
            dA_k = np.zeros(tuple(gridShape))
            for sr, si in enumerate(sampleIndices):
                A_k[si] = A_ps[i,sr]
                dA_k[si] = dA_ps[i,sr]
            A_pk.append(A_k)
            dA_pk.append(dA_k)
        A_pk = np.array(A_pk)
        dA_pk = np.array(dA_pk)
        A_pk = A_pk.reshape((numberProps, A_k.size))
        dA_pk = dA_pk.reshape((numberProps, dA_k.size))        
        self.EA_pk = A_pk
        self.dEA_pk = dA_pk
        return (A_pk, dA_pk)
        
    def computeExpectations(self, A_psn, I_s, gridShape):
        """
        Parameters:
        -----------
        A_psn :    bi-dimensional list P x S of (np.ndarray, float, shape N(S))
                      P - number of properties
                      S - number of sampled states
                      N - number of uncorrelated configurations
                      A_psn[p,s,n] - property 'p' for configuration 'n' of sampled state 's'

        I_s   :    np.ndarray of shape (S)
                      I_s[s] - linear index of sampled state 's'
        
        gridShape : tuple (n_1, n_2, ..., n_D), int
                      n_d - number of partitions for edge 'd'

        Returns:
        --------
        (EA_pk, dEA_pk) :    EA_pk :    np.ndarray, float, shape (P,K)
                            dEA_pk :    np.ndarray, float, shape (P,K)

        """
        # mean and std for sampled states
        (A_ps, dA_ps) = self._computeAvgStd(A_psn)
        return self.computeExpectationsFromAvgStd(A_ps, dA_ps, I_s, gridShape)

class MBAR (SurrogateModel):

    """Note: legacy is False for MBAR, but true for the subclass MBARLegacy."""
    legacy = False
    reweight = True
    corners = False
    kind = 'mbar'

    def __init__(self):
        self.Eff_k = []
        self.W_k = []
        self.EA_pk = []
        self.dEA_pk = []
        return

    def computeExpectations(self, A_pkn, u_kn, N_k): 
        """
        Computes expectations (mean values and uncertainties) for the properties.
        These values are returned as a tuple and also stored internally.

        Parameters:
        -----------
        A_pkn :    np.ndarray, float, shape (P, K, N)
                      P - number of properties
                      K - number of states
                      N - number of uncorrelated configurations
                      A_pkn[p,k,n] - property 'p' for configuration 'n' at state 'k'

        u_kn  :    np.ndarray, float, shape (K, N)
                      u_kn[k,n] - reduced energy for configuration 'n' at state 'k'

        N_k   :    np.ndarray of shape (K)
                      N_k[k] - number of reduced energy samples collected from state 'k'

        Returns:
        --------
        (EA_pk, dEA_pk) :    EA_pk :    np.ndarray, float, shape (P,K)
                            dEA_pk :    np.ndarray, float, shape (P,K)

        """
        n_props = A_pkn.shape[0]
        n_states = A_pkn.shape[1]
        n_confs = A_pkn.shape[2]
        self.EA_pk = np.zeros((n_props, n_states))
        self.dEA_pk = np.zeros((n_props, n_states))

        # consistency check
        if (n_states != u_kn.shape[0]):
            raise ValueError("Number of states inconsistent between properties and reduced energies.")
        if (n_states != N_k.shape[0]):
            raise ValueError("Number of states inconsistent between properties and collected samples.")

        mbar = pymbar.MBAR(u_kn, N_k)
        self.Eff_k = mbar.computeEffectiveSampleNumber()
        self.W_k = mbar.getWeights()
        self.N_k = N_k

        # loop over properties
        for i in range(n_props):
            mbarComputeOutput = mbar.computeExpectations(A_pkn[i,:,:], state_dependent=True)
            # MBAR API is different depending on ``legacy'' option
            if (self.legacy):
                self.EA_pk[i,:] = mbarComputeOutput[0]
                self.dEA_pk[i,:] = mbarComputeOutput[1]
            else:
                self.EA_pk[i,:] = mbarComputeOutput['mu']
                self.dEA_pk[i,:] = mbarComputeOutput['sigma']

        return (self.EA_pk, self.dEA_pk)

    def writeLogToDirectory(self, dir_path):
        """
        Parameters:
        -----------
        dir_path :    path of directory where files will be stored

        The output files are:

        dir_path/Eff_k.dat :    efficiencies of the estimates
        dir_path/N_k.dat   :    number of samples collected from each state
        dir_path/W_k.dat   :    MBAR weights 

        If not existant, dir_path is created.
        """
        system('mkdir -p ' + dir_path)
        np.savetxt(dir_path + '/Eff_k.dat', self.Eff_k)
        np.savetxt(dir_path + '/N_k.dat', self.N_k)
        np.savetxt(dir_path + '/W_k.dat', self.W_k)

class MBARLegacy (MBAR):
    """
    A specific class for MBAR legacy which differs from MBAR by the value of the
    'legacy' attribute. All other attributes and methods are inherited from the
    base MBAR class. Note that even the principal ``computeExpectations`` method
    is inherited from MBAR.
    """
    legacy = True

class Interpolation (SurrogateModel):

    reweight = False
    allowed_kinds = ['nearest', 'linear', 'cubic']

    def __init__(self, kind): 
        if kind in self.allowed_kinds:
            self.kind = kind
            if (kind == 'nearest'):
                self.corners = False
            else:
                self.corners = True
        else:
            raise ValueError ("Unknown interpolation type: {}.".format(kind))

    def _computeAvgStd(self, A_psn):
        numberProps = len(A_psn)
        numberStates = len(A_psn[0])
        A_ps = []
        dA_ps = []
        for i in range(numberProps):
            A_s = []
            dA_s = []
            for s in range(numberStates):
                numberOfConfs = len(A_psn[i][s])
                if (numberOfConfs < 2):
                    raise Exception
                elif (numberOfConfs == 2):
                    # In this case, first line is average and second line is error.
                    A_s.append(A_psn[i][s][0])
                    dA_s.append(A_psn[i][s][1])
                else:
                    A_s.append(np.mean(A_psn[i][s]))
                    dA_s.append(np.std(A_psn[i][s], ddof=1) / np.sqrt(len(A_psn[i][s])))
            A_ps.append(A_s)
            dA_ps.append(dA_s)
        A_ps = np.array(A_ps)
        dA_ps = np.array(dA_ps)
        return (A_ps, dA_ps)

    def computeExpectationsFromAvgStd(self, A_ps, dA_ps, I_s, gridShape):
        # make grid from shape
        gridDomain = np.meshgrid(*[np.arange(s) for s in gridShape], indexing='ij')
        # determine grid indices for sampled states
        sampleIndices = []
        linearIdx = 0
        for idx in np.ndindex(gridShape):
            if linearIdx in I_s:
                sampleIndices.append(idx)
            linearIdx += 1
        # interpolate property data
        A_pk = []
        dA_pk = []
        numberProps = A_ps.shape[0]
        for i in range(numberProps):
            A_k  = griddata(sampleIndices,  A_ps[i,:], tuple(gridDomain), method=self.kind)
            dA_k = griddata(sampleIndices, dA_ps[i,:], tuple(gridDomain), method=self.kind)
            A_pk.append(A_k)
            dA_pk.append(dA_k)
        A_pk = np.array(A_pk)
        dA_pk = np.array(dA_pk)
        A_pk = A_pk.reshape((numberProps, A_k.size))
        dA_pk = dA_pk.reshape((numberProps, dA_k.size))        
        self.EA_pk = A_pk
        self.dEA_pk = dA_pk
        return (A_pk, dA_pk)
        
    def computeExpectations(self, A_psn, I_s, gridShape):
        """
        Parameters:
        -----------
        A_psn :    bi-dimensional list P x S of (np.ndarray, float, shape N(S))
                      P - number of properties
                      S - number of sampled states
                      N - number of uncorrelated configurations
                      A_psn[p,s,n] - property 'p' for configuration 'n' of sampled state 's'

        I_s   :    np.ndarray of shape (S)
                      I_s[s] - linear index of sampled state 's'
        
        gridShape : tuple (n_1, n_2, ..., n_D), int
                      n_d - number of partitions for edge 'd'

        Returns:
        --------
        (EA_pk, dEA_pk) :    EA_pk :    np.ndarray, float, shape (P,K)
                            dEA_pk :    np.ndarray, float, shape (P,K)

        """
        # mean and std for sampled states
        (A_ps, dA_ps) = self._computeAvgStd(A_psn)
        return self.computeExpectationsFromAvgStd(A_ps, dA_ps, I_s, gridShape)


class GaussianProcessRegressionInterpolation (Interpolation):
    # GPR inherits from Interpolation but overrides the actual
    # interpolation operation.
    reweight = False
    corners = True
    kind = 'gpr'

    def __init__(self): 
        return

    @staticmethod
    def _featureNormalization(Xsamples, Xgrid):
        cp = Xsamples.copy().astype('float64')
        for i in range(Xsamples.ndim):
            max = np.max(Xgrid[:, i])
            min = np.min(Xgrid[:, i])
            cp[:, i] = (Xsamples[:, i] - min)/(max - min)
        return cp

    @staticmethod
    def _GPfit(Xarray, Yarray, Xpred, coreKernel=None, noiseLevel=None, optimizationCycles=20):
        # Normalize each feature to [0,1] interval
        X_ = GaussianProcessRegressionInterpolation._featureNormalization(Xarray, Xpred)
        Xpred_ = GaussianProcessRegressionInterpolation._featureNormalization(Xpred, Xpred)
        # If there is only one feature, reshape the arrays Xarray and Xpred
        if (Xarray.ndim == 1):
            X_ = Xarray.reshape(-1, 1)
        if (Xpred.ndim == 1):
            Xpred_ = Xpred.reshape(-1, 1)
        if (coreKernel is None):
            # Loop over test core Kernels
            kernels = [ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * RBF(),
                       ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * Matern(),
                       ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * DotProduct(),
                       ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * ExpSineSquared(),
                       ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * RationalQuadratic()]
        else:
            # Use given kernel only
            kernels = [coreKernel,]
        if (noiseLevel is not None):
            # Add noise level
            for i, k in enumerate(kernels):
                kernels[i] = k + WhiteKernel(noise_level=noiseLevel, noise_level_bounds='fixed')
        maxLogLikelihood = -1.0e+23
        optGP = None
        # Loop over kernels
        for k in kernels:
            # Initialize GPR process
            gp = gaussian_process.GaussianProcessRegressor(kernel=k,
                                                           n_restarts_optimizer=optimizationCycles)
            # Fit model
            gp.fit(X_, Yarray)
            # Test log of marginal likelihood
            if (gp.log_marginal_likelihood_value_ > maxLogLikelihood):
                optGP = gp
                maxLogLikelihood = gp.log_marginal_likelihood_value_
        print('Optimal GPR is ', optGP)
        print('LogLikelihood = ', maxLogLikelihood)
        return optGP.predict(Xpred_, return_std=True) # ypred, sigma

    def computeExpectationsFromAvgStd(self, A_ps, dA_ps, I_s, gridShape):
        # make grid from shape
        gridDomain = np.meshgrid(*[np.arange(s) for s in gridShape], indexing='ij')
        # determine grid indices for sampled states
        sampleIndices = []
        gridIndices = []
        linearIdx = 0
        for idx in np.ndindex(gridShape):
            gridIndices.append(np.array(idx))
            if linearIdx in I_s:
                sampleIndices.append(np.array(idx))
            linearIdx += 1
        gridIndices = np.array(gridIndices)
        sampleIndices = np.array(sampleIndices)
        # interpolate property data
        A_pk = []
        dA_pk = []
        numberProps = A_ps.shape[0]
        for i in range(numberProps):
            noise_level = np.mean(dA_ps[i,:]**2)
            A_k, dA_k  = GaussianProcessRegressionInterpolation._GPfit(sampleIndices, A_ps[i,:], gridIndices, noiseLevel=noise_level)
            A_pk.append(A_k)
            dA_pk.append(dA_k)
        A_pk = np.array(A_pk)
        dA_pk = np.array(dA_pk)
        A_pk = A_pk.reshape((numberProps, A_k.size))
        dA_pk = dA_pk.reshape((numberProps, dA_k.size))        
        self.EA_pk = A_pk
        self.dEA_pk = dA_pk
        return (A_pk, dA_pk)
        
    def computeExpectations(self, A_psn, I_s, gridShape):
        """
        Parameters:
        -----------
        A_psn :    bi-dimensional list P x S of (np.ndarray, float, shape N(S))
                      P - number of properties
                      S - number of sampled states
                      N - number of uncorrelated configurations
                      A_psn[p,s,n] - property 'p' for configuration 'n' of sampled state 's'

        I_s   :    np.ndarray of shape (S)
                      I_s[s] - linear index of sampled state 's'
        
        gridShape : tuple (n_1, n_2, ..., n_D), int
                      n_d - number of partitions for edge 'd'

        Returns:
        --------
        (EA_pk, dEA_pk) :    EA_pk :    np.ndarray, float, shape (P,K)
                            dEA_pk :    np.ndarray, float, shape (P,K)

        """
        # mean and std for sampled states
        (A_ps, dA_ps) = self._computeAvgStd(A_psn)
        return self.computeExpectationsFromAvgStd(A_ps, dA_ps, I_s, gridShape)
