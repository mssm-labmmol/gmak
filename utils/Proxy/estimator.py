from abc import ABC, abstractmethod
from grid import Grid
import numpy as np


class BaseEstimator(ABC):
    @abstractmethod
    def get_estimate(self, xs):
        pass

class SampleBasedMixin(ABC):
    @abstractmethod
    def replace_samples(self, new_samples):
        pass


class FunctionalEstimator(BaseEstimator):
    def __init__(self, func):
        """
        Parameters:

        func: xs -> double
        """
        self.func = func

    def get_estimate(self, xs):
        return self.func(xs)

    def get_err_estimate(self, xs):
        return 0


class NoisyEstimator(BaseEstimator):
    def __init__(self, func_avg, func_std, nsamples):
        """
        Parameters:

        func_avg: xs -> double parameterizes avg
        func_std: xs -> double parameterizes std dev
        nsamples: int, number of uncorrelated samples
        """
        self.func_avg = func_avg
        self.func_std = func_std
        self.nsamples = nsamples
        self.ncalls   = 0

    def get_estimate(self, xs):
        self.ncalls += 1
        out = np.random.normal(loc=self.func_avg(xs),
                                scale=self.func_std(xs))
        return out

    def get_err_estimate(self, xs):
        print("WARNING: My get_err_estimate function is WRONG!"
              " Rethink it with (n-1)s^2/sigma^2 ~ chisq_{n-1}")
        n = self.nsamples
        std_factor = np.sqrt(1 - (1.0/np.e) * (1 + 1.0/n)**n)
        std = self.func_std(xs)
        return np.random.normal(loc=std, scale=std_factor*std)
        

class Interpolator(BaseEstimator, SampleBasedMixin):
    def __init__(self, grid, samples, sampler, method):
        self.sampler = sampler
        self.samples = samples
        self.method  = method # linear or cubic
        self.grid = grid
        self.up_to_date = False

    def update_model(self):
        from scipy.interpolate import griddata
        if not self.up_to_date:
            # get sampled values from sampler.
            values = [self.sampler.get_estimate(xs) for xs in self.samples]
            self.interp = griddata(self.samples,
                                   values,
                                   self.grid.to_griddata(),
                                   method=self.method)
        self.up_to_date = True

    def replace_samples(self, new_samples):
        """Parameters:

        new_samples: list of x points for new_samples
        """
        self.samples = new_samples
        self.up_to_date = False

    def _get_estimate_idx(self, idxs):
        return self.interp[tuple(idxs)]

    def get_estimate(self, xs):
        # locate `xs` in grid and get estimate.
        self.update_model()
        loc = self.grid.locate(xs)
        out = self._get_estimate_idx(loc)
        return out


class GPR(BaseEstimator, SampleBasedMixin):
    def __init__(self, samples, sampler, normalizer):
        self.sampler = sampler
        self.samples = samples
        self.normalizer = normalizer
        self.fit()

    def replace_samples(self, new_samples):
        """Parameters:

        new_samples: list of x points for new_samples
        """
        self.samples = new_samples
        self.fit()
        
    def fit(self):
        from sklearn import gaussian_process
        from sklearn.gaussian_process.kernels import (RBF,
                                                      ConstantKernel,
                                                      WhiteKernel,
                                                      Matern,
                                                      DotProduct,
                                                      ExpSineSquared,
                                                      RationalQuadratic)
        # get sampled values
        values_avg = np.array([self.sampler.get_estimate(xs)
                               for xs in self.samples])
        values_err = np.array([self.sampler.get_err_estimate(xs)
                               for xs in self.samples])
        Xs = np.array([self.normalizer(xs) for xs in self.samples])
        # reshape if feature is 1D
        if Xs.ndim == 1:
            Xs = Xs.reshape(-1, 1)
        # Kernels to be tested.
        kernels = [ConstantKernel() * RBF(),
                ConstantKernel() * Matern(),
                ConstantKernel() * DotProduct(),
                ConstantKernel() * ExpSineSquared(),
                ConstantKernel() * RationalQuadratic()]
        # Add noise level
        for i, k in enumerate(kernels):
            kernels[i] = k + WhiteKernel(noise_level=np.mean(values_err),
                                         noise_level_bounds='fixed')

        maxLogLikelihood = -1.0e+23
        optGP = None
        # Loop over kernels
        for k in kernels:
            # Initialize GPR process
            gp = gaussian_process.GaussianProcessRegressor(kernel=k,
                                                           n_restarts_optimizer=20)
            # Fit model
            gp.fit(Xs, values_avg)
            # Test log of marginal likelihood
            if (gp.log_marginal_likelihood_value_ > maxLogLikelihood):
                optGP = gp
                maxLogLikelihood = gp.log_marginal_likelihood_value_
        self.gpr = optGP

    def get_estimate(self, xs):
        xn = self.normalizer(xs)
        out = self.gpr.predict([xn])[0]
        return out

    def get_err_estimate(self, xs):
        xn = self.normalizer(xs)
        return self.gpr.predict([xn], return_std=True)[1][0]
