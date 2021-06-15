import numpy as np
import sklearn
import sys
from . import normalization
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split

from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel, Matern, DotProduct, ExpSineSquared, RationalQuadratic

def plot_train_test_data(fn, Xtrain, Xtest, Ypred_train, Ystd_train, Ytrain, Ypred_test, Ystd_test, Ytest):
    
    fig = plt.figure()
    ax  = plt.gca()

    ax.errorbar(Ytrain, Ypred_train, yerr=Ystd_train, label='Train', fmt='o')

    ax.errorbar(Ytest,  Ypred_test,  yerr=Ystd_test,  label='Test', fmt='o')

    Yfull = np.append(Ytrain, Ytest, axis=0)

    limits = [np.min(Yfull), np.max(Yfull)]

    ax.plot(limits, limits, 'k--')

    plt.legend()

    plt.savefig(fn)

def save_train_test_data(fn, Xtrain, Xtest, Ypred_train, Ystd_train, Ytrain, Ypred_test, Ystd_test, Ytest):

    fp = open(fn, 'w')

    fp.write("# Train\n")
    for yt, yp, yerr in zip(Ytrain, Ypred_train, Ystd_train):
        fp.write("%14.3f%14.3f%14.3f\n" % (yt, yp, yerr))

    fp.write("# Test\n")
    for yt, yp, yerr in zip(Ytest, Ypred_test, Ystd_test):
        fp.write("%14.3f%14.3f%14.3f\n" % (yt, yp, yerr))

    fp.close()

class RegressorBase:

    def __init__(self, normalization=None, train_test_split_funct=None):
        if (normalization is None):
            self.normalization = EmptyNormalization()
        else:
            self.normalization = normalization

        if (train_test_split_funct is None):
            self.train_test_split_funct = train_test_split
        else:
            self.train_test_split_funct = train_test_split_funct

    # Virtual Functions
    # -----------------

    def train(self, x, y):
        pass

    # Virtual. Assumes normalization.
    def estimate(self, x):
        pass

    def printInfo(self, stream=sys.stdout):
        pass

    def getQuality(self):
        pass

    # Interface
    # ---------

    def predict(self, x, normalize=True):
        if (normalize):
            xn_ = self.normalization.normalize(x)
            return self.estimate(xn_)
        else:
            return self.estimate(x)

    def _fit(self, x, y):
        if (self.normalization is None):
            xn_ = x
        else:
            xn_ = self.normalization.normalize(x)
        Xtrain_, Xtest_, Ytrain_, Ytest_ = self.train_test_split_funct(xn_, y)

        self.train(Xtrain_, Ytrain_)

        Ypred_train_, Ystd_train_ = self.predict(Xtrain_, normalize=False)
        Ypred_test_, Ystd_test_ = self.predict(Xtest_, normalize=False)

        self.train_rmse_ = np.sqrt(np.mean((Ypred_train_ - Ytrain_)**2))
        self.test_rmse_ = np.sqrt(np.mean((Ypred_test_ - Ytest_)**2))

        # compute hits in train
        self.train_hits_percent_sigma_ = 0.0
        self.train_hits_percent_2sigma_  = 0.0
        for y, yest, ystd in zip(Ytrain_, Ypred_train_, Ystd_train_):
            if (ystd == 0.0):
                self.train_hits_percent_sigma_  += 1.0
                self.train_hits_percent_2sigma_ += 1.0
            elif abs(y - yest) <= ystd:
                self.train_hits_percent_sigma_  += 1.0
                self.train_hits_percent_2sigma_ += 1.0
            elif abs(y - yest) <= 2*ystd:
                self.train_hits_percent_2sigma_ += 1.0
        self.train_hits_percent_sigma_ /= (1.0 * Xtrain_.shape[0])
        self.train_hits_percent_2sigma_ /= (1.0 * Xtrain_.shape[0])

        # compute hits in test
        self.test_hits_percent_sigma_ = 0.0
        self.test_hits_percent_2sigma_  = 0.0
        for y, yest, ystd in zip(Ytest_, Ypred_test_, Ystd_test_):
            if (ystd == 0.0):
                self.test_hits_percent_sigma_  += 1.0
                self.test_hits_percent_2sigma_ += 1.0
            elif abs(y - yest) <= ystd:
                self.test_hits_percent_sigma_  += 1.0
                self.test_hits_percent_2sigma_ += 1.0
            elif abs(y - yest) <= 2*ystd:
                self.test_hits_percent_2sigma_ += 1.0
        self.test_hits_percent_sigma_ /= (1.0 * Xtest_.shape[0])
        self.test_hits_percent_2sigma_ /= (1.0 * Xtest_.shape[0])

        return Xtrain_, Xtest_, Ypred_train_, Ystd_train_, Ytrain_, Ypred_test_, Ystd_test_, Ytest_

    def fit(self, x, y, plot=None, data=None):

        fit_results = self._fit(x, y)

        if (plot is not None):
            plot_train_test_data(plot, *fit_results)

        if (data is not None):
            save_train_test_data(data, *fit_results)

class GaussianProcessRegressor(RegressorBase):

    def __init__(self, normalization=None, train_test_split_funct=None, noise_level=0):
        super().__init__(normalization, train_test_split_funct)
        self.noise_kernel = WhiteKernel(noise_level=noise_level, noise_level_bounds='fixed')

    def train(self, x, y):
        # test core Kernels
        kernels = [ConstantKernel() * RBF(),
                ConstantKernel() * Matern(),
                ConstantKernel() * DotProduct(),
                ConstantKernel() * ExpSineSquared(),
                ConstantKernel() * RationalQuadratic()]
        for i, k in enumerate(kernels):
            kernels[i] = k + self.noise_kernel

        maxLogLikelihood = -1.0e+23
        optGP = None

        # Loop over kernels
        for k in kernels:
            # Initialize GPR process
            gp = sklearn.gaussian_process.GaussianProcessRegressor(kernel=k,
                                                                   n_restarts_optimizer=5)
            # Fit model
            gp.fit(x, y)
            # Test log of marginal likelihood
            if (gp.log_marginal_likelihood_value_ > maxLogLikelihood):
                optGP = gp
                maxLogLikelihood = gp.log_marginal_likelihood_value_

        self.log_marginal_likelihood_ = optGP.log_marginal_likelihood_value_
        self.model = optGP


    def estimate(self, x):
        return self.model.predict(x, return_std=True)

    def getQuality(self):
        return self.log_marginal_likelihood_

    def printResults(self, stream=sys.stdout):
        stream.write(self.model.__repr__())
        stream.write("\n")
        stream.write("Log marginal likelihood: {:<.4f}".format(self.log_marginal_likelihood_))
        stream.write("\n")

# Tests

class GaussianProcessRegressorTest:

    def _fitCore(self, X, Y, normalization, noiseLevel=0):

        GPR = GaussianProcessRegressor(normalization=normalization, noise_level=noiseLevel)
        GPR.fit(X, Y)

        # Assert the model performs OK for this simple data.
        assert(GPR.train_rmse_ < 1.0)
        assert(GPR.test_rmse_ < 1.0)
        assert(GPR.train_hits_percent_sigma_ == 1.0)

    def fit1D(self, noiseLevel=0):
        X = np.linspace(0, 1, num=10).reshape(-1, 1)
        Y = 1.15 * np.linspace(0, 1, num=10)
        norm = normalization.DefaultNormalization(np.array([0.0]), np.array([1.00]))
        self._fitCore(X, Y, norm, noiseLevel=noiseLevel)

    def fit1DNoise(self):
        self.fit1D(noiseLevel=0.1)

    def fitMultiD(self):
        X = np.array([[0.00, 1.00],
                      [0.25, 0.75],
                      [0.50, 0.50],
                      [0.75, 0.25],
                      [1.00, 0.00]])
        Y = np.array([1.00, 0.75, 0.50, 0.75, 0.00])
        norm = normalization.DefaultNormalization(np.array([0.0, 0.0]), np.array([1.00, 1.00]))
        self._fitCore(X, Y, norm)

    def runTests(self):
        self.fit1D()
        self.fit1DNoise()
        self.fitMultiD()
