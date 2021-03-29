import numpy as np
import sys

from scipy.optimize import minimize

class OptimizerSimplex:
    def optimize(self, func, x_0):
        """
        Parameters
        ----------
        func : callable
            Function to be minimized, called as func(x), where x is 
            a numpy.ndarray.
        x_0 : numpy.ndarray
            Initial guess for minimization.
        
        Returns
        -------
        (x, y) : tuple
            Tuple with results of minimization.
        """
        self.results = minimize(func, x_0, method='Nelder-Mead')
        return self.results.x, self.results.fun

    def printResults(self, stream=sys.stdout):
        stream.write(self.results.__repr__())
        stream.write('\n')

    def getResults(self):
        return self.results.x, self.results.fun

class OptimizerSimplexTest:

    def quadratic_fun(self, x):
        out = 0.0
        for k in x:
            out += (k - 1.0)**2
        return out

    def testQuadratic3D(self):
        opt = OptimizerSimplex()
        x_0 = np.array([0.1, 1.4, 4.0])
        x_min, y_min = opt.optimize(self.quadratic_fun, x_0)
        np.testing.assert_allclose(x_min, np.array([1.0, 1.0, 1.0]), rtol=1e-4)

    def testQuadratic1D(self):
        opt = OptimizerSimplex()
        x_0 = np.array([4.0])
        x_min, y_min = opt.optimize(self.quadratic_fun, x_0)
        np.testing.assert_allclose(x_min, np.array([1.0]), rtol=1e-4)

    def runTests(self):
        self.testQuadratic1D()
        self.testQuadratic3D()
