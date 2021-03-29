import numpy as np

# Extension: A new Normalization class should contain the method
# `normalize', which receives an input values `x' and returns a
# normalized copy of `x'.

class DefaultNormalization:

    def __init__(self, inf_lim, sup_lim):
        self.inf_lim = inf_lim
        self.sup_lim = sup_lim

    def normalize(self, x):
        _x = x.copy()

        for i in range(_x.shape[1]):
            max = self.sup_lim[i]
            min = self.inf_lim[i]

            _x[:,i] = (_x[:,i] - min)/(max - min)
        
        return _x

class EmptyNormalization:

    def __init__(self):
        return

    def normalize(self, x):
        return x.copy()

# Unit Testing

class DefaultNormalizationTest:

    def __init__(self):
        return

    def preserveNormalizationOnNormalized(self):
        X = np.array([[0.0, 0.2, 0.3],
                      [0.4, 1.0, 0.1],
                      [0.1, 0.1, 1.0],
                      [1.0, 0.0, 0.0]])
        norm = DefaultNormalization(np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 1.0]))
        X_n  = norm.normalize(X) 
        np.testing.assert_array_equal(X, X_n)

    def checkSingleDimension(self):
        X = np.array([[0.0],
                      [0.4],
                      [0.1],
                      [1.0]])
        norm = DefaultNormalization(np.array([0.0]), np.array([1.0]))
        X_n  = norm.normalize(X)
        np.testing.assert_array_equal(X, X_n)

    def normalizeLinear(self):
        X_base = np.array([[0.0, 0.2, 0.3],
                           [0.4, 1.0, 0.1],
                           [0.1, 0.1, 1.0],
                           [1.0, 0.0, 0.0]])
        X = 10.0 * X_base
        norm = DefaultNormalization(np.array([0.0, 0.0, 0.0]), np.array([10.0, 10.0, 10.0]))
        X = norm.normalize(X)
        np.testing.assert_array_equal(X_base, X)

    def checkNormalizationOnBiggerRange(self):
        X_base = np.array([[0.0, 0.2, 0.3],
                           [0.4, 1.0, 0.1],
                           [0.1, 0.1, 1.0],
                           [1.0, 0.0, 0.0]])
        Xrange = 10.0 * X_base
        X_test = np.array([[5.0, 4.0, 3.0]])
        norm = DefaultNormalization(np.array([0.0, 0.0, 0.0]), np.array([10.0, 10.0, 10.0]))
        X = norm.normalize(X_test)
        np.testing.assert_equal(X, np.array([[0.5, 0.4, 0.3]]))

    def runTests(self):
        self.preserveNormalizationOnNormalized()
        self.checkSingleDimension()
        self.normalizeLinear()
        self.checkNormalizationOnBiggerRange()
