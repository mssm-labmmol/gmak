import numpy as np
import random
import pandas as pd
import sys
from   . import normalization
from   . import regression
from   . import optimizer
from   functools import partial
from sklearn.model_selection import train_test_split

class IOptimizationDriverNormalizationFactory:

    # virtual
    def getNormalization(self):
        pass

class OptimizationDriverDefaultNormalizationFactory(IOptimizationDriverNormalizationFactory):

    def __init__(self, optDriver):
        self.driver = optDriver

    def getNormalization(self, output_name):
        x = self.driver.getInputData(output_name)
        norms = {}
        inf_lim = np.zeros((x.shape[1]))
        sup_lim = np.zeros((x.shape[1]))
        for i in range(x.shape[1]):
            inf_lim[i] = np.min(x[:,i])
            sup_lim[i] = np.max(x[:,i])
        norm = normalization.DefaultNormalization(inf_lim=inf_lim, sup_lim=sup_lim)
        return norm

class OptimizationDriver:

    def __init__(self, test_size=0.25, prefix=None):
        # TODO: Specify factories for normalization, regressor, 
        # optimizer and function based on input strings.
        self._data = dict()
        self.test_size = test_size
        self.prefix = prefix
        return

    def minimFunc(self, x):
        chisq = 0.0
        wei_sum = 0.0
        prediction = self.predict(x)
        for k in prediction.keys():
            this_chisq = self._data[k]['weight'] * ((prediction[k][0][0] - self._data[k]['reference_value'])/self._data[k]['reference_value'])**2
            chisq += this_chisq
            wei_sum += self._data[k]['weight']
        chisq /= wei_sum
        return np.sqrt(chisq)
        
    def predict(self, x):
        return {k: self._data[k]['regressor'].predict(np.array([x,])) for k in self._data.keys()}

    def getInputData(self, output_name):
        return self._data[output_name]['input_values']

    def addPropertyFromCSV(self, csv_fn, weight, referenceValue, sep='\s+'):
        _df = pd.read_csv(csv_fn, sep=sep)
        _cols = list(_df.columns.values)
        self._data[_cols[-2]] = {
            'input_names': _cols[:-2],
            'input_values': _df[_cols[:-2]].to_numpy(),
            'output_values': _df[_cols[-2]].to_numpy(),
            'output_error_values': _df[_cols[-1]].to_numpy(),
            'reference_value': referenceValue,
            'weight': weight,
        }
        
        norm = OptimizationDriverDefaultNormalizationFactory(self).getNormalization(_cols[-2])

        self._data[_cols[-2]]['regressor'] = regression.GaussianProcessRegressor(
            normalization=norm,
            train_test_split_funct=partial(train_test_split, test_size=self.test_size),
            noise_level=np.mean(self._data[_cols[-2]]['output_error_values']),
        )

        self.optimizer = optimizer.OptimizerSimplex()

    def fit(self, format='pdf'):
        for k in self._data.keys():
            print(f"Fitting model for property `{k}'...")
            if (self.prefix is not None):
                self._data[k]['regressor'].fit(self._data[k]['input_values'], self._data[k]['output_values'], plot=(self.prefix + "_" + k + "." + format), data=(self.prefix + "_" + k + ".dat"))
            else:
                self._data[k]['regressor'].fit(self._data[k]['input_values'], self._data[k]['output_values'])
            print("Done.")

    def optimize(self):
        print("Optimizing...")
        name = list(self._data.keys())[0]
        out = self.optimizer.optimize(self.minimFunc, x_0=random.choice(self._data[name]['input_values']))
        print("Done.")
        return out

    def getAverageQuality(self):
        outs = []
        for k in self._data.keys():
            outs.append( self._data[k]['regressor'].getQuality() )
        return np.mean(outs)

    def printResults(self, stream=sys.stdout):
        """
        Prints the following to `stream':
            - optimal parameter set
            - expected values of the properties (and errors) for the optimal parameter set
            - information about the regression model, which depends on the model used
            - average quality of the models
        """
        #self.optimizer.printResults(stream)

        #for prop in self._data.keys():
        #    self._data[prop]['regressor'].printResults(stream)

        stream.write("----------------------------------------------------------------------\n")

        stream.write("Average quality metric: {}\n".format(self.getAverageQuality()))

        stream.write("----------------------------------------------------------------------\n")

        # Expected values and errors for each property
        x, y = self.optimizer.getResults()

        stream.write("%-14s = " % ("X"))
        for xi in x:
            stream.write("%14.3e" % xi)
        stream.write("\n")

        prediction = self.predict(x)

        for prop in prediction.keys():
            stream.write("%-14s = %14.3e +/- %-14.3e\n" % (prop,
                                                          prediction[prop][0][0],
                                                          prediction[prop][1][0]))

        stream.write("%-14s = %14.3e\n" % ("Score", y))
        

class OptimizerDriverTest:

    def runTests(self):
        X = np.array([
                [0.0, 0.0],
                [0.0, 0.25],
                [0.0, 0.50],
                [0.0, 0.75],
                [0.0, 1.00],
                [0.25, 0.0],
                [0.25, 0.25],
                [0.25, 0.50],
                [0.25, 0.75],
                [0.25, 1.00],
                [0.50, 0.0],
                [0.50, 0.25],
                [0.50, 0.50],
                [0.50, 0.75],
                [0.50, 1.00],
                [0.75, 0.0],
                [0.75, 0.25],
                [0.75, 0.50],
                [0.75, 0.75],
                [0.75, 1.00],
                [1.00, 0.0],
                [1.00, 0.25],
                [1.00, 0.50],
                [1.00, 0.75],
                [1.00, 1.00],
            ])
        Y = 0.5 * (X[:,0] + X[:,1])
        driver = OptimizationDriver()
        driver._data = {
            'A': {
                'input_values': X.copy(),
                'output_values': 2*Y,
                'reference_value': 1.0,
                'weight': 1.0,
                'regressor': regression.GaussianProcessRegressor(
                    normalization=normalization.DefaultNormalization(inf_lim=np.array([0,0]), sup_lim=np.array([1.0, 1.0])),
                    noise_level=0.01,
                )
            },
            'B': {
                'input_values': X.copy(),
                'output_values': 1.2*Y,
                'reference_value': 0.6,
                'weight': 1.0,
                'regressor': regression.GaussianProcessRegressor(
                    normalization=normalization.DefaultNormalization(inf_lim=np.array([0,0]), sup_lim=np.array([1.0, 1.0])),
                    noise_level=0.02,
                )
            },
            'C': {
                'input_values': X.copy(),
                'output_values': 1.5*Y,
                'reference_value': 0.75,
                'weight': 1.0,
                'regressor': regression.GaussianProcessRegressor(
                    normalization=normalization.DefaultNormalization(inf_lim=np.array([0,0]), sup_lim=np.array([1.0, 1.0])),
                    noise_level=0.01,
                )
            }
        }
        driver.optimizer = optimizer.OptimizerSimplex()
        driver.fit()
        driver.optimize()
