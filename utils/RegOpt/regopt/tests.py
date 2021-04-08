#!/usr/bin/env python3

from normalization import DefaultNormalizationTest
from regression import GaussianProcessRegressorTest
from optimizer import OptimizerSimplexTest
from regopt import OptimizerDriverTest

if __name__ == '__main__':
    #DefaultNormalizationTest().runTests()
    #GaussianProcessRegressorTest().runTests()
    #OptimizerSimplexTest().runTests()
    OptimizerDriverTest().runTests()
