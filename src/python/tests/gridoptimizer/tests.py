import sys
sys.path.append('../..')

import os
import unittest
import gridoptimizer
from   RunFromInput import ParameterGrid, GridPoint, Density
from   copy import deepcopy
import numpy as np

dirname = os.path.dirname(__file__)

# optimizer attributes
opt = gridoptimizer.gridOptimizer(maxSteps=0, percentCutoff=.45)
opt.referenceValues = {'first': 515.00, 'second': 350.00}
opt.referenceTolerances = {'first': 0.4, 'second': 0.8}
opt.referenceWeights = {'first': 0.5, 'second': 1.0}

# test grid
dummy_grid = ParameterGrid()
dummy_grid.dim = 1
dummy_grid.size = [5,]
dummy_grid.linear_size = 5
dummy_grid.grid_points = [
    GridPoint("", 0),
    GridPoint("", 1),
    GridPoint("", 2),
    GridPoint("", 3),
    GridPoint("", 4)]
# add some property estimates -- e.g. of Density
# for property id 'first'
dummy_grid.grid_points[0].estimated_properties['first'] = Density(500.00, 0.2)
dummy_grid.grid_points[1].estimated_properties['first'] = Density(510.00, 0.4)
dummy_grid.grid_points[2].estimated_properties['first'] = Density(520.00, 0.3)
dummy_grid.grid_points[3].estimated_properties['first'] = Density(530.00, 0.2)
dummy_grid.grid_points[4].estimated_properties['first'] = Density(540.00, 0.4)
# for property id 'second'
dummy_grid.grid_points[0].estimated_properties['second'] = Density(400.00, 0.5)
dummy_grid.grid_points[1].estimated_properties['second'] = Density(410.00, 1.4)
dummy_grid.grid_points[2].estimated_properties['second'] = Density(320.00, 0.5)
dummy_grid.grid_points[3].estimated_properties['second'] = Density(630.00, 0.3)
dummy_grid.grid_points[4].estimated_properties['second'] = Density(240.00, 0.8)

expected_ranked_scores = [
    (2, 0.070209536726226),
    (0, 0.117848296337329),
    (1, 0.140083034655634),
    (4, 0.258139180921841),
    (3, 0.653413686296367)]

class TestOptimizer(unittest.TestCase):

    def test_fillWithScores(self):
        my_opt = deepcopy(opt)
        my_grid = deepcopy(dummy_grid)
        my_opt.fillWithScores(my_grid)
        self.assertListEqual([x[0] for x in expected_ranked_scores], [x[0] for x in my_opt.stateScores])
        self.assertTrue(np.allclose(np.array([x[1] for x in expected_ranked_scores]),
                                    np.array([x[1] for x in my_opt.stateScores])))

    def test_determineNextSample(self):
        my_opt = deepcopy(opt)
        my_grid = deepcopy(dummy_grid)
        my_opt.fillWithScores(my_grid)
        self.assertEqual(my_opt.determineNextSample(my_grid), 1)

if __name__ == '__main__':
    unittest.main()

