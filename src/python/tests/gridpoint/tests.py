import sys
sys.path.append('../..')

import os
import unittest
import gridoptimizer
from   RunFromInput import GridPoint, PropertyBase, Density, LiquidProtocol
from   copy import deepcopy
import numpy as np

dirname = os.path.dirname(__file__)

point = GridPoint("", 32)
propexample = Density(550.00, 0.2)
point.add_property_estimate('id-of-prop', 'density', propexample)
point.add_itp("/some/example/of/path")

prot = LiquidProtocol('my-liquid', 1500, '/path/to/coords', [], [2.3,2.3,2.3], [])
other_prot = LiquidProtocol('other-liquid', 1500, '/path/to/coords', [], [2.3,2.3,2.3], [])

class TestGP(unittest.TestCase):

    def test_internal_struct(self):
        self.assertEqual(point.itp_path, "/some/example/of/path")
        self.assertEqual(point.get_property_estimate('id-of-prop'), 550.00)
        self.assertEqual(point.get_property_err('id-of-prop'), 0.2)
        self.assertFalse(point.is_sample)
        point.set_as_sample()
        self.assertTrue(point.is_sample)

    def test_protocol_outputs(self):
        output_dict = {'xtc': '/path/to/xtc', 'tpr': '/path/to/tpr'}
        point.add_protocol_output(prot, output_dict)
        expected_dict = {'my-liquid': {'xtc': '/path/to/xtc', 'tpr': '/path/to/tpr'}}
        self.assertEqual(expected_dict, point.protocol_outputs)
        point.add_protocol_output(other_prot, output_dict)
        expected_dict = {'my-liquid': {'xtc': '/path/to/xtc', 'tpr': '/path/to/tpr'},
                         'other-liquid': {'xtc': '/path/to/xtc', 'tpr': '/path/to/tpr'}}
        self.assertEqual(expected_dict, point.protocol_outputs)

    def test_add_properties(self):
        output_prop = '/path/to/property.xvg'
        point.add_atomic_property_output('dummy-atom-prop', prot, output_prop)
        expected_dict = {'my-liquid': {'dummy-atom-prop': '/path/to/property.xvg'}}
        self.assertEqual(expected_dict, point.atomic_properties)
        output_prop = '/path/to/other/property.xvg'
        point.add_atomic_property_output('other-dummy-atom-prop', prot, output_prop)
        expected_dict = {'my-liquid': {'dummy-atom-prop': '/path/to/property.xvg',
                                       'other-dummy-atom-prop': '/path/to/other/property.xvg' }}
        self.assertEqual(expected_dict, point.atomic_properties)
        point.add_atomic_property_output('dummy-atom-prop', other_prot, output_prop)
        expected_dict = {'my-liquid': {'dummy-atom-prop': '/path/to/property.xvg',
                                       'other-dummy-atom-prop': '/path/to/other/property.xvg' },
                         'other-liquid': {'dummy-atom-prop': '/path/to/other/property.xvg'}}
        self.assertEqual(expected_dict, point.atomic_properties)
        return

if __name__ == '__main__':
    unittest.main()

