import sys
sys.path.append('../..')

import os
import unittest
import surrogate_model 
import numpy as np

dirname = os.path.dirname(__file__)

class TestMBAR(unittest.TestCase):
    
    def compute_and_write(self, nprops, bool_legacy):
        nstates = 1089
        # mount A_pkn from data
        A_pkn = [None] * nprops
        for i in range(nprops):
            prop_file = os.path.join(dirname, 'data/prop_{}.dat'.format(i))
            A_pkn[i] = np.loadtxt(prop_file)
        A_pkn = np.array(A_pkn)
        # mount u_kn from data/u_kn.dat
        u_kn = np.loadtxt(os.path.join(dirname, 'data/u_kn.dat'))
        # mount N_k from data/N_k.dat
        N_k = np.loadtxt(os.path.join(dirname, 'data/N_k.dat'))
        # compute
        if (bool_legacy):
            mbar = surrogate_model.MBARLegacy()
        else:
            mbar = surrogate_model.MBAR()
        mbar.computeExpectations(A_pkn, u_kn, N_k)
        # write log files
        if (bool_legacy):
            output_dir = os.path.join(dirname, 'output_mbar_legacy_{}'.format(nprops))
        else:
            output_dir = os.path.join(dirname, 'output_mbar_{}'.format(nprops))
        ref_dir = os.path.join(dirname, 'reference')
        mbar.writeLogToDirectory(output_dir)
        for i in range(nprops):
            # write output data
            mbar.writeExpectationsToFile(output_dir + '/prop_{}_EA_k.dat'.format(i),
                    output_dir + '/prop_{}_dEA_k.dat'.format(i), i)
            # compare output data with reference based on file values
            # average
            this_val = np.loadtxt(output_dir + '/prop_{}_EA_k.dat'.format(i))
            ref_val  = np.loadtxt(ref_dir + '/prop_{}_EA_k.dat'.format(i))
            self.assertTrue(np.allclose(this_val, ref_val))
            # std
            this_val = np.loadtxt(output_dir + '/prop_{}_dEA_k.dat'.format(i))
            ref_val  = np.loadtxt(ref_dir + '/prop_{}_dEA_k.dat'.format(i))
            self.assertTrue(np.allclose(this_val, ref_val))
        # compare log data with reference
        self.assertTrue(np.array_equal(np.loadtxt(output_dir + '/N_k.dat'),
                np.loadtxt(ref_dir + '/N_k.dat')))
        self.assertTrue(np.allclose(np.loadtxt(output_dir + '/Eff_k.dat'),
                np.loadtxt(ref_dir + '/Eff_k.dat')))
        self.assertTrue(np.allclose(np.loadtxt(output_dir + '/W_k.dat'),
                np.loadtxt(ref_dir + '/W_k.dat')))
        # remove created files
        os.system('rm -rf ' + output_dir)
   
    def test_nonlegacy_compute_and_write_single_property(self):
        self.compute_and_write(1, False)

    def test_nonlegacy_compute_and_write_multiple_properties(self):
        self.compute_and_write(2, False)

    def test_legacy_compute_and_write_single_property(self):
        self.compute_and_write(1, True)

    def test_legacy_compute_and_write_multiple_properties(self):
        self.compute_and_write(2, True)

#class TestInterpolate(unittest.TestCase):
#
#    def test_nearest(self):
#        return
#
#    def test_linear(self):
#        return
#
#    def test_cubic(self):
#        return

if __name__ == '__main__':
    unittest.main()

