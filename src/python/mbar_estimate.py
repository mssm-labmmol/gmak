#!/usr/bin/python
#

import numpy as np
import pymbar 
import os

def prepare_mbar_no_pv (u_matrix, sampled_states, temp, nk_out, ukn_out):

    u_matrix = np.array(u_matrix)

    nk_out = os.path.abspath(nk_out)
    ukn_out = os.path.abspath(ukn_out)

    R = 8.3144598e-3

    M = len(u_matrix[:,0])
    N = len(u_matrix[0,:])

    N_k = np.zeros((M))
    u_kn = []

    for i in range(M):
        tmp_u_kn = []
        for j in range(N):
            U  = np.loadtxt(u_matrix[i,j], usecols=(1,), comments=['@','#'])
            pV = np.array([0.0 for u in U])
            for x,y in zip(U,pV):
                tmp_u_kn.append( (x+y)/(R*temp) )
            N_k[sampled_states[j]] = len(U)
        u_kn.append(tmp_u_kn)
    
    u_kn = np.array(u_kn)
        
    np.savetxt(nk_out, N_k)
    np.savetxt(ukn_out, u_kn)

    return (N_k, u_kn)

def prepare_mbar (u_matrix, pv_matrix, sampled_states, temp, nk_out, ukn_out):

    u_matrix = np.array(u_matrix)
    pv_matrix = np.array(pv_matrix)

    nk_out = os.path.abspath(nk_out)
    ukn_out = os.path.abspath(ukn_out)

    R = 8.3144598e-3

    M = len(u_matrix[:,0])
    N = len(u_matrix[0,:])

    N_k = np.zeros((M))
    u_kn = []

    for i in range(M):
        tmp_u_kn = []
        for j in range(N):
            U  = np.loadtxt(u_matrix[i,j], usecols=(1,), comments=['@','#'])
            pV = np.loadtxt(pv_matrix[i,j], usecols=(1,), comments=['@','#'])
            for x,y in zip(U,pV):
                tmp_u_kn.append( (x+y)/(R*temp) )
            N_k[sampled_states[j]] = len(U)
        u_kn.append(tmp_u_kn)
    
    u_kn = np.array(u_kn)
        
    np.savetxt(nk_out, N_k)
    np.savetxt(ukn_out, u_kn)

    return (N_k, u_kn)

# modified yMHG qua abr 29 19:47:46 -03 2020
# now prop_matrix is NPROPS x NSTATES x NSAMPLEDSTATES
# prop_matrix[i,j,k] = trajectory of property i for sampled state k evaluated at state j
# in practice something like reweighted_properties/property[i]_k_j.xvg
def estimate_properties (u_matrix, pv_matrix, sampled_states, temp, nk_out, ukn_out, eff_out, eig_out, mat_out, prop_matrix, out_preffixes):

    u_matrix = np.array(u_matrix)
    if (pv_matrix != ""):
        pv_matrix = np.array(pv_matrix)
    prop_matrix = np.array(prop_matrix)

    for i in range(len(out_preffixes)):
        out_preffixes[i] = os.path.abspath(out_preffixes[i])
        path_of_preffix = '/'.join(out_preffixes[i].split('/')[0:-1])
        # create path if it does not exist
        os.system("mkdir -p " + path_of_preffix)

    eff_out = os.path.abspath(eff_out)
    path_of_preffix = '/'.join(eff_out.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)

    path_of_preffix = '/'.join(nk_out.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)

    path_of_preffix = '/'.join(ukn_out.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)

    M = prop_matrix.shape[0] # number of properties
    K = prop_matrix.shape[1] # total number of states
    N = prop_matrix.shape[2] # number of sampled states
    
    (N_k, u_kn) = prepare_mbar (u_matrix, pv_matrix, sampled_states, temp, nk_out, ukn_out)

    mbar = pymbar.MBAR (u_kn, N_k)
    Eff_k = mbar.computeEffectiveSampleNumber()

    over = mbar.computeOverlap()
    Eig_k = over['eigenvalues']
    Mat_k = over['matrix']

    print "MBAR WEIGHTS"
    print mbar.getWeights()

    np.savetxt(eff_out, Eff_k)
    np.savetxt(eig_out, Eig_k)
    np.savetxt(mat_out, Mat_k)

    # 2d-list holding the paths of the files containing the results (average, uncertainty and effective samples)
    out = [["" for j in range(2)] for i in range(M)]
    out_values = [[0.0 for j in range(2)] for i in range(M)]

    # yMHG qua abr 29 19:56:39 -03 2020
    # each A_kn matrix used in MBAR must be NSTATES x NCONFIGURATION
    # i.e. A_kn[i,j] is property A evaluated at state i for j-th configuration
    # for each property
    for i in range(M):
        A_kn = []
        # for each state
        for j in range(K):
            P = []
            # for each sampled state
            for k in range(N):
                P += np.loadtxt(prop_matrix[i,j,k], usecols=(1,), comments=['#','@']).tolist()
            # now P contains all reweighted samples concatenated
            A_kn.append(P)
            
        A_kn = np.array(A_kn)
        mbarComputeOutput = mbar.computeExpectations(A_kn, state_dependent=True)

        EA_k = mbarComputeOutput['mu']
        dEA_k = mbarComputeOutput['sigma']

        # write files with results
        np.savetxt (out_preffixes[i] + "_EA_k.dat", EA_k)
        np.savetxt (out_preffixes[i] + "_dEA_k.dat", dEA_k)

        # fill output array
        out[i][0] = out_preffixes[i] + "_EA_k.dat"
        out[i][1] = out_preffixes[i] + "_dEA_k.dat"
        out_values[i][0] = EA_k
        out_values[i][1] = dEA_k

    return np.array(out_values)


# modified for reweighting of properties by yMHG qua abr 29 21:08:30 -03 2020
def estimate_properties_no_pv (u_matrix, sampled_states, temp, nk_out, ukn_out, eff_out, eig_out, mat_out, prop_matrix, out_preffixes):

    u_matrix = np.array(u_matrix)
    prop_matrix = np.array(prop_matrix)

    for i in range(len(out_preffixes)):
        out_preffixes[i] = os.path.abspath(out_preffixes[i])
        path_of_preffix = '/'.join(out_preffixes[i].split('/')[0:-1])
        # create path if it does not exist
        os.system("mkdir -p " + path_of_preffix)

    eff_out = os.path.abspath(eff_out)
    path_of_preffix = '/'.join(eff_out.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)

    path_of_preffix = '/'.join(nk_out.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)

    path_of_preffix = '/'.join(ukn_out.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)

    M = prop_matrix.shape[0] # number of properties
    K = prop_matrix.shape[1] # total number of states
    N = prop_matrix.shape[2] # number of sampled states 

    (N_k, u_kn) = prepare_mbar_no_pv (u_matrix, sampled_states, temp, nk_out, ukn_out)

    mbar = pymbar.MBAR (u_kn, N_k)
    Eff_k = mbar.computeEffectiveSampleNumber()
    over = mbar.computeOverlap()
    Eig_k = over['eigenvalues']
    Mat_k = over['matrix']

    print "MBAR WEIGHTS"
    print mbar.getWeights()

    np.savetxt(eff_out, Eff_k)
    np.savetxt(eig_out, Eig_k)
    np.savetxt(mat_out, Mat_k)

    # 2d-list holding the paths of the files containing the results (average, uncertainty and effective samples)
    out = [["" for j in range(2)] for i in range(M)]
    out_values = [[0.0 for j in range(2)] for i in range(M)]

    # yMHG qua abr 29 19:56:39 -03 2020
    # each A_kn matrix used in MBAR must be NSTATES x NCONFIGURATION
    # i.e. A_kn[i,j] is property A evaluated at state i for j-th configuration
    # for each property
    for i in range(M):
        A_kn = []
        # for each state
        for j in range(K):
            P = []
            # for each sampled state
            for k in range(N):
                P += np.loadtxt(prop_matrix[i,j,k], usecols=(1,), comments=['#','@']).tolist()
            # now P contains all reweighted samples concatenated
            A_kn.append(P)
              
        A_kn = np.array(A_kn)
        mbarComputeOutput = mbar.computeExpectations(A_kn, state_dependent=True)

        EA_k = mbarComputeOutput['mu']
        dEA_k = mbarComputeOutput['sigma']

        # write files with results
        np.savetxt (out_preffixes[i] + "_EA_k.dat", EA_k)
        np.savetxt (out_preffixes[i] + "_dEA_k.dat", dEA_k)

        # fill output array
        out[i][0] = out_preffixes[i] + "_EA_k.dat"
        out[i][1] = out_preffixes[i] + "_dEA_k.dat"
        out_values[i][0] = EA_k
        out_values[i][1] = dEA_k

    return np.array(out_values)

