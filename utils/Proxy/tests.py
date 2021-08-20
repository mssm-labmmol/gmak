#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np
import custom
import objective_function
import gridoptimizer
import estimator
import grid
import ext_optimizers

if __name__ == "__main__":

    opt_grid = grid.Grid.from_limits_spacing((2., 2.),
                                             (4., 4.),
                                             (0.1, 0.1))
    samples = opt_grid.extract_pattern()

    #sampler = estimator.NoisyEstimator(func_avg=lambda x: np.sum(np.array(x)**2),
    #                                   func_std=lambda x: 0.01*np.sum(np.array(x)**2),
    #                                   nsamples=200)
    sampler = estimator.NoisyEstimator(func_avg=lambda x: 900+0.26*np.sum(np.array(x)**2) - 0.48*x[0]*x[1],
                                       func_std=lambda x: 0.30, #*(0.26*np.sum(np.array(x)**2) - 0.48*x[0]*x[1]),
                                       nsamples=200)

    #estim = estimator.Interpolator(opt_grid,
    #                               samples,
    #                               sampler,
    #                               'linear')
    #estim = estimator.GPR(samples, sampler, lambda x: x)
    #obj_func = objective_function.RealObjectiveFunction(estim)

    obj_func_sampled = objective_function.RealObjectiveFunction(sampler)


    # GMAK
    #gmak = gridoptimizer.StandardGridMaker(grid=opt_grid,
    #                                       objective_function=obj_func,
    #                                       shift_cut=0.01,
    #                                       maxshifts=10,
    #                                       conv_margins=[(0.15, 0.85),
    #                                                     (0.15, 0.85)])



    #gmak.run()

    #print(gmak.get_min())
    #print("ncalls=", gmak.objective_function.get_estimators()[0].sampler.ncalls)

    # SIMPLEX
    simplex = ext_optimizers.SimplexOptimizer(obj_func_sampled, (3.0, 3.0))
    simplex_out = simplex.run()
    print(simplex_out)

    ## CG
    #cg = ext_optimizers.CGOptimizer(obj_func_sampled, (3.0, 3.0))
    #cg_out = cg.run()
    #print(cg_out)

