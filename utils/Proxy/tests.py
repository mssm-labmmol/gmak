#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import custom
import objective_function
import gridoptimizer
import estimator
import grid
import ext_optimizers


def sampler_factory(id: int, std_dev: float) -> estimator.BaseEstimator:
    if id == 0:
        return estimator.NoisyEstimator(func_avg=lambda x: 900 + np.sum(np.array(x)**2),
                                        func_std=lambda x: std_dev,
                                        nsamples=200)
    elif id == 1:
        return estimator.NoisyEstimator(func_avg=lambda x: 900 + 0.26*np.sum(np.array(x)**2) - 0.48*x[0]*x[1],
                                           func_std=lambda x: std_dev,
                                           nsamples=200)
    else:
        raise ValueError




if __name__ == "__main__":

    std_devs = [0.0, 0.3, 3.0]

    # Part 1: [2,4] x [2,4], comparing different GMAKs for different Std. Devs.
    fp = open("test_results.log", "w")
    for func_id in [0, 1]:
        for std_dev in std_devs:
            for method in ["linear", "cubic", "gpr"]:
                min_x_arr = []
                min_val_arr = []
                ncalls_arr = []
                for repeat in range(5):
                    grid_lower = (2.0, 2.0)
                    grid_upper = (4.0, 4.0)
                    grid_spacing = (0.1, 0.1)
                    opt_grid = grid.Grid.from_limits_spacing(grid_lower,
                                                             grid_upper,
                                                             grid_spacing)
                    samples = opt_grid.extract_pattern()
                    sampler = sampler_factory(func_id, std_dev)
                    if method == "gpr":
                        estim = estimator.GPR(samples, sampler, lambda x: x)
                    else:
                        estim = estimator.Interpolator(opt_grid,
                                                       samples,
                                                       sampler,
                                                       method)
                    obj_func = objective_function.RealObjectiveFunction(estim)
                    gmak = gridoptimizer.StandardGridMaker(grid=opt_grid,
                                                           objective_function=obj_func,
                                                           shift_cut=0.01,
                                                           maxshifts=10,
                                                           conv_margins=[(-1.15, 0.85),
                                                                         (0.15, 0.85)])
                    gmak.run()
                    min_x, min_val = gmak.get_min()
                    min_x_arr.append(min_x)
                    min_val_arr.append(min_val)
                    ncalls_arr.append(gmak.get_ncalls())

                min_x_arr = np.array(min_x_arr)
                min_val_arr = np.array(min_val_arr)
                ncalls_arr = np.array(ncalls_arr)
                fp.write("%10d%10.2f%10s%10.4f+/-%10.4f%10.4f+/-%10.4f%10.4f+/-%10.4f%10.4f+/-%10.4f\n" % (func_id, std_dev, method,
                                                                                                                          np.mean(min_x_arr[:,0]), np.std(min_x_arr[:,0]),
                                                                                                                          np.mean(min_x_arr[:,1]), np.std(min_x_arr[:,1]),
                                                                                                                          np.mean(min_val_arr), np.std(min_val_arr),
                                                                                                                          np.mean(ncalls_arr), np.std(ncalls_arr)))

    fp.close()




    #obj_func_sampled = objective_function.RealObjectiveFunction(sampler)

    ## GMAK



    ##gmak.run()

    ##print(gmak.get_min())
    ##print("ncalls=", gmak.objective_function.get_estimators()[0].sampler.ncalls)

    ### SIMPLEX
    ##simplex = ext_optimizers.SimplexOptimizer(obj_func_sampled, (3.0, 3.0))
    ##simplex_out = simplex.run()
    ##print(simplex_out)

    ### CG
    ##cg = ext_optimizers.CGOptimizer(obj_func_sampled, (3.0, 3.0))
    ##cg_out = cg.run()
    ##print(cg_out)

    ## BasinHop
    #basinhop = ext_optimizers.BasinHopOptimizer(obj_func_sampled, (3.0, 3.0))
    #basinhop_out = basinhop.run()
    #print("BasinHop results", basinhop_out)

    ## DE
    #de = ext_optimizers.DEOptimizer(obj_func_sampled, 6, [(-4,4),(-4,4)])
    #de_out = de.run()
    #print("DE results", de_out)


