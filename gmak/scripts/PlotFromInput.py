#!/usr/bin/env python3
import gmak.runcmd as runcmd

# General imports.
import pickle
import sys
import os

# Custom imports.
from gmak.gridbase       import  ParameterGrid
from gmak.subgrid        import  ParameterSubgrid
from gmak.gridoptimizer  import  gridOptimizer
from gmak.read_input     import  initialize_from_input

def pickle_load_from_file(fn):
    fp = open(fn, 'rb')
    obj = pickle.load(fp)
    fp.close()
    return obj

if __name__ == '__main__':
    if (len(sys.argv) == 1):
        print("# USAGE: ./PlotFromInput.py input.inp grid.bin [--save-data]")
        exit(0)

    (base_workdir,
     grid,
     protocols,
     properties,
     protocolsHash,
     optimizer,
     surrogateModelHash,
     subgridHash) = initialize_from_input (sys.argv[1], False)
    
    grid_path = sys.argv[2]
    dirname = base_workdir
    save_data_bool = ('--save-data' in sys.argv)

    grid = pickle_load_from_file(grid_path)
    grid.xlabel = "$X$"
    grid.ylabel = "$Y$"
    optimizer = grid.optimizer

    grid.resetWorkdir(dirname)

    for prop in properties:
        reference_value = optimizer.referenceValues[prop]
        if type(grid) is ParameterGrid:
            grid.plot_property_to_file(prop, optimizer)
            grid.plot_property_err_to_file(prop, optimizer)
            grid.plot_property_diff_to_file(prop, properties[prop], reference_value, optimizer)
            if (save_data_bool):
                grid.save_property_values_to_file(prop, optimizer)
                grid.save_property_err_to_file(prop, optimizer)
                grid.save_property_diff_to_file(prop, reference_value, optimizer)
        elif type(grid) is ParameterSubgrid:
            grid.plot_property_to_file(prop)
            grid.plot_property_err_to_file(prop)
            grid.plot_property_diff_to_file(prop, properties[prop], reference_value)
            if (save_data_bool):
                grid.save_property_values_to_file(prop)
                grid.save_property_err_to_file(prop)
                grid.save_property_diff_to_file(prop, reference_value)

    if (save_data_bool):
        optimizer.printToFile(grid, grid.makeStepPropertiesdir(optimizer) + "/optimizer_data.dat")
        optimizer.printToFile (grid, grid.makeStepPropertiesdir(optimizer) + "/full_data.dat", sorted=False)

    optimizer.plotToPdf(grid, grid.makeStepPropertiesdir(optimizer) + "/optimizer_score.pdf")
