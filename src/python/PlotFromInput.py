#!/usr/bin/env python3

# General imports.
import pickle
import sys
import os

# Custom imports.
from  gridbase       import  ParameterGrid
from  subgrid        import  ParameterSubgrid
from  gridoptimizer  import  gridOptimizer
from  read_input     import  initialize_from_input

def pickle_load_from_file(fn):
    fp = open(fn, 'rb')
    obj = pickle.load(fp)
    fp.close()
    return obj

if __name__ == '__main__':
    if (len(sys.argv) == 1):
        print("# USAGE: ./PlotFromInput.py input.inp grid.bin optimizer.bin")
        exit(0)

    (base_workdir,
     grid,
     protocols,
     properties,
     protocolsHash,
     optimizer,
     gaCoverInterface,
     surrogateModelHash,
     subgridHash) = initialize_from_input (sys.argv[1], False)
    
    grid_path = sys.argv[2]
    optimizer_path = sys.argv[3]
    dirname = os.path.dirname(grid_path)

    grid = pickle_load_from_file(grid_path)
    optimizer = pickle_load_from_file(optimizer_path)

    for prop in properties:
        reference_value = optimizer.referenceValues[prop]
        if type(grid) is ParameterGrid:
            grid.plot_property_to_file(prop, optimizer)
            grid.plot_property_err_to_file(prop, optimizer)
            grid.plot_property_diff_to_file(prop, properties[prop], reference_value, optimizer)
        elif type(grid) is ParameterSubgrid:
            grid.plot_property_to_file(prop)
            grid.plot_property_err_to_file(prop)
            grid.plot_property_diff_to_file(prop, properties[prop], reference_value)

    optimizer.plotToPdf(grid, dirname + "/optimizer_score.pdf")
