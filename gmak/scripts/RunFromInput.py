#!/usr/bin/env python3

import numpy as np
import os
import sys
import copy
import pickle

import gmak.logger as logger
from gmak.gridshifter import *
from gmak.parameter_space import *
from gmak.grid_ana import *
from gmak.read_input import *
from gmak.mdputils import mdpUtils
from gmak.gridbase import GridPoint, ParameterGrid
from gmak.gridoptimizer import gridOptimizer
from gmak.state import *
from gmak.write_input import *
from gmak.results_assembler import *
from gmak.gridpoint_selector import *

import gmak.config
import gmak.runcmd as runcmd
import gmak.simulate as simulate

if ('--legacy' in sys.argv):
    bool_legacy = True
else:
    bool_legacy = False

# Run-time custom files.
if os.path.isfile("custom.py"):
    sys.path.insert(0, '.')
    import custom

#plotFlag = not ('--no-plot' in sys.argv)
plotFlag = False # never plot
validateFlag = ('--validate' in sys.argv)
writeNewInput = ('--write-input' in sys.argv)

# -----
# Gmx Path
# -----
try:
    gmak.config.ConfigVariables.gmx = sys.argv[sys.argv.index('--gmx') + 1]
except ValueError:
    try:
        import subprocess
        # get gmx path from 'which gmx'
        gmak.config.ConfigVariables.gmx = subprocess.check_output(['which', 'gmx']).strip().decode('utf-8')
    except:
        # Maybe gmx is not needed for the job?
        pass


try:
    cmdNprocs = sys.argv[sys.argv.index('-np') + 1]
    simulate.mdrun_nprocs = int(cmdNprocs)
except ValueError:
    pass

def main():
    if ('--restart' in sys.argv): # restart from a given binary file
        binFilename = sys.argv[sys.argv.index('--restart') + 1]
        binFile = open(binFilename, 'rb')
        globalState.setFromBinary(binFile) # load state
        binFile.close()

        (base_workdir, # read variables as if reading from an input file,
         grid,         # but actually reading from state
         protocols,    #
         properties,   # this way the rest of this program can remain the same
         protocolsHash,
         optimizer,
         surrogateModelHash) = globalState.getInitializationState()

        logger.globalLogger.putMessage('Restarting from {}'.format(binFilename))

        # Merge new settings from input file, if it is given
        if not sys.argv[1].startswith('-'):
            # Set a merging state
            mergeState = State()
            mergeState.setFromInput(initialize_from_input(sys.argv[1], bool_legacy, validateFlag))
            globalState.merge(mergeState)

    else:
        globalState.setFromInput(initialize_from_input(sys.argv[1], bool_legacy, validateFlag))

        (base_workdir, # read variables from input file
         grid,
         protocols,
         properties,
         protocolsHash,
         optimizer,
         surrogateModelHash) = globalState.getInitializationState()

    # Initialize logger.globalLogger with output file inside the workdir
    try:
        os.mkdir(base_workdir)
    except FileExistsError:
        pass
    logger.globalLogger = logger.Logger(
        os.path.join(
            base_workdir,
            "gmak_{}.dat".format(os.getpid())))

    # Initialize results assembler - hard-coded for now, will go into input file
    resultsAssembler = ResultsAssembler(grid.getParameterNames())

    # Initialize selector - hard-coded for now, will go into input file
    selector = createSelector('best', gridOptimizer=optimizer, howMany=9)

    logger.globalLogger.putMessage('BEGIN MAINLOOP', dated=True)
    logger.globalLogger.indent()
    grid.run(protocols, optimizer, surrogateModelHash, properties, protocolsHash, resultsAssembler, plotFlag)
    logger.globalLogger.unindent()
    logger.globalLogger.putMessage('END MAINLOOP', dated=True)

    # DEBUG
    resultsAssembler.print(sys.stdout)

    # Write assembled results to a file
    resultsAssembler.writeToFile(grid.makeCurrentWorkdir() + "/assembled")


    # *********************** Write new input **********************************

    if (writeNewInput):
        mod = InputModifier(sys.argv[1])
        #
        mod.set_workdir(base_workdir + "_best_points")
        mod.set_main_variation(grid.parSpaceGen)
        mod.set_samples(selector.selectPoints())
        mod.write_to_file(sys.argv[1] + "_best_points")

    # *********************** Subgrid part           **************************        

    #subgridHash['parspacegen'] = copy.deepcopy(subgridHash['parspacegen'])
    #subgridHash['parspacegen'].rescale(subgridHash['factors'])
    #subgridObj = ParameterSubgrid(subgridHash['parspacegen'], grid, properties, subgridHash['method'], bool_legacy)
    #for prop in properties:
    #    refValue = optimizer.referenceValues[prop]
    #    subgridObj.save_property_values_to_file(prop)
    #    subgridObj.save_property_diff_to_file(prop, refValue)
    #    subgridObj.save_property_err_to_file(prop)
    #    if plotFlag:
    #        subgridObj.plot_property_to_file(prop)
    #        subgridObj.plot_property_diff_to_file(prop, properties[prop], refValue)
    #        subgridObj.plot_property_err_to_file(prop)
    #subgridObj.printAndPlotScores(optimizer, plotFlag)
    #subgridObj.writeParameters()
    #subgridObj.save_to_binary(optimizer)

if __name__ == "__main__":
    main()
