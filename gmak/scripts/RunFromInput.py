#!/usr/bin/env python3

from .. import _version
__version__ = _version.get_versions()['version']

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
import argparse
import subprocess

def find_gmx():
    try:
        return subprocess.check_output(['which', 'gmx']).strip().decode('utf-8')
    except:
        return None

program_description = """
gmak is a tool for the optimization of force-field parameters aided by a
grid-based mapping of the parameter-search space and by the use of surrogate
models. The optimization targets the reproduction of reference (experimental or
theoretical) values of (bio)physical properties.

gmak also serves as a tool for the visualization and quantification of the
influence of the force-field parameters on the target properties.
"""

def gmak_run(argumentNamespace):
    plotFlag = False
    bool_legacy = False
    writeNewInput = False
    binFilename = argumentNamespace.restart
    inputFilename = argumentNamespace.input
    validateFlag = argumentNamespace.validate
    customFlag = argumentNamespace.custom

    # Run-time custom files.
    if customFlag:
        if os.path.isfile("custom.py"):
            sys.path.insert(0, '.')
            import custom
        else:
            raise OSError("Could not find file custom.py.")

    # Set module/class attributes
    gmak.config.ConfigVariables.gmx = argumentNamespace.gmx
    if argumentNamespace.gnp is not None:
        simulate.mdrun_nprocs = argumentNamespace.gnp

    if binFilename:
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

        # Merge new settings from input file
        # Set a merging state
        mergeState = State()
        mergeState.setFromInput(initialize_from_input(inputFilename, bool_legacy, validateFlag))
        globalState.merge(mergeState)
    else:
        globalState.setFromInput(initialize_from_input(inputFilename, bool_legacy, validateFlag))

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
            "gmak_{}.log".format(os.getpid())))

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
    # resultsAssembler.print(sys.stdout)

    # Write assembled results to a file
    # resultsAssembler.writeToFile(grid.makeCurrentWorkdir() + "/assembled")


    # *********************** Write new input **********************************

    if (writeNewInput):
        mod = InputModifier(inputFilename)
        #
        mod.set_workdir(base_workdir + "_best_points")
        mod.set_main_variation(grid.parSpaceGen)
        mod.set_samples(selector.selectPoints())
        mod.write_to_file(inputFilename + "_best_points")

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

def main():
    # Read arguments
    parser = argparse.ArgumentParser(description=program_description)
    parser.add_argument("input",
                        metavar="INPUT",
                        type=str,
                        help="The path of the input file.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument("--validate",
                        action="store_true",
                        help="Activate validation mode.")
    parser.add_argument("--gmx",
                        metavar="GMXPATH",
                        type=str,
                        default=find_gmx(),
                        help="Path of the gmx binary. If not given, it is inferred using the command `which'.")
    parser.add_argument("--gnp",
                        metavar="NPROCS",
                        type=int,
                        help="Number of parallel threads used in GROMACS simulations (option `-nt' of `mdrun').")
    parser.add_argument("--restart",
                        metavar="BINPATH",
                        type=str,
                        help="Path of the binary state file.")
    parser.add_argument("--custom",
                        action="store_true",
                        help="Read the customization settings from a file named `custom.py' in the current directory.")
    args = parser.parse_args()
    # Run the job
    gmak_run(args)

if __name__ == "__main__":
    main()
