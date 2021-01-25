#!/usr/bin/env python3

from logger import *
from gridshifter import * 
from parameters import * 
from grid_ana import * 
from coverage import *
from read_input import *
from mdputils import mdpUtils
from protocols import LiquidProtocol, GasProtocol, SlabProtocol
from gridbase import GridPoint, ParameterGrid
from gridoptimizer import gridOptimizer
from subgrid import *

import numpy as np
import os
import sys
import copy
import runcmd

if ('--legacy' in sys.argv):
    bool_legacy = True
else:
    bool_legacy = False

plotFlag = not('--no-plot' in sys.argv)

if __name__ == "__main__":

    (base_workdir,
     grid,
     protocols,
     properties,
     protocolsHash,
     optimizer,
     gaCoverInterface,
     surrogateModelHash,
     subgridHash) = initialize_from_input (sys.argv[1], bool_legacy)

    # push nearest copies into optimizer
    optimizer.pushNearest(surrogateModelHash)

    # SET RANDOM NUMBER SEED - TODO

    # In each protocol, check writing frequencies of energy and trajectory to avoid any problems
    mdp_ut = mdpUtils()
    for protocol in protocols:
        # Check only production run.
        mdp_ut.parse_file(protocol.mdps[-1])
        efreq = mdp_ut.get_writefreq_energy()
        if (protocol.type == 'slab'):
            xfreq = mdp_ut.get_writefreq_pos()
        else:
            xfreq = mdp_ut.get_writefreq_compressedpos()
        if (efreq != xfreq):
            raise ValueError("For protocol {}, efreq ({}) and xfreq ({}) don't match.".format(
                protocol.name, efreq, xfreq))

    # *********************** End of checks for run  **************************        

    globalLogger.putMessage('BEGIN MAINLOOP', dated=True)
    globalLogger.indent()
    grid.run(protocols, optimizer, surrogateModelHash, properties, protocolsHash, plotFlag)
    globalLogger.unindent()
    globalLogger.putMessage('END MAINLOOP', dated=True)

    # *********************** Subgrid part           **************************        

    subgridHash['parspacegen'] = copy.deepcopy(subgridHash['parspacegen'])
    subgridHash['parspacegen'].rescale(subgridHash['factors'])
    subgridObj = ParameterSubgrid(subgridHash['parspacegen'], grid, properties, subgridHash['method'], bool_legacy)
    for prop in properties:
        refValue = optimizer.referenceValues[prop]
        subgridObj.save_property_values_to_file(prop)
        subgridObj.save_property_diff_to_file(prop, refValue)
        subgridObj.save_property_err_to_file(prop)
        if plotFlag:
            subgridObj.plot_property_to_file(prop)
            subgridObj.plot_property_diff_to_file(prop, properties[prop], refValue)
            subgridObj.plot_property_err_to_file(prop)
    subgridObj.printAndPlotScores(optimizer, plotFlag)
    subgridObj.writeParameters()
    subgridObj.save_to_binary(optimizer)
