############
Command-line
############

.. code-block::

    usage: gmak [-h] [--validate] [--gmx GMXPATH] [--gnp NPROCS]
                [--restart BINPATH] [--custom]
                INPUT
    
    gmak is a tool for the optimization of force-field parameters aided by a grid-
    based mapping of the parameter-search space and by the use of surrogate
    models. The optimization targets the reproduction of reference (experimental
    or theoretical) values of (bio)physical properties. gmak also serves as a tool
    for the visualization and quantification of the influence of the force-field
    parameters on the target properties.
    
    positional arguments:
      INPUT              The path of the input file.
    
    optional arguments:
      -h, --help         show this help message and exit
      --validate         Activate validation mode.
      --gmx GMXPATH      Path of the gmx binary. If not given, it is inferred
                         using the command `which'.
      --gnp NPROCS       Number of parallel threads used in GROMACS simulations
                         (option `-nt' of `mdrun').
      --restart BINPATH  Path of the binary state file.
      --custom           Read the customization settings from a file named
                         `custom.py' in the current directory.
