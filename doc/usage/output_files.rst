############
Output Files
############

Basic file structure
====================

The file structure below depicts the output files of a ``gmak`` job.
The highlighted tokens that start with the ``%`` character may
represent a unique value (e.g., ``%workdir``) or several possible
values (e.g., ``%protocol``), each yielding an individual file or
directory. The asterisk indicates files that are created only under
special circumstances. The tokens and files are explained in more
details below.

.. code-block:: tree

    %workdir
    ├── %system 
    │   ├── ...
    ├── gmak_%jobid.log
    ├── state_%jobid.bin
    ├── grid_%iter
    │   ├── %protocol
    │   │   ├── %surrogatemodel-in-protocol
    │   │   │   └── estimated_properties
    │   │   │       ├── %componentprop-in-surrogatemodel_dEA_k.dat
    │   │   │       ├── %componentprop-in-surrogatemodel_EA_k.dat
    │   │   └── simu
    │   │       ├── %sampledgridpoint
    │   │       │   ├── %componentprop-in-protocol.xvg
    │   │       │   ├── *filtered_%componentprop-in-protocol.xvg
    │   │       │   └── %coords-in-protocol
    │   ├── parameters_%variation.dat
    │   ├── samples_0.dat
    │   └── step_0
    │       ├── %compositeprop_EA_k.dat
    │       ├── %compositeprop_dEA_k.dat
    │       ├── %compositeprop_diff.dat
    │       ├── full_data.dat
    │       ├── *full_data.dat.cis
    │       ├── optimizer_data.dat
    │       └── *optimizer_data.dat.cis

Tokens
------

``%workdir``
    The name of the work directory of the job (the ``workdir``
    parameter in the :doc:`input file </usage/input_file>`).

``%system``
    The name of a system (defined in a :doc:`system block
    </usage/blocks/system>` of the :doc:`input file
    </usage/input_file>`). In the file structure, it corresponds to a
    directory designed to contain the topology files of the
    corresponding system, for all grid-shift iterations and grid
    points.

``%jobid``
    The PID of the ``gmak`` job.

``%iter``
    A grid-shift iteration.

``%protocol``
    The name of a protocol (defined in a :doc:`protocol block
    </usage/blocks/protocol>` of the :doc:`input file
    </usage/input_file>`).

``%surrogatemodel-in-protocol``
    The name of a surrogate model associated with a protocol, i.e.,
    used to calculate at least one component property associated with
    the protocol (defined in a :doc:`compute block
    </usage/blocks/compute>` of the :doc:`input file
    </usage/input_file>`).

``%componentprop-in-surrogatemodel``
    The name of a component property associated a surrogate model
    (defined in a :doc:`compute block </usage/blocks/compute>` of the
    :doc:`input file </usage/input_file>`).

``%sampledgridpoint``
    The linear index of a sampled grid point.

``%componentprop-in-protocol``
    The name of a component property associated with a protocol.

``%coords-in-protocol``
    The path of the initial coordinates associated with a protocol
    (defined in a :doc:`coordinates block </usage/blocks/coordinates>`
    of the :doc:`input file </usage/input_file>`).

``%compositeprop``
    The name of a composite property (defined in a :doc:`compute block
    </usage/blocks/compute>` of the :doc:`input file
    </usage/input_file>`).

``%variation``
    The name of a variation (defined in a :doc:`variation block
    </usage/blocks/variation>` of the :doc:`input file
    </usage/input_file>`).

Files
-----

``gmak_%jobid.log``
    The :ref:`log file <usage/output_files:Log file>` of the job.

``state_%jobid.bin``
    The state of the job stored as a binary file. This contains
    several data that can be used to :doc:`restart a run
    </usage/restarting_a_run>` or perform :doc:`post-processing analysis </usage/post_processing>`.

``%componentprop-in-surrogatemodel_dEA_k.dat``
    The statistical uncertainties of the component property for all
    grid points, as estimated by the surrogate model.

``%componentprop-in-surrogatemodel_EA_k.dat``
    The expected values of the component property for all grid points,
    as estimated by the surrogate model.

``%componentprop-in-protocol.xvg``
    The values of the component property, as estimated from the
    protocol simulations. This can be a timeseries or a file
    containing the expected value and uncertainty of the property.

``*filtered_%componentprop-in-protocol.xvg``
    If the component property is timeseries-based, this file contains
    the subsambled data after statistical processing to remove
    auto-correlation.

``%coords-in-protocol``
    The path of the initial coordinates associated with a protocol
    (defined in a :doc:`coordinates block </usage/blocks/coordinates>`
    of the :doc:`input file </usage/input_file>`).

``parameters_%variation.dat``
    The values of the force-field parameters of the variation for each
    grid point.

``samples_0.dat``
    The linear indexes of the sampled grid points. The suffix ``_0``
    has no special meaning in the current version of the program.

``%compositeprop_EA_k.dat``
    The expected values of the composite property for all grid points.

``%compositeprop_dEA_k.dat``
    The uncertainties in the estimate of the composite property for
    all grid points.

``%compositeprop_diff.dat``
    The difference between the expected values and the reference value
    of the composite property for all grid points.

``full_data.dat``
    A summary file containing the expected values and uncertainties
    of all composite properties for all grid points, as well as the
    value of the score function.

``*full_data.dat.cis``
    If a score-uncertainty function is provided, this file contains
    the corresponding confidence intervals of the score for all grid
    points.

``optimizer_data.dat``
    This is the summary file ``full_data.dat`` ordered from smallest
    to largest value of the score function.

``*optimizer_data.dat.cis``
    This is the file ``full_data.dat.cis`` ordered from smallest to
    largest value of the score function.


GROMACS-compatible Systems
--------------------------

For GROMACS-compatible systems, the ``%system`` directory contains
files ``%system_%iter_%sampledgridpoint.top``, where the tokens
``%system``, ``%iter`` and ``%sampledgridpoint`` are explained above.


GROMACS-compatible General protocol
-----------------------------------

For the GROMACS-compatible general protocol, the ``%sampledgridpoint``
directories contain additional directories corresponding to each
simulation in the sequence of simulations of the protocol. These
directories are identified by an index, ``%stage``, and each one
contains the output files of the simulations, e.g., ``%stage.tpr``,
``%stage.xtc``, ``%stage.edr``, etc.


GROMACS-compatible Alchemical protocol
--------------------------------------

As explained :ref:`here <overview/protocols:gromacs-compatible
alchemical protocol>`, the GROMACS-compatible Alchemical protocol
implements a sequence of GROMACS-compatible General (sub)protocols for
each state in the alchemical transformation. Each of these
subprotocols corresponds to a protocol directory with name
``%protocol-%state``, where ``%protocol`` is the name of the
alchemical protocol and ``%state`` is the corresponding state of the
alchemical transformation. 


Log file
========

The progress of the job is reported in real time in the log file.  It
registers the beginning and end of all simulations, the results of
checking whether they need to be extended or not, the values of the
extended lengths (if existent), the result of the grid-shifting
procedure, and, for some surrogate models, complementary information
regarding the fitting of the model.


Stdout and Stderr
=================

In general, ``gmak`` itself does not write any information in the
stdout and stderr streams, except for error messages---logging is done
entirely in the :ref:`log file <usage/output_files:Log file>`, and
output data is written in the output files. Other processes
initialized by the program (e.g., GROMACS binaries), however, may
write in those streams.
