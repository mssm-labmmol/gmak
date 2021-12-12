.. gridmaker documentation master file, created by
   sphinx-quickstart on Fri Nov  5 16:38:43 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to gridmaker's documentation!
=====================================

.. note:: While conceived as a tool for parameter optimization,
   ``gmak`` also serves as a tool for visualizing and quantifying the
   influence of the parameters on the desired properties, as well as a
   general driver tool for running simulations and calculating
   properties--what I mean by this is that the user can create a
   customized simulation and analysis environment--define its own
   property-calculation routines, simulation routines, etc., and
   easily recycle them via ``gmak`` to simply carry out simulations
   and analyzes, and not necessarily for parameter optimization (e.g.,
   using a 1x1 grid to ignore parameter modifications; or using
   explicit grids to compare different models of the same system).

.. toctree::
   :maxdepth: 3
   :caption: Overview:

   overview/general_workflow
   overview/grid
   overview/interaction_parameters
   overview/systems_and_topologies
   overview/coordinates
   overview/protocols
   overview/properties
   overview/surrogate_model
   overview/score
   overview/grid_shifting

.. toctree::
   :maxdepth: 3
   :caption: Program Usage:

   usage/command-line
   usage/customization_api
   usage/input_file
   usage/output_files
   usage/post_processing
   usage/restarting_a_run
   usage/validation_mode

.. toctree::
   :maxdepth: 3
   :caption: Example Applications:

   examples/custom_property
   examples/custom_package


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


