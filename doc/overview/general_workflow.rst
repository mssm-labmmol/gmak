################
General Workflow
################

Optimization problem
====================

The general optimization problem to which ``gmak`` has been designed
as a solution strategy is:

.. math::

    \text{minimize}\quad (s \circ \mathbf{P}) (\mathbf{x; \mathbf{r}})

where

- :math:`\mathbf{x} \in \mathbb{R}^D` is the independent variable and
  is associated with the values of the adjusted force-field parameters
- :math:`s: \mathbb{R}^{N_P} \to \mathbb{R}` is the so-called score
  function, which quantifies the deviation of values
  :math:`(p_1,\ldots,p_{N_P})` of :math:`N_P` distinct target
  properties from some given reference values
  :math:`\mathbf{r}=(r_1,\ldots,r_{N_P})`

- :math:`\mathbf{P}:\mathbb{R}^D\to\mathbb{R}^{N_P}` is a  function
  with random component functions
  :math:`P_i:\mathbb{R}^D\to\mathbb{R}`, :math:`i=1,\ldots,N_P`

In the context of molecular mechanics, the value
:math:`P_i(\mathbf{x})` is typically associated with the expected
value of the :math:`i`-th target property for the
force-field-parameter values :math:`\mathbf{x}`.  It is a random
variable, because it is estimated from simulations, which have
statistical uncertainties associated with them. As a consequence,
:math:`(s\circ\mathbf{P})(\mathbf{x};\mathbf{r})` is a random variable
that quantifies the accuracy of the parameter values
:math:`\mathbf{x}` in reproducing the reference values of the
:math:`N_P` target properties of which the calculation is embedded in the
functions :math:`P_i`.

Solution strategy
=================

The solution strategy adopted in ``gmak`` stands on the following key
aspects.

Discretization of the Parameter Space
    The admissible values of the variable :math:`\mathbf{x}` are
    restricted to a universe that is discrete in the underlying
    parameter space.

Floating Search-space Boundaries
    At any given point in execution time, only a finite, bounded
    subset of the variable universe is available to the program---the
    parameter-search grid. However, the accessible subset is
    periodically updated, effectively moving over the parameter space
    in the general direction that is favorable to solving the
    optimization problem.

Surrogate Model
    The program relies on a surrogate model to numerically estimate
    the target properties for all accessible parameter search-space
    points based on the simulation of only a subset of those points.


Important terms
---------------

For clarity, some important terms are given an introductory definition
prior to describing the general workflow of ``gmak``.  Those are:

Input File
    The input file contains the instructions and data necessary to run
    the program from scratch. It follows the syntax described in 
    :doc:`/usage/input_file`.

Restart File
    The restart file contains the program state at a given point in
    time saved in binary format. It can be used to restart an
    interrupted program run.
    
Customization File
    The customization file provides the implementation of customized
    functions and should include those in the program by means of the
    :doc:`/usage/customization_api`.

(Parameter-search) Grid
    The parameter-search grid is the outcome of discretizing the
    parameter-search space. It is a multidimensional array where the
    elements are the explored values of the force-field parameters
    chosen as variables in the parameter-optimization problem.

    The parameter-search grid is explained more formally in
    :doc:`/overview/grid`.

Grid Point
    The grid points (or grid cells) are the elements of the
    parameter-search grid, i.e. vectors where each element is 
    a component of the variable :math:`\mathbf{x}`.

Grid Origin
    The grid origin is the grid point that is the zeroth element of
    the parameter-search grid.

Grid Shifting
    The grid-shifting procedure modifies the parameter-search grid
    such that its origin is altered without modification of the
    underlying rule that binds the grid points together. It is the
    manner in which the floating behavior of the parameter-search
    space boundaries is implemented. It is explained more formally in
    :doc:`/overview/grid`.

    The grid-shifting procedure can be customized by means of the
    :ref:`customization API <usage/customization_api:grid shifting>`.

Score Function
    The score function is the outmost component function :math:`s` of
    the objective function of the parameter-optimization problem (as
    explained above).

    The score function can be customized by means of the
    :ref:`customization API <usage/customization_api:score function>`.

Composite Property
    The composite properties are those properties of which the values
    are directly used in the computation of the score function.

    Routines for the calculation of composite properties can be
    customized by means of the :ref:`customization API
    <usage/customization_api:properties>`.

Component Property
    The component properties (for a given composite property) are any
    intermediary properties of which the values are used to compose
    the value of the corresponding composite property by means of a
    mathematical function.

    Routines for the calculation of component properties can be
    customized by means of the :ref:`customization API
    <usage/customization_api:properties>`.

(Composite) Property Tolerance
    The tolerance of a composite property is an upper bound for the
    uncertainty in its estimated value. The violation of the tolerance
    signals the need to extend the length of the simulations that
    underlie the calculation of the property.

Surrogate Model
    The surrogate models compute estimates, and their corresponding
    uncertainties, for the values of component properties for the
    entire parameter-search grid, based on their corresponding
    expected values and statistical errors for simulated grid points.

    Surrogate models can be customized by means of the
    :ref:`customization API <usage/customization_api:surrogate
    model>`.


Workflow steps
--------------


*Step 0:* Preparation
    In this phase, the issued command is parsed for the path of the
    input file and restart files, as well as for program options and
    their values.

    Furthermore, the customization file ``custom.py`` is searched for
    and, if existent, interpreted by the program.

*Step 1:* Object Initialization
    In this phase, the input or restart file is read, and the objects
    are initialized with their types and data as therein indicated.
    
    As far as the program is concerned, the run is fully determined
    after this phase has completed.

*Step 2:* Topology Writing
    In this phase, the topology files are written to disk, reflecting
    the values of the explored parameters at each grid point.

*Step 3:* Simulations and Extensions
    In this phase, all necessary simulations, including extensions,
    are launched in sequence.
    
    A given simulation is deemed as not necessary if the restart file
    indicates so, if it is flagged as complete for the run at
    hand, or if the program finds evidence in disk that the simulation
    has already completed.

*Step 4:* Component-property Calculation
    In this phase, the component properties are extracted from the
    simulations and statistically treated to yield estimates for their
    expected values and statistical errors.

*Step 5:* Surrogate-model Estimation
    In this phase, the estimates for the expected values of the
    component properties, as well as the corresponding uncertainties,
    are computed for all grid points by applying the surrogate models
    to the estimates and statistical errors for the sampled grid
    points.

*Step 6:* Composite-property Calculation
    In this phase, the estimates for the expected values of the
    composite properties, as well as the corresponding uncertainties,
    are computed based on the surrogate-model estimates for the
    component properties.

*Step 7:* Score Calculation
    In this phase, the score function is calculated for all grid
    points.

*Step 8:* Extension Check
    In this phase, it is checked whether the uncertainties in the
    estimates of the composite properties for the sampled grid points
    are all less than their pre-established tolerances. If that is not
    the case, new simulation lengths are attributed to the
    transgressing protocols, and the program flow returns to Step 3
    with all the other simulations flagged as complete.

*Step 9:* Grid Shifting
    In this phase, it is verified whether the parameter-search grid 
    needs to be shifted, and, if that is the case, the shifting is
    carried out, and the program flow returns to Step 2.
    
    If the grid does not need to be shifted, the program ends.



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
