#################
Customization API
#################

The customization API provides functions that one can use to define
custom behavior in the program, whereby specific parts of the workflow
can be given new implementations (e.g., the function used to calculate
the score of a grid point). The user-defined implementation is given a
label, which can be used in the input file to request it.  The
customization API is provided in the module ``gmak.api``.

The functions in the customization API are named with the pattern
``add_SPEC_custom_OBJ``. ``OBJ`` indicates the modified part of the
program workflow. ``SPEC`` is used to particularize some contexts when
``OBJ`` is not sufficient. Roughly speaking, calling such a function
includes a new implementation for the ``OBJ`` part of the workflow.
The label of this new implementation, as well as the actual
implementation, are provided in the function parameters.  This
frequently involves passing user-defined functions as these functions
parameters. These function-parameter functions have signature
requirements that are specific to ``SPEC`` and ``OBJ``.  To properly
document them, dummy implementations are provided in the
``gmak.api_signatures`` module and referred to in this documentation.

Here, these modules are documented in sections separated by context.
Each one includes also a brief explanation of how the corresponding
functions fit in the program workflow. Auxiliary functions and classes
that are involved in the customization API but not part of it are
documented in the last section.

Usage
=====

The module ``gmak.api`` is made available by installing ``gmak``.
However, simply using it in a Python script or module does not
immediately affect the program. To make the modifications take effect,
the script or module should be saved in the working directory (i.e.,
the directory where the program is deployed) with the name
``custom.py``, and the command-line option ``--custom`` must be passed
to the program. It will then be interpreted at runtime, incorporating
the user-defined implementations into the program.

.. note:: We recommend using the `template file`_ provided with the
   program as a starting point to write the ``custom.py`` file.

.. note:: Customization examples are provided in the
   ``gmak.api_examples`` package. They are discussed in the Example
   Applications part of the documentation.

Systems
=======

In ``gmak``, the type of a system determines the construction of its
topology. There are two aspects to this. First, a suitable
representation for the topology must be chosen, so that it is passed
around in the program and used as input for the simulation and
property-calculation routines. This representation will be referred to
as the topology-output object. Secondly, the way to apply the values
of the parameters modified by the program to these topology-output
objects must be specified.

To add a custom system type to the program, one should use the
function :py:func:`~gmak.api.add_custom_system`. It receives as
parameters the name of the custom system type and two user-implemented
functions which are chained together. The first one creates a
topology-output object for some given system, grid-shift iteration and
grid point. The second one applies the values of the parameters
modified by the program to the topology-output object returned by the
first one.

Alternatively, one can customize the behavior of the default ``gmx``
system type, used for GROMACS-based systems. While nonbonded and
macro-type parameters (see :ref:`overview/interaction_parameters`) are
handled automatically by the program, other types of parameters
(falling into the umbrella category of custom parameters) are not. For
those, the user can
:py:func:`~gmak.api.add_gmx_custom_parameter_writer`.

.. autofunction:: gmak.api.add_custom_system

.. autofunction:: gmak.api.add_gmx_custom_parameter_writer

.. autofunction:: gmak.api_signatures.topo_out_creator

.. autofunction:: gmak.api_signatures.topo_out_writer

.. autofunction:: gmak.api_signatures.gmx_custom_parameter_writer

Protocols
=========

A protocol type determines how the :doc:`protocol simulations
</overview/protocols>` are carried out for a given system.  This
includes the choice of the software package used to carry out the
simulations, the deployment of replicas, the sequence of simulations,
etc.  To add a custom protocol type to the program, one should use the
function :py:func:`~gmak.api.add_custom_protocol`.  Support for
:ref:`extending the simulations <overview/protocols:simulation
extensions>` is also provided, as well as for setting up a
:ref:`followable protocol<overview/protocols:following protocols>`.

Up to four functions can be implemented by the user and supplied to
the protocol-customization function. The first and mandatory one is
the :py:func:`~gmak.api_signatures.simulator` function, which is
responsible for carrying out the simulations. If extending the
simulations is desired, this function has to handle not only setting
up the initial simulation, but also extending it when appropriate.
These two situations can be distinguished based on the value of the
``ext`` parameter, which evaluates to ``True`` for extension
simulations.

The next two functions are optional and should be specified only when
extending the simulations is desired.  The first one of those is
:py:func:`~gmak.api_signatures.calc_initial_len`, which returns the
length of the initial simulation (prior to extending it). It is a
function, and not a variable, to allow for recycling the protocol type
for different applications. The second one is
:py:func:`~gmak.api_signatures.calc_extend`, which contains the logic
for determining the length of an extended simulation.

Finally, :py:func:`~gmak.api_signatures.get_last_frame`, if provided,
makes the protocol type followable. This function returns the
configuration file that the user wants to use as a starting point for
the following protocol. Typically, this will be the last frame of the
production run of the protocol.

.. autofunction:: gmak.api.add_custom_protocol

.. autofunction:: gmak.api_signatures.simulator

.. autofunction:: gmak.api_signatures.calc_initial_len

.. autofunction:: gmak.api_signatures.calc_extend

.. autofunction:: gmak.api_signatures.get_last_frame


Properties
==========

The customization of properties occurs on the levels of component and
composite properties (see :ref:`overview/properties:properties`). For
this reason, to set up custom properties, one must actually set up
custom composite and component properties.  For the trivial case where
the composite property is associated with a single component property,
the API provides a way to automatically derive the former from the
latter one.

To add a custom component-property type to the program, one should use
the function :py:func:`~gmak.api.add_custom_component_property`. It
receives as arguments a calculator function and a boolean parameter
that indicates whether the latter function returns (a) samples of the
property for each simulation frame; or (b) the expected value and
uncertainty of the property. In the former case, the program
automatically takes care of carrying out the statistical treatment of
the sampled values in order to produce the expected value and
uncertainty. In the latter case, this is the responsibility of the
user implementation. The calculator function is called individually
for each sampled grid point and protocol after the corresponding
simulations are performed.

To add a custom composite-property type to the program, one should use
the function :py:func:`~gmak.api.add_custom_composite_property`. It
receives as argument the expected values and uncertainties of the
corresponding component properties for a particular system and grid
point and returns the corresponding expected value and
uncertainty of the composite property.

.. note:: There is no need to define a custom property type for a
   component property that can be extracted with the ``gmx energy``
   program.  A component property of type ``gmx_PROPNAME`` in the
   input file is automatically interpreted as one that is obtained
   with ``gmx energy`` with ``PROPNAME`` as the input string.

.. autofunction:: gmak.api.add_custom_component_property

.. autofunction:: gmak.api.add_custom_composite_property

.. autofunction:: gmak.api_signatures.component_calculator

.. autofunction:: gmak.api_signatures.composite_calculator


Surrogate Model
===============

The surrogate-model type defines the function used to estimate the
expected values and uncertainties of a *component* property for *all*
grid points based on these quantities for the sampled grid points.

To add a custom surrogate-model type to the program, one should use
the function :py:func:`~gmak.api.add_custom_surrogate_model`.

.. autofunction:: gmak.api.add_custom_surrogate_model

.. autofunction:: gmak.api_signatures.compute


Score Function
==============

The score-function type defines how the estimates of the values and
uncertainties of the *composite* properties are used to calculate a
score for the grid point.

To add a custom score-function type, one should use the function
:py:func:`~gmak.api.add_custom_score`.

.. autofunction:: gmak.api.add_custom_score

.. autofunction:: gmak.api_signatures.calc_score

.. autofunction:: gmak.api_signatures.calc_score_err

Grid Shifting
=============

The grid-shifting type defines how the new origin of the grid is calculated 
for grid-shifting procedure.

To add a custom grid-shifting type, one should use the function
:py:func:`~gmak.api.add_custom_gridshifter`.

.. autofunction:: gmak.api.add_custom_gridshifter

.. autofunction:: gmak.api_signatures.calculator


Utilitary Functions and Classes
===============================

.. autofunction:: gmak.cartesiangrid.flat2tuple

.. autoclass:: gmak.custom_attributes.CustomizableAttributesMixin.InputParameters
   :members: 

.. autoclass:: gmak.systems.TopologyOutput
   :members:

.. autoclass:: gmak.gmx_system.GmxTopologyOutput
   :members:
   :show-inheritance:

.. autoclass:: gmak.interaction_parameter.InteractionParameter()
   :members:

.. autoclass:: gmak.interaction_parameter.InteractionParameterType
   :show-inheritance:
   :undoc-members:
   :inherited-members:
   :member-order: bysource

.. autoclass:: gmak.interaction_parameter.InteractionAtom
   :members:

.. autoclass:: gmak.interaction_parameter.InteractionPair
   :members:

.. _template file:
   https://github.com/mssm-labmmol/gridmaker/main/share/custom.py
