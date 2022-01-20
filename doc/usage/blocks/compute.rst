
######################
The ``$compute`` block
######################

The ``$compute`` block sets up a composite property and the associated component properties (see :ref:`overview/properties:component and composite properties`).
This block can appear multiple times in the input file.


The input parameters that are set in this block are listed below.

Block parameters
================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - name
     - string
     -  The name for the composite property.
     - 
   * - type
     - string
     -  The type of the composite property.
     - 
   * - protocols
     - List[object reference]
     - *(Type-specific)* The names of the protocols on which the calculation of the component properties is based.
     - Depends on the protocols to which it refers. Each composite-property type requires different protocols; see below. 
   * - surrogate_model
     - object reference
     -  The type of surrogate model used in the estimation of the component properties.
     - 
   * - components
     - List[string]
     - *(Required for custom properties only)* The types of the component properties involved in the custom composite property.
     - 

Type ``density``
================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - protocols
     - object reference
     -  The name of the protocol for which the property is calculated.
     - Depends on the protocol to which it refers. Read the protocol requirements in :ref:`overview/properties:density`. 

Type ``gamma``
==============

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - protocols
     - object reference
     -  The name of the protocol for which the property is calculated.
     - Depends on the protocol to which it refers. Read the protocol requirements in :ref:`overview/properties:surface-tension coefficient`. 

Type ``dg``
===========

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - protocols
     - object reference
     -  The name of the protocol for which the property is calculated.
     - Depends on the protocol to which it refers. Read the protocol requirements in :ref:`overview/properties:free-energy difference`. 
   * - temperature
     - numerical
     - *(Optional)* The value of the reference temperature of the alchemical-transformation simulations. By default, it is inferred from the ``.mdp`` file of the production run.
     - 

Type ``dhvap``
==============

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - protocols
     - List[object reference, object reference, object reference]
     -  The names of the protocols for which the liquid-phase potential energy, gas-phase potential energy and polarization-energy correction, respectively, are calculated.
     - Depends on the protocols to which it refers. Read the protocol requirements in :ref:`overview/properties:enthalpy of vaporization`. The calculation of the gas-phase potential energy and/or the polarization-energy correction can be ignored by using ``none`` as the corresponding protocol names. 
   * - nmols
     - numerical
     -  The number of molecules in the liquid phase.
     - 
   * - mu
     - numerical
     - *(Optional)* The experimental value of the molecular dipole moment in the gas phase (in cubic nanometers).
     - Required for polarization-energy correction. 
   * - alpha
     - numerical
     - *(Optional)* The value of the molecular isotropic polarizability in the gas-phase (in Debye).
     - Required for polarization-energy correction. 
   * - C
     - numerical
     - *(Optional)* The value of the additional constant corrections (in kJ per mole; default is 0.0).
     - 
   * - temperature
     - numerical
     - *(Optional)* The value of the reference temperature. By default, it is inferred from the ``.mdp`` file of the production run of the liquid simulation.
     - 

.. note:: Parameters that are not listed above can also be supplied.
   They are not recognized by the program in any special way, but are
   parsed and made available in the :doc:`/usage/customization_api`,
   together with all the other block parameters, as an
   :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
   object.

Example
=======

See :doc:`/examples/tutorial` for a commented example.