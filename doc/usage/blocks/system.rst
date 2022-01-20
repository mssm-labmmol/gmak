
#####################
The ``$system`` block
#####################

The ``$system`` block sets up a routine for constructing the topology of a system.
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
     -  The name for the system/topology-construction routine.
     - 
   * - type
     - string
     -  The type of the system/topology-construction routine (see :doc:`/overview/systems_and_topologies`).
     - 

Type ``gmx``
============

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - template
     - path
     -  The path of the topology file to be used as a template.
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