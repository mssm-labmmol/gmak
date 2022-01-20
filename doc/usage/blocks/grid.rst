
###################
The ``$grid`` block
###################

The ``$grid`` block sets up the :ref:`sampled points <overview/grid:the list of sampled points>`.
This block can appear only once in the input file.
It intrinsically depends on the :doc:`$variation <variation>`, :doc:`$gridshift <gridshift>` blocks.

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

   * - samples
     - List[numerical]
     -  The linear indexes of the sampled points of the grid.
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