
########################
The ``$gridshift`` block
########################

The ``$gridshift`` block sets up the routine for the grid-shifting procedure.
This block can appear only once in the input file.


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

   * - type
     - string
     - *(Optional)* The type of grid-shifting procedure (default value is ``default`` - see :doc:`/overview/grid_shifting`).
     - 
   * - maxshifts
     - numerical
     - *(Mergeable)*  The maximum number of grid-shifting iterations.
     - Setting it to zero disables grid-shifting. 

Type ``default``
================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - margins
     - List[numerical]
     - *(Mergeable)*  The margins :math:`\delta_i` and :math:`\Delta_i` for each grid dimension :math:`i`. The values should be listed in the order :math:`\delta_1`, :math:`\Delta_1`, :math:`\delta_2`, :math:`\Delta_2`, :math:`\cdots`, :math:`\delta_D`, :math:`\Delta_D`.
     - 
   * - ncut
     - numerical
     - *(Mergeable)*  The value of :math:`n_\text{cut}`.
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