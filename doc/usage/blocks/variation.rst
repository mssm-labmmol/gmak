
########################
The ``$variation`` block
########################

The ``$variation`` block sets up the variations (see :doc:`/overview/grid`).
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
     -  The name for the variation object.
     - The :ref:`main variation <overview/grid:main variation>` *must* be named ``main``. 
   * - type
     - string
     -  The type of the variation object. The allowed values are cartesian, coupled and explicit.
     - 
   * - pars
     - List[string]
     -  !FIX_MISSING_INCLUDE_EXCLUDE! The names of the parameters associated with the variation (see :doc:`/overview/interaction_parameters`).
     - 

Type ``cartesian``
==================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - function
     - List[expression]
     - *(Optional)* A list of mathematical expressions that correspond to the components of the enhanced-cartesian-variation scaling function (see :ref:`overview/grid:cartesian variation`).
     - Each expression accepts numbers, the usual operations ``+,-,*,/,**`` as well as the functions ``exp`` and ``log``. The original components of a variation member can be referred to as ``x[0]``, ``x[1]``, etc. The expressions cannot contain spaces. 
   * - start
     - List[numerical]
     -  The components of the cartesian-variation origin.
     - 
   * - step
     - List[numerical]
     -  The components of the cartesian-variation step.
     - 
   * - size
     - List[numerical]
     -  The dimensions of the cartesian variation.
     - 

Type ``coupled``
================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - using
     - object reference
     -  The name of the variation to which the variation at hand is coupled.
     - Depends on the variation to which it refers. 
   * - function
     - List[expression]
     -  The list of mathematical expressions that correspond to the components of the coupling function (see :ref:`overview/grid:coupled variation`).
     - Each expression accepts numbers, the usual operations ``+,-,*,/,**`` as well as the functions ``exp`` and ``log``. The original components of a variation member can be referred to as ``x[0]``, ``x[1]``, etc. The expressions cannot contain spaces. 

Type ``explicit``
=================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - dim
     - numerical
     -  The number of dimensions of the grid members, i.e. how many parameters are defined by each grid point.
     - This value must match with the number of parameters specified in the ``pars`` option. 
   * - values
     - List[numerical]
     -  The explicit parameter values associated with each grid member. The values are read in batches of size ``dim`` for each grid point.
     - The resulting grid is uni-dimensional with size equal to the number of entries divided by ``dim``. 

.. note:: Parameters that are not listed above can also be supplied.
   They are not recognized by the program in any special way, but are
   parsed and made available in the :doc:`/usage/customization_api`,
   together with all the other block parameters, as an
   :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
   object.

Example
=======

.. code-block:: gmi

    $variation
    TO_BE_REPLACED_BY_TUTORIAL
    $end


Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam
nonumy eirmod tempor invidunt ut labore et dolore magna aliquyam
erat, sed diam voluptua. At vero eos et accusam et justo duo dolores
et ea rebum.  Stet clita kasd gubergren, no sea takimata sanctus est
Lorem ipsum dolor sit amet.

Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam
nonumy eirmod tempor invidunt ut labore et dolore magna aliquyam
erat, sed diam voluptua. At vero eos et accusam et justo duo dolores
et ea rebum.  Stet clita kasd gubergren, no sea takimata sanctus est
Lorem ipsum dolor sit amet.
