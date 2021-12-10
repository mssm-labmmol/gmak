
#######################
The ``$optimize`` block
#######################

The ``$optimize`` block sets up the composite properties that enter in the calculation of the score function, as well as their weights and uncertainty tolerances.
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
     - *(Mergeable)* *(Optional)* The type of the score function (for the default type, see :doc:`/overview/score`).
     - 
   * - properties
     - List[object reference]
     -  The names of the composite properties that enter in the score.
     - Depends on the properties to which it refers. 
   * - references
     - List[numerical]
     - *(Mergeable)*  The reference/target values of the properties in the score.
     - The value for each property should be compatible with the order in which it is listed in the ``properties`` option. 
   * - weights
     - List[numerical]
     - *(Mergeable)*  The weights of the properties in the score.
     - The value for each property should be compatible with the order in which it is listed in the ``properties`` option. 
   * - tolerances
     - List[numerical]
     - *(Mergeable)*  The uncertainty tolerances of the properties.
     - The value for each property should be compatible with the order in which it is listed in the ``properties`` option. 

.. note:: Parameters that are not listed above can also be supplied.
   They are not recognized by the program in any special way, but are
   parsed and made available in the :doc:`/usage/customization_api`,
   together with all the other block parameters, as an
   :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
   object.

Example
=======

.. code-block:: gmi

    $optimize
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
