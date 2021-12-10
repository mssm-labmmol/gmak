
##########################
The ``$coordinates`` block
##########################

The ``$coordinates`` block sets up a routine for constructing initial configurations.
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
     -  The name chosen for the configuration-initialization routine.
     - 
   * - type
     - string
     -  The type of the configuration-initialization routine (see :doc:`/overview/coordinates`).
     - 

Type ``any``
============

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - coords
     - path
     -  The path of the recycled configuration file.
     - 

Type ``folllow``
================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - follow
     - object reference
     -  The name of the protocol from which the recycled configuration is retrieved.
     - Depends on the protocol to which it refers and on the :doc:`$grid <grid>` block. 

Type ``gmx_liquid``
===================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - coords
     - path
     -  The path of a configuration file containing a single copy of the desired molecule.
     - 
   * - nmols
     - numerical
     -  The desired number of molecules.
     - 
   * - box
     - List[string, numerical] OR List[string, numerical, numerical, numerical]
     -  The type of box (``cubic`` or ``rectangular``) followed by the desired box dimensions (one length or a triple of lengths).
     - 

Type ``gmx_solvation``
======================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - coords
     - List[path, path]
     -  A pair of configuration-file paths for the solute and solvent molecules.
     - 
   * - nmols
     - List[numerical, numerical]
     -  A pair of integers specifying the desired number of solute and solvent molecules, respectively.
     - 
   * - box
     - List[string, numerical] OR List[string, numerical, numerical, numerical]
     -  The type of box (``cubic`` or ``rectangular``) followed by the desired box dimensions (one length or a triple of lengths).
     - 

Type ``gmx_slab``
=================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - coords
     - path
     -  The path of a configuration file containing a single copy of the desired molecule.
     - 
   * - nmols
     - numerical
     -  The desired number of molecules.
     - 
   * - box
     - List[string, numerical] OR List[string, numerical, numerical, numerical]
     -  The type of box (``cubic`` or ``rectangular``) followed by the desired box dimensions (one length or a triple of lengths).
     - 
   * - axis
     - string
     - *(Optional)* The axis extended in preparing the simulation box. Possible values are x, y or z (default).
     - 
   * - factor
     - numerical
     - *(Optional)* The factor by which the length of the extended axis is multiplied (default is 5.0).
     - 

Type ``gmx_slab_follow``
========================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - follow
     - object reference
     -  The name of the protocol from which the recycled pure-liquid configuration is retrieved.
     - Depends on the protocol to which it refers and on the :doc:`$grid <grid>` block. 
   * - factor
     - numerical
     - *(Optional)* The factor by which the length of the extended axis is multiplied (default is 5.0).
     - 

.. note:: Parameters that are not listed above can also be supplied.
   They are not recognized by the program in any special way, but are
   parsed and made available in the :doc:`/usage/customization_api`,
   together with all the other block parameters, as an
   :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
   object.

Example
=======

.. code-block:: gmi

    $coordinates
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
