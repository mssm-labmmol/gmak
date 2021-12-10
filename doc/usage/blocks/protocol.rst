
#######################
The ``$protocol`` block
#######################

The ``$protocol`` block sets up a :doc:`protocol </overview/protocols>`.
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
     -  The name for the protocol.
     - 
   * - type
     - string
     -  The type of the protocol (see :doc:`/overview/protocols`).
     - 
   * - system
     - object reference
     -  The name of the system with which the protocol is associated.
     - Depends on the system to which it refers. 
   * - coords
     - object reference
     -  The name of the coordinates with which the protocol is associated.
     - Depends on the coordinates to which it refers. 

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

   * - mdps
     - List[path]
     -  The list of file paths of the input-parameter files (extension .mdp) for each step of the simulation workflow.
     - 
   * - maxsteps
     - numerical
     - *(Mergeable)*  The maximum number of steps allowed for the length of the production run.
     - 
   * - minfactor
     - numerical
     - *(Mergeable)* *(Optional)* Controls the minimum value allowed for the updated length of a simulation that has been extended (see :ref:`simulation extensions <overview/protocols:simulation extensions>`; default is 1.1).
     - 

Type ``gmx_alchemical``
=======================

 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks

   * - mdps
     - List[path]
     -  The list of file paths of the template input-parameter files (extension .mdp) for each step of the simulation workflow.
     - 
   * - maxsteps
     - numerical
     - *(Mergeable)*  The maximum number of steps allowed for the length of the production runs.
     - 
   * - minfactor
     - numerical
     - *(Mergeable)* *(Optional)* Controls the minimum value allowed for the updated length of a simulation that has been extended (see :ref:`simulation extensions <extensions2>`; default is 1.1).
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

    $protocol
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
