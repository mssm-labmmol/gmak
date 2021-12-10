######################
Systems and Topologies
######################

In ``gmak``, defining a system type means specifying how to 
write a topology file, given

- the interaction-parameter objects (see
  :doc:`/overview/interaction_parameters`)
- additional data given in the input file in the form of
  :py:class:`input
  parameters<gmak.custom_attributes.CustomizableAttributesMixin.InputParameters>`.

``gmak`` provides a system type suitable for use with the GROMACS
software package---the ``gmx`` system type, described below.

For flexibility, ``gmak`` also provides a function in the 
:ref:`customization API<usage/customization_api:systems>` that allows
one to define a custom system type---either suitable for GROMACS or
for other simulation packages.

GROMACS-compatible Systems
==========================

Template Topology
-----------------

The ``gmx`` system type requires the path of a topology file with
extension ``.top`` to be provided as an input parameter named
``template``. This file works as a template, where the values of the
interaction parameters modified by ``gmak`` are replaced for at each
grid point.

.. note:: When the template topology is read, all ``#include``
   directives are expanded and replaced by the contents of the
   corresponding file.  If the directive does not point to an absolute
   path, ``gmak`` looks for it recursively in the local directory and
   in the GROMACS library directory as inferred from the path of the
   ``gmx`` binary or the ``GMXLIB`` environment variable.


Replacement of Interaction Parameters
-------------------------------------

The replacement of the values of the interaction parameters in a
template is carried out recursively. This means that each interaction
parameter :math:`I_n` is replaced for in a template :math:`T_{n-1}`
that results from replacing the previous parameter :math:`I_{n-1}` on
a previous template :math:`T_{n-2}`, and so on. The starting point
:math:`T_0` is, of course, the template file provided as an input
parameter in the
input file. Schematically:

.. math::
    T_0 \xrightarrow{I_1} T_1 \xrightarrow{I_2} T_2 \xrightarrow{I_3}
    \ldots T_{n-1} \xrightarrow{I_n} T_n


The handling of the interaction parameters for the ``gmx`` system type
is described below.


Lennard-Jones Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

Atomtypes
"""""""""

The values of :ref:`atomtype interaction
parameters<overview/interaction_parameters:atomtype parameters>` are
replaced for in a template by substituing them for the values
of the corresponding atomtype parameters inside the ``[ atomtypes ]``
directive. The original template values are stored for further use.

After that, the standard pairs inside the ``[ nonbond_params ]``
directives are scanned and those that derive from the atomtype based
on the inclusion/exclusion regular-expression rules explained in
:ref:`overview/interaction_parameters:atomtype parameters` have their
interaction-parameter values scaled to reflect the change in the
parameters of the atomtype.

The value of the scaling factor depends on 

- the Lennard-Jones potential representation
  (:math:`C_6`-:math:`C_{12}` or :math:`\sigma`-:math:`\epsilon`)
- the combination rule
- the type of the parameter in the context of the Lennard-Jones
  potential (:math:`C_6/\sigma` or :math:`C_{12}/\epsilon`) 

In GROMACS, the Lennard-Jones representation and the combination rule
are `tied
together <https://manual.gromacs.org/current/reference-manual/topologies/parameter-files.html#non-bonded-parameters>`_
and chosen in the ``[ defaults ]`` directive of the topology. Using
:math:`\alpha` to denote the type of the interaction parameter being
adjusted in the context of the Lennard-Jones potential, and
:math:`\alpha_i` and :math:`\alpha_f` to denote the template (in the
recursive sense described above) and target values, respectively, the
scaling factor :math:`f` is given by

.. math::
   f = 1 + \frac{1}{2}\left(\frac{\alpha_f}{\alpha_i} - 1 \right)

for :math:`\alpha \equiv \sigma` and combination rule 2, and

.. math::
   f = \sqrt{\alpha_f/\alpha_i}

for all other cases.

.. note::
   We highly encourage reading the `parameter-files section
   <https://manual.gromacs.org/current/reference-manual/topologies/parameter-files.html>`_
   of the GROMACS user guide for more information on the combination
   rules.

Standard Pairtypes
""""""""""""""""""

The values of :ref:`standard-pairtype interaction
parameters<overview/interaction_parameters:standard-pairtype
parameters>` are replaced for in a template by substituting
them for the values of the corresponding standard-pairtype parameters
inside the ``[ nonbond_params ]`` directive.

1,4 Pairtypes
"""""""""""""

The values of :ref:`1,4-pairtype interaction
parameters<overview/interaction_parameters:1,4-pairtype parameters>`
are replaced for in a template by substituting them for the
values of the corresponding 1,4-pairtype parameters inside the ``[
pairtypes ]`` directive.


.. warning::
   
   Under no circumstances will ``gmak`` *add* a parameter type to the
   template file---it only *modifies* the values of existing parameter
   types.


Macro-based Parameters
~~~~~~~~~~~~~~~~~~~~~~

The values of :ref:`macro-based interaction
parameters<overview/interaction_parameters:macro-based parameters>`
are replaced for a template by replacing every instance of
the interaction-parameter name (a string token) by the
interaction-parameter value.

Custom Parameters
~~~~~~~~~~~~~~~~~

The values of :ref:`custom interaction
parameters<overview/interaction_parameters:custom parameters>` are
replaced for in a template by relying on a user-implemented
function (see :ref:`usage/customization_api:systems` in the
customization API).

