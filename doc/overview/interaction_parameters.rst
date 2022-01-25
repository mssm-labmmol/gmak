######################
Interaction Parameters
######################

In the context of ``gmak``, force-field parameters are represented as
interaction-parameter objects. An interaction parameter is initialized
from a string (its *name*) that may identify the type of the
interaction in the context of the force-field and also the particles
involved. The name of the interaction parameter is also used to
unambiguously associate it to one of :ref:`the
parameter-search-grid<overview/grid:the parameter-search grid>` axes.
The value of the interaction parameter is set based on this
association for each grid point.

The following interaction-parameter types are available in ``gmak``.

.. seealso::

   :py:class:`~gmak.interaction_parameter.InteractionParameter`
       The class of the interaction-parameter object.


Lennard-Jones parameters
========================


Atomtype Parameters
-------------------

Atomtype parameters are identified by a name with pattern
``PARAMETER_ATOMTYPE``, where ``ATOMTYPE`` is a string identifying the
type of the atom within the forcefield, and ``PARAMETER`` is either
the character ``V`` or the character ``W``.

If the Lennard-Jones potential is represented in the
:math:`C_6`-:math:`C_{12}` form, ``V`` is interpreted as :math:`C_6`
and ``W``, as :math:`C_{12}`. If it is in the
:math:`\sigma`-:math:`\epsilon` form, ``V`` is interpreted as
:math:`\sigma`, and ``W``, as :math:`\epsilon`.

The value of an atomtype parameter directly affects that of the
interaction between self pairs, i.e. those involving two atoms of the
given atomtype. The effect over standard (i.e., non-1,4) interaction
pairs involving one atom of the given atomtype (referred to as the
*reference* atomtype) and a second atom of a different atomtype
(referred to as the *other* atomtype) is controlled by two regular
expressions:

pairs_include
    The regular expression controlling which standard pairtypes may
    have parameters affected by those of the reference atomtype. The
    regular expression itself matches with the name of the other
    atomtype. The default value is ``.*`` (matches other atomtypes
    with any name).

pairs_exclude
    The regular expression controlling which standard pairtypes do
    *not* have their parameters affected by those of the reference
    atomtype. The regular expression itself matches with the name of
    the other atomtype. The default value is ``^$`` (matches
    nothing).

Only the standard pairtypes that match against the regular expression
``pairs_include`` but not against ``pairs_exclude`` are affected by
the reference atomtype.

.. seealso::
   :py:class:`~gmak.interaction_parameter.InteractionParameterType.LJ_V`, :py:class:`~gmak.interaction_parameter.InteractionParameterType.LJ_W`
       Members of the enumeration :py:class:`~gmak.interaction_parameter.InteractionParameterType` 
       that identify standard Lennard-Jones parameters in ``gmak``.

.. seealso::
   :py:class:`~gmak.interaction_parameter.InteractionAtom`
       The class representing interaction atoms.

.. seealso::
   :py:meth:`~gmak.interaction_parameter.InteractionPair.derives_from_atom`
       Method of the
       :py:class:`~gmak.interaction_parameter.InteractionPair` class used to
       verify if a pair derives from an atom following the
       regular-expression rules described above.
    

Standard-pairtype Parameters
----------------------------

Standard-pairtype parameters are identified by a name with pattern
``PARAMETER_ATOMTYPEI_ATOMTYPEJ``, where ``ATOMTYPEI`` and
``ATOMTYPEJ`` are strings identifying the types of the atoms involved
in the interaction within the force field, and ``PARAMETER`` is either
the character ``V`` or the character ``W``, as described above.

.. seealso::
   :py:class:`~gmak.interaction_parameter.InteractionParameterType.LJ_V`, :py:class:`~gmak.interaction_parameter.InteractionParameterType.LJ_W`
       Members of the enumeration :py:class:`~gmak.interaction_parameter.InteractionParameterType` 
       that identify standard Lennard-Jones parameters in ``gmak``.

.. seealso::
   :py:class:`~gmak.interaction_parameter.InteractionPair`
       The class representing interaction pairs.


1,4-pairtype Parameters
-----------------------

1,4-pairtype parameters are identified by a name with pattern
``14_PARAMETER_ATOMTYPEI_ATOMTYPEJ``, where ``ATOMTYPEI`` and
``ATOMTYPEJ`` are strings identifying the types of the atoms involved
in the interaction within the force field, and ``PARAMETER`` is either
the character ``V`` or the character ``W``, as described above.

The 1,4-pairtype parameters differ from the standard-pairtype
parameters in that they are applied exclusively to particles in the
same molecule and separated by exactly three consecutive bonds.

.. seealso::
   :py:class:`~gmak.interaction_parameter.InteractionParameterType.LJ_14_V`, :py:class:`~gmak.interaction_parameter.InteractionParameterType.LJ_14_W`
       Members of the enumeration :py:class:`~gmak.interaction_parameter.InteractionParameterType` 
       that identify 1,4 Lennard-Jones parameters in ``gmak``.

.. seealso::
   :py:class:`~gmak.interaction_parameter.InteractionPair`
       The class representing interaction pairs.


Macro-based parameters
======================

Macro-based parameters are identified by a name that starts with the
``@`` character.

.. seealso::
   
   :py:class:`~gmak.interaction_parameter.InteractionParameterType.MacroParameter`
       Member of the enumeration :py:class:`~gmak.interaction_parameter.InteractionParameterType` 
       that identifies macro-based parameters in ``gmak``.


Custom parameters
=================

Custom parameters are those whose name does not match any of the
patterns above. They are an umbrella type to accomodate parameters
that are not processed in any special way by the program, leaving this
job to the user.

.. seealso::
   
   :py:class:`~gmak.interaction_parameter.InteractionParameterType.CustomParameter`
       Member of the enumeration :py:class:`~gmak.interaction_parameter.InteractionParameterType` 
       that identifies custom parameters in ``gmak``.

