from enum import Enum
from typing import Optional
import re


class InteractionParameterType(Enum):
    """
    An enumeration that identifies the type of interaction parameter.
    """

    MacroParameter = 1
    """
    The type assigned to parameters of which the name starts with the ``@``
    character.  For the GROMACS-based system type, it identifies parameters
    that are applied to the topology file via simple textual replacement. For
    custom system types, this has no special behavior.
    """

    LJ_V = 2
    """
    The type assigned to :math:`\sigma` or :math:`C_6` Lennard-Jones
    parameters (see :ref:`overview/interaction_parameters:lennard-jones
    parameters` for naming conventions).
    """

    LJ_W = 3
    """
    The type assigned to :math:`\epsilon` or :math:`C_{12}` Lennard-Jones
    parameters  (see :ref:`overview/interaction_parameters:lennard-jones
    parameters` for naming conventions).
    """

    LJ_14_V = 4
    """
    The type assigned to :math:`\sigma` or :math:`C_6` Lennard-Jones parameters
    when applied to 1,4 (third-neighbor) pairs (see
    :ref:`overview/interaction_parameters:lennard-jones parameters` for naming
    conventions).
    """

    LJ_14_W = 5
    """
    The type assigned to :math:`\epsilon` or :math:`C_{12}` Lennard-Jones
    parameters when applied to 1,4 (third-neighbor) pairs (see
    :ref:`overview/interaction_parameters:lennard-jones parameters` for naming
    conventions).
    """

    CustomParameter = 6
    """
    The type assigned to parameters which are not of the types above.  For the
    GROMACS-based system type, parameters with this type are not handled
    automatically and require that the user
    :py:func:`~gmak.gmx_system.add_gmx_custom_parameter_writer`.
    """


class InteractionParticles:
    pass


class InteractionAtom(InteractionParticles):
    """
    The class representing an interaction atom.
    """
    def __init__(self, name, pairs_include=None, pairs_exclude=None):
        """
        :param name: the name of the atom
        :type name: str
        :param pairs_include: regular expression controlling the standard pairtypes affected by this atom
        :param pairs_include: str
        :param pairs_include: regular expression controlling the standard pairtypes which are not affected by this atom
        :param pairs_include: str
        """
        self.name = name
        if pairs_include is None:
            self.pairs_include = r'.*'
        else:
            self.pairs_include = pairs_include
        if pairs_exclude is None:
            self.pairs_exclude = r'^$'
        else:
            self.pairs_exclude = pairs_exclude

    def __eq__(self, other):
        return self.name == other.name

    @property
    def name(self):
        """
        The name of the atom type. It is set automatically from the name of the
        corresponding :py:class:`InteractionParameter`.

        :type: str
        """
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def pairs_include(self):
        """
        The regular expression controlling the standard pairtypes affected by the atom.

        :type: str
        """
        return self._pairs_include

    @pairs_include.setter
    def pairs_include(self, s):
        self._pairs_include = s

    @property
    def pairs_exclude(self):
        """
        The regular expression controlling the standard pairtypes that are not affected by the atom.

        :type: str
        """
        return self._pairs_exclude

    @pairs_exclude.setter
    def pairs_exclude(self, s):
        self._pairs_exclude = s


class InteractionPair(InteractionParticles):
    """
    The class representing an interaction pair.

    Two pairs can be compared for equality regardless of the order of their
    atom types::

        >>> InteractionPair("CH3", "CH1") == InteractionPair("CH1", "CH3")
        True
        >>> InteractionPair("CH3", "CH1") == InteractionPair("OA", "CH3")
        False

    """
    def __init__(self, name_i, name_j):
        self.name_i = name_i
        self.name_j = name_j

    @property
    def name_i(self):
        """
        The name of the first atom type of the pair. It is set automatically
        from the name of the corresponding :py:class:`InteractionParameter`.

        :type: str
        """
        return self._name_i

    @name_i.setter
    def name_i(self, n):
        self._name_i = n

    @property
    def name_j(self):
        """
        The name of the second atom type of the pair. It is set automatically
        from the name of the corresponding :py:class:`InteractionParameter`.

        :type: str
        """
        return self._name_j

    @name_j.setter
    def name_j(self, n):
        self._name_j = n

    def __eq__(self, other):
        my_at_i = self.name_i
        my_at_j = self.name_j
        at_i = other.name_i
        at_j = other.name_j
        if (my_at_i == at_i) and (my_at_j == at_j):
            return True
        if (my_at_j == at_i) and (my_at_i == at_j):
            return True
        return False

    def to_tuple(self):
        """
        Returns the component atom types as a tuple.

        :rtype: tuple

        .. code-block:: python

            >>> ip = InteractionPair("CH3", "CH1")
            >>> ip.to_tuple()
            ("CH3", "CH1")
            >>> "CH1" in ip.to_tuple()
            True
            >>> "CH3" in ip.to_tuple()
            True
            >>> "OW" in ip.to_tuple()
            False

        """
        return (self.name_i, self.name_j)

    def derives_from_atom(self, atom):
        """
        Verifies if the pair derives from an atom based on the properties
        :py:attr:`~InteractionAtom.pairs_include` and
        :py:attr:`~InteractionAtom.pairs_exclude` of the latter.

        :param atom: The atom that is verified if the pair derives from
        :type atom: :py:class:`~InteractionAtom`
        :return: True if the the pair derives from the atom; False otherwise.
        :rtype: bool
        """
        incl = atom.pairs_include
        excl = atom.pairs_exclude
        other = self.get_other(atom)
        if other is None:
            return False
        if re.match(incl, other) and not re.match(excl, other):
            return True
        return False

    def get_other(self, atom: InteractionAtom) -> InteractionAtom:
        try:
            if self.to_tuple().index(atom.name) == 0:
                other = self.name_j
            else:
                other = self.name_i
        except ValueError:
            return None
        return other


class InteractionParameter:
    """
    The class representing an interaction parameter that is set by the program.
    It contains its type and value and, when adequate, the types of the
    particles involved in the interaction.
    """

    def __init__(self,
                 name: str,
                 type: InteractionParameterType,
                 particles: Optional[InteractionParticles] = None,
                 value: Optional[float] = None):
        self.name = name
        self.type = type
        self.particles = particles
        self.value = value

    @property
    def name(self):
        """
        The string that identifies the interaction parameter. It is the same
        one used in the ``$variation`` block of the input file.

        :type: str
        """
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def type(self):
        """
        The type of interaction parameter

        :type: :py:class:`InteractionParameterType`
        """
        return self._type

    @type.setter
    def type(self, t):
        self._type = t

    @property
    def value(self):
        """
        The value of the interaction parameter

        :type: float
        """
        return self._value

    @value.setter
    def value(self, v):
        self._value = v

    @property
    def particles(self):
        """
        The particles involved in the interaction

        :type: :py:class:`None`, :py:class:`InteractionAtom` or
            :py:class:`InteractionPair`
        """
        return self._particles

    @particles.setter
    def particles(self, p):
        self._particles = p

    @classmethod
    def from_string(cls, string, exclude_pairs=None, include_pairs=None):
        m = re.match(r'^@(\S+)', string)
        if m:
            return cls(m.group(0),
                       InteractionParameterType.MacroParameter)
        m = re.match(r'^(14){,1}_{,1}([VW])_(\S+)_(\S+)', string)
        if m:
            if m.group(1):
                if m.group(2) == 'V':
                    return cls(m.group(0),
                               InteractionParameterType.LJ_14_V,
                               particles=InteractionPair(m.group(3),
                                                         m.group(4)))
                else:
                    return cls(m.group(0),
                               InteractionParameterType.LJ_14_W,
                               particles=InteractionPair(m.group(3),
                                                         m.group(4)))
            else:
                if m.group(2) == 'V':
                    return cls(m.group(0),
                               InteractionParameterType.LJ_V,
                               particles=InteractionPair(m.group(3),
                                                         m.group(4)))
                else:
                    return cls(m.group(0),
                               InteractionParameterType.LJ_W,
                               particles=InteractionPair(m.group(3),
                                                         m.group(4)))

        m = re.match(r'^(14){,1}_{,1}([VW])_(\S+)$', string)
        if m:
            if m.group(1):
                if m.group(2) == 'V':
                    return cls(m.group(0),
                               InteractionParameterType.LJ_14_V,
                               particles=InteractionAtom(m.group(3),
                                                         pairs_include=include_pairs,
                                                         pairs_exclude=exclude_pairs))
                else:
                    return cls(m.group(0),
                               InteractionParameterType.LJ_14_W,
                               particles=InteractionAtom(m.group(3),
                                                         pairs_include=include_pairs,
                                                         pairs_exclude=exclude_pairs))
            else:
                if m.group(2) == 'V':
                    return cls(m.group(0),
                               InteractionParameterType.LJ_V,
                               particles=InteractionAtom(m.group(3),
                                                         pairs_include=include_pairs,
                                                         pairs_exclude=exclude_pairs))
                else:
                    return cls(m.group(0),
                               InteractionParameterType.LJ_W,
                               particles=InteractionAtom(m.group(3),
                                                         pairs_include=include_pairs,
                                                         pairs_exclude=exclude_pairs))

        # If no match until here, it must be a custom parameter.
        return cls(string, InteractionParameterType.CustomParameter)


class CombinationRuleType(Enum):
    C6_C12_Geometric = 1
    Sigma_Epsilon_Mixed = 2
    Sigma_Epsilon_Geometric = 3


class CombinationRule:
    def __init__(self, combtype: int):
        self.combtype = CombinationRuleType(combtype)

    def rescale_param_value(self,
                            new_atomic_param: InteractionParameter,
                            old_pair_param_value: float,
                            old_atom_param_value: float) -> float:
        if self.combtype == CombinationRuleType.C6_C12_Geometric:
            return old_pair_param_value * (new_atomic_param.value/old_atom_param_value) ** .5
        elif self.combtype == CombinationRuleType.Sigma_Epsilon_Mixed:
            if new_atomic_param.type == InteractionParameterType.LJ_V:
                return old_pair_param_value + 0.5 * (new_atomic_param.value - old_atom_param_value)
            elif new_atomic_param == InteractionParameterType.LJ_W:
                return old_pair_param_value * (new_atomic_param.value/old_atom_param_value) ** .5
            else:
                raise ValueError(f"Can't apply combination rule to {new_atomic_param.type}.")
        elif self.combtype == CombinationRuleType.Sigma_Epsilon_Geometric:
            return old_pair_param_value * (new_atomic_param.value/old_atom_param_value) ** .5
        else:
            raise ValueError(f"Unknown combination rule {self.combtype}.")
