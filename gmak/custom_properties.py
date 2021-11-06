from gmak.atomic_properties import add_custom_property

def calculator(topology, protocol_output, property_attrs):
    """
    The function used to calculate the custom property.

    :param topology: The topology file(s) considered in the simulations
    :type topology: str or list of str
    :param protocol_output: The output files of the simulation. For custom
        protocols, it is the ``dict`` returned by the
        :py:func:`~gmak.custom_protocols.simulator` function. For GROMACS-based
        protocols, it is the simulation-output ``dict`` described in
        :ref:`protocols:gromacs-based protocols`.
    :type protocol_output: dict
    :param property_attrs: The property attributes defined in the input file
    :type property_attrs:
        :py:class:`~CustomizableAttributesMixin.CustomizableAttributesData`
    :return: A tuple ``(EA, dEA)`` with the expected value and uncertainty of
        the property, or a list with the values of the property for each frame of
        the simulation. In the latter case, the program will automatically
        calculate the corresponding average and statistical uncertainty.
    :rtype: tuple or list
    """
    pass

