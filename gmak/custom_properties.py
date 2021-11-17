from gmak.component_properties import add_custom_component_property

def component_calculator(topology, protocol_output, property_attrs):
    """
    The function used to calculate the custom component property.

    :param topology: The topology considered in the simulations
    :type topology: :py:class:`~gmak.systems.TopologyOutput`
    :param protocol_output: The output files of the simulation. For custom
        protocols, it is the ``dict`` returned by the
        :py:func:`~gmak.custom_protocols.simulator` function. For GROMACS-based
        protocols, it is the simulation-output ``dict`` described in
        :ref:`protocols:gromacs-based protocols`.
    :type protocol_output: dict
    :param property_attrs: The property attributes defined in the input file
    :type property_attrs:
        :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.CustomizableAttributesData`
    :return: A tuple ``(EA, dEA)`` with the expected value and uncertainty of
        the property, or a list with the values of the property for each frame of
        the simulation. In the latter case, the program will automatically
        calculate the corresponding average and statistical uncertainty.
    :rtype: tuple or list
    """
    pass

def composite_calculator(values, errs, property_attrs):
    """
    :param values: The list with the expected values of each component
        property. Each member is a tuple ``(PTYPE, VALUE)`` where
        ``PTYPE`` is the type of component property and ``VALUE`` is the
        expected value.
    :type values: list
    :param errs: The list with the uncertainties of each component
        property. Each member is a tuple ``(PTYPE, VALUE)`` where
        ``PTYPE`` is the type of component property and ``VALUE`` is the
        estimated uncertainty.
    :type errs: list
    :param property_attrs: The property attributes defined in the input file
    :type property_attrs:
        :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.CustomizableAttributesData`
    :return: A tuple ``(EA, dEA)`` with the expected value and uncertainty of
        the composite property
    :rtype: tuple
    """
    pass
