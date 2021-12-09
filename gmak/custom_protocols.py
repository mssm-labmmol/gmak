from gmak.protocols import add_custom_protocol

def simulator(length, topology, coords, ext, protocol_attrs, workdir):
    """
    The main simulator function for a given grid point.

    :param length: The length of the current production run
    :type length: int or float
    :param topology: The topology considered in the simulations
    :type topology: :py:class:`~gmak.systems.TopologyOutput`
    :param coords: The initial configuration file for the simulations
    :type coords: str or list of str
    :param ext: Indicates whether this is an extension or the first simulation.
    :type ext: bool
    :param protocol_attrs: The protocol attributes defined in the input file
    :type protocol_attrs: :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    :param workdir: The directory where the simulations are run.
    :type workdir: str
    :return: Output files of the simulation
    :rtype: dict
    """
    pass

def calc_initial_len(protocol_attrs):
    """
    Calculates the initial length of the production runs for a given grid
    point.

    :param protocol_attrs: The protocol attributes defined in the input file
    :type protocol_attrs: :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    :return: The initial length of the production run
    :rtype: int or float
    """
    pass


def calc_extend(errs_tols, last_length, protocol_attrs):
    """
    Returns the new length of the simulation based on the uncertainties and
    tolerances of the properties involving the protocol for the given grid
    point. Return :py:obj:`None` to indicate that no more extensions are
    needed.

    :param errs_tols: The dictionary of uncertainties and tolerances. The keys
        are the property names. The values are dictionaries
        ``{'tol': TOL, 'err': ERR}``, where ``TOL`` (float) is the tolerance and
        ``ERR`` (float) is the uncertainty for the property.
    :type errs_tols: dict
    :param last_length: The current length of the production run.
    :type last_length: int or float
    :param protocol_attrs: The protocol attributes defined in the input file
    :type protocol_attrs: :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    :return: The new length of the production run, or None
    :rtype: int or float
    """
    pass

def get_last_frame(protocol_output):
    """
    Returns the path of the configuration file corresponding to the last frame
    of the production run.

    :param protocol_output: The output files of the simulation (returned by
        :py:func:`~simulator`)
    :type protocol_output: dict
    :return: The path of the configuration file corresponding to the last frame
        of the production run
    :rtype: str
    """
    pass

