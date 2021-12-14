"""
Signatures for the parameter functions of the customization API.
"""

def calculator(tuple_indexes,
               scores,
               propnames,
               averages,
               uncertainties):
    """
    The grid-shifting procedure.

    :param tuple_indexes: The tuple indexes of all grid points, ordered by
        linear index.
    :type tuple_indexes: list
    :param scores: A list with the scores of the grid points, ordered by linear
        index.
    :type scores: list
    :param propnames: The names of the composite properties.
    :type propnames: list
    :param averages: A 2D array with the estimated values of the composite
        properties for all grid points.  The first index is the linear index of
        the grid point.  The second index corresponds to the property index in
        ``propnames``.
    :type averages: numpy.ndarray
    :param uncertainties: A 2D array with the estimated uncertainties of the
        composite properties for all grid points.  The first index is the
        linear index of the grid point.  The second index corresponds to the
        property index in ``propnames``.
    :type uncertainties: numpy.ndarray
    :return: The :ref:`shifting tuple <overview/parameter-search grid:Variation Shifting>`,
        or :py:obj:`None` to complete the run. For example, a return value of
        (-5, 5) indicates that the origin is shifted by -5 grid cells in the
        first coordinate and +5 grid cells in the second.
    :rtype: tuple or None
    """
    pass


def component_calculator(topology, protocol_output, property_pars):
    """
    The function used to calculate the custom component property.

    :param topology: The topology considered in the simulations
    :type topology: :py:class:`~gmak.systems.TopologyOutput`
    :param protocol_output: The output files of the simulation. For custom
        protocols, it is the ``dict`` returned by the
        :py:func:`~gmak.api_signatures.simulator` function. For GROMACS-based
        protocols, it is the simulation-output ``dict`` described in
        :ref:`protocols:gromacs-based protocols`.
    :type protocol_output: dict
    :param property_pars: The property input parameters defined in the input file
    :type property_pars:
        :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    :return: A tuple ``(EA, dEA)`` with the expected value and uncertainty of
        the property, or a list with the values of the property for each frame of
        the simulation. In the latter case, the program will automatically
        calculate the corresponding average and statistical uncertainty.
    :rtype: tuple or list
    """
    pass


def composite_calculator(values, errs, property_pars):
    """
    The function used to calculate the custom composite property.

    :param values: The list with the expected values of each component
        property. Each member is a tuple ``(PTYPE, VALUE)`` where ``PTYPE`` (str)
        is the type of component property and ``VALUE`` (float) is the expected
        value.
    :type values: list
    :param errs: The list with the uncertainties of each component
        property. Each member is a tuple ``(PTYPE, VALUE)`` where ``PTYPE``
        (str) is the type of component property and ``VALUE`` (float) is the
        estimated uncertainty.
    :type errs: list
    :param property_pars: The property input parameters defined in the input file
    :type property_pars:
        :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    :return: A tuple ``(EA, dEA)`` with the expected value and uncertainty of
        the composite property
    :rtype: tuple
    """
    pass


def simulator(length, topology, coords, ext, protocol_pars, workdir):
    """
    The main simulator function.

    :param length: The length of the current production run
    :type length: int or float
    :param topology: The topology considered in the simulations
    :type topology: :py:class:`~gmak.systems.TopologyOutput`
    :param coords: The path of the initial configuration file for the simulations
    :type coords: str
    :param ext: Indicates whether this is an extension or the first simulation.
    :type ext: bool
    :param protocol_pars: The protocol input parameters defined in the input file
    :type protocol_pars: :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    :param workdir: The directory where the simulations are run.
    :type workdir: str
    :return: Output files of the simulation
    :rtype: dict
    """
    pass


def calc_initial_len(protocol_pars):
    """
    Calculates the initial length of the production runs.

    :param protocol_pars: The protocol input parameters defined in the input file
    :type protocol_pars: :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    :return: The initial length of the production run
    :rtype: int or float
    """
    pass


def calc_extend(errs_tols, last_length, protocol_pars):
    """
    Returns the new length of the simulation based on the uncertainties and
    tolerances of the properties involving the protocol.  Return :py:obj:`None`
    to indicate that no more extensions are needed.

    :param errs_tols: The dictionary of uncertainties and tolerances. The keys
        are the property names. The values are dictionaries
        ``{'tol': TOL, 'err': ERR}``, where ``TOL`` (float) is the tolerance and
        ``ERR`` (float) is the uncertainty for the property.
    :type errs_tols: dict
    :param last_length: The current length of the production run.
    :type last_length: int or float
    :param protocol_pars: The protocol input parameters defined in the input file
    :type protocol_pars: :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    :return: The new length of the production run, or None
    :rtype: int or float or None
    """
    pass


def get_last_frame(protocol_output):
    """
    Returns the path of the configuration file to be followed. It typically
    corresponds to the last frame of the production run.

    :param protocol_output: The output files of the simulation (returned by
        :py:func:`~simulator`)
    :type protocol_output: dict
    :return: The path of the configuration file to be followed
    :rtype: str
    """
    pass


def calc_score(estimates, errors, weis, refs):
    """
    Calculates the score for a grid point.

    :param estimates: The dictionary of composite-property estimates for a grid
        point. The keys are the property names. Each key maps to the expected
        value of the property at a given grid point.
    :type estimates: dict
    :param errors: The dictionary of composite-property uncertainties for a grid
        point. The keys are the property names. Each key maps to the value of
        the uncertainty of the property at a given grid point.
    :type errors: dict
    :param weis: The dictionary of composite-property weights for a grid point.
        The keys are the property names. Each key maps to the weight of the
        property in the optimization.
    :type weis: dict
    :param refs: The dictionary of composite-reference values of the properties
        for a grid point. The keys are the property names. Each key maps to the
        reference values of the property in the optimization.
    :type refs: dict
    :return: The value of the score
    :rtype: float
    """
    pass


def calc_score_err(estimates, errors, weis, refs):
    """
    Calculates the score uncertainty for a grid point.

    :param estimates: The dictionary of composite-property estimates for a grid
        point. The keys are the property names. Each key maps to the expected
        value of the property at a given grid point.
    :type estimates: dict
    :param errors: The dictionary of composite-property uncertainties for a grid
        point. The keys are the property names. Each key maps to the value of
        the uncertainty of the property at a given grid point.
    :type errors: dict
    :param weis: The dictionary of composite-property weights for a grid point.
        The keys are the property names. Each key maps to the weight of the
        property in the optimization.
    :type weis: dict
    :param refs: The dictionary of composite-reference values of the properties
        for a grid point. The keys are the property names. Each key maps to the
        reference values of the property in the optimization.
    :type refs: dict
    :return: A tuple ``(MIN, MAX)`` that represents a confidence interval for
        the score
    :rtype: tuple
    """
    pass


def compute(EA_s, dEA_s, I_s, gridshape, X_ki):
    """
    Computes the expected values and uncertainties of a component property for the
    entire grid.

    :param EA_s: A 1D array of shape ``(NSAMP,)`` with the average values of the
        property for each sampled point.
    :type EA_s: numpy.ndarray
    :param dEA_s: A 1D array of shape ``(NSAMP,)`` with the uncertainties of the
        property for each sampled point.
    :type dEA_s: numpy.ndarray
    :param I_s: A list of length ``NSAMP`` with the linear indexes of the
        sampled grid points.
    :type I_s: list
    :param gridshape: A tuple with the grid dimensions.
    :type gridshape: tuple
    :param X_ki: A 2D array with the parameter-space values for each grid point
        (more specifically, the parameters of the :ref:`main variation`). The
        first index is the linear index, and the second index is the coordinate
        of the parameter-space point.  For example, ``X_ki[1,0]`` is the value
        of the first parameter of the :ref:`main variation` for the grid point
        with linear index equal to 1.
    :type X_ki: numpy.ndarray
    :return: A tuple ``(EA_k, dEA_k)`` containing the estimated values and
        uncertainties, respectively, of the property for each grid point.
        ``EA_k`` and ``dEA_k`` are 1D arrays indexed by the linear index of
        the grid point.
    :rtype: tuple


    .. note:: The function :py:func:`gmak.cartesiangrid.flat2tuple` can be
       used to convert a linear to a tuple index if desired.
    """
    pass


def gmx_custom_parameter_writer(param, istream, ostream):
    """
    Example of the signature of a function used to apply the values of custom
    parameters to a GROMACS topology file. Such a function is called for all
    custom parameters, so if any distinction needs to be made between them, it
    should be done in the body of this function.

    :param param: The custom interaction parameter.
    :type param: :py:class:`~gmak.interaction_parameter.InteractionParameter`
    :param istream: A readable input stream with the content of an intermediate
        topology file to which the value of ``param`` is not yet applied. This
        file is obtained from the template file given in the ``$system`` block
        by expanding the ``#include`` directives in it and applying to the
        resulting content the values of the interaction parameters specified
        before ``param`` in the input file.
    :type istream: io.StringIO
    :param ostream: A writable stream in which the user must put the content
        of the input stream with the value of ``param`` properly applied to it.
    :type ostream: io.StringIO
    """
    pass


def topo_out_creator(workdir, name, grid, state, system_pars):
    """
    Returns a :py:class:`~gmak.systems.TopologyOutput` that encodes a topology.

    This instance can be initialized with a simple ``TopologyOutput()``
    statement, and, after that, the user can set the instance attributes as
    desired. The returned instance is the same one passed as a parameter in
    :py:func:`~gmak.api_signatures.topo_out_writer`. It is also passed to the
    :py:func:`~gmak.api_signatures.simulator` and
    :py:func:`~gmak.api_signatures.component_calculator` functions.

    :param workdir: The directory where the topology files should be written
    :type workdir: str
    :param name: The name of the current system
    :type name: str
    :param grid: The current grid-shift iteration (see
        :ref:`overview/general_workflow:general workflow`)
    :type grid: int
    :param state: The linear index of the current grid point (see
        :ref:`overview/grid_variations:indexing`)
    :type state: int
    :param system_pars: The system input parameters defined in the input file
    :type system_pars:
        :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    :return: The topology-output object
    :rtype: :py:class:`~gmak.systems.TopologyOutput`

    .. note:: If your custom system type is to be used with other GROMACS
       compatible objects (e.g., with the :ref:`overview/protocols:general
       protocol[GmxProtocol]`), you should use
       :py:class:`~gmak.gmx_system.GmxTopologyOutput` and set this function to::

           import gmak.gmx_system.GmxTopologyOutput as GmxTopologyOutput
           return GmxTopologyOutput(workdir, name, grid, state)

    .. note:: For even more flexibility, you can create your own customized
       topology-output class, inheriting from :py:class:`~gmak.systems.TopologyOutput`
       (as done e.g. for :py:class:`~gmak.gmx_system.GmxTopologyOutput`).
       In this way, you can provide custom methods as well as data to the topology.

    """
    pass


def topo_out_writer(params, topo_out, system_pars):
    """
    Applies the interaction-parameter values to the topology-output object.  If
    any topology files need to be written, this should also be done in this
    function.

    :param params: The list of interaction parameters to be applied to the
        topology
    :type params: list of
        :py:class:`~gmak.interaction_parameter.InteractionParameter`
    :param topo_out: The :py:class:`~gmak.systems.TopologyOutput` instance
        returned by the :py:func:`~gmak.api_signatures.topo_out_creator` function.
    :type topo_out: :py:class:`~gmak.systems.TopologyOutput`
    :param system_pars: The system input parameters defined in the input file
    :type system_pars:
        :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
    """
    pass
