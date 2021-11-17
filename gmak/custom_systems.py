def gmx_custom_parameter_writer(param, istream, ostream):
    """
    Example of the signature of a function used to apply the values of custom
    parameters to a GROMACS topology file.

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


def topo_out_creator(workdir, name, grid, state, system_attrs):
    """
    Returns a :py:class:`~gmak.systems.TopologyOutput` instance that stores as
    attributes the data necessary to write the topology file. It is intended to
    be used as a generalization of the path of the topology file when more
    flexibility is needed.

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
    :param system_attrs: The system attributes defined in the input file
    :type system_attrs:
        :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.CustomizableAttributesData`
    :return: A :py:class:`~gmak.systems.TopologyOutput` instance
        containing any data the user considers necessary to write the topology
        files later on in the workflow. This instance can be initialized with a
        simple ``TopologyOutput()`` statement, and, after that, the user can
        set the instance attributes as desired. The returned instance is the
        same one passed as a parameter in
        :py:func:`~gmak.custom_systems.topo_out_writer`. It is also passed to
        the :py:func:`~gmak.custom_protocols.simulator` function.
    :rtype: :py:class:`~gmak.systems.TopologyOutput`

    .. note:: If your custom system type is to be used with other GROMACS
       compatible objects (e.g., with the :ref:`overview/protocols:general
       protocol[GmxProtocol]`), you should use
       :py:class:`~gmak.gmx_system.GmxTopologyOutput` and set this function
       to::

           import gmak.gmx_system.GmxTopologyOutput as GmxTopologyOutput
           return GmxTopologyOutput(workdir, name, grid, state)
    """
    pass


def topo_out_writer(params, topo_out, system_attrs):
    """
    Writes the topology files with the interaction-parameter values properly
    applied to it. The topology file paths should be somehow encoded in the
    input topology-output object.

    :param params: The list of interaction parameters to be applied to the
        topology.
    :type params: list of
        :py:class:`~gmak.interaction_parameter.InteractionParameter`
    :param topo_out: The :py:class:`~gmak.systems.TopologyOutput` instance
        returned by the :py:func:`~gmak.custom_systems.topo_out_creator` function.
    :type topo_out: :py:class:`~gmak.systems.TopologyOutput`
    :param system_attrs: The system attributes defined in the input file
    :type system_attrs:
        :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.CustomizableAttributesData`
    """
    pass
