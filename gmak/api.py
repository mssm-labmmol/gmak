from gmak.component_properties import CustomAtomicProperty
from gmak.gridoptimizer import ScoreFactory
from gmak.gridshifter import CustomGridShifterFactory
from gmak.property import CustomPropertyFactory
from gmak.protocols import CustomProtocolFactory
from gmak.surrogate_model import CustomSurrogateModelFactory
from gmak.systems import CustomSystemFactory
from gmak.gmx_system import GmxCustomReplace


def add_custom_component_property(type_name, component_calculator, is_timeseries):
    """
    Adds a custom component-property type to the program. In the input file, it
    can be referenced with the type ``type_name``.

    :param type_name: Name of the type of the custom component property.
    :type type_name: str
    :param component_calculator: The function used to calculate the custom
        component property (see
        :py:func:`~gmak.api_signatures.component_calculator`)
    :type component_calculator: callable
    :param is_timeseries: ``True`` indicates that the component property is
        obtained as a timeseries; ``False``, as a tuple ``(EA, dEA)`` with the
        expected value and statistical uncertainty.
    :type is_timeseries: bool
    """
    CustomAtomicPropertyFactory.add_custom_component_property(type_name,
                                                           component_calculator,
                                                           is_timeseries)


def add_gmx_custom_parameter_writer(funct):
    """
    Adds to the program a function to apply the values of custom parameters
    (see :ref:`overview/systems_and_topologies:custom parameters`) to a GROMACS
    topology file.

    :param funct: The function used to apply the values of custom parameters
        (see :py:meth:`~gmak.api_signatures.gmx_custom_parameter_writer`)
    :type funct: callable
    """
    GmxCustomReplace.add_gmx_custom_parameter_writer(funct)


def add_custom_score(type_name, calc_score, calc_score_err=None):
    """
    Adds a custom score-function type to the program. In the input file, it can
    be referenced with the type ``type_name``.

    :param type_name: Name of the type of the custom function
    :type type_name: str
    :param calc_score: The score function (see :py:func:`~gmak.api_signatures.calc_score`)
    :type calc_score: callable
    :param calc_score_err: (optional) The score uncertainty function (see
        :py:func:`~gmak.api_signatures.calc_score_err`). Defaults to
        :py:obj:`None`, which means that the uncertainties are not calculated.
    :type calc_score_err: callable
    """
    ScoreFactory.add_custom_score(type_name, calc_score, calc_score_err)


def add_custom_gridshifter(type_name, calculator):
    """
    Adds a custom grid-shifting procedure to the program. In the input file, it
    can be referenced with the type ``type_name``.

    :param type_name: Name of the type of the custom protocol.
    :type type_name: str
    :param calculator: The custom grid-shifting function (see
        :py:func:`~gmak.api_signatures.calculator`).
    :type calculator: callable
    """
    CustomGridShifterFactory.add_custom_gridshifter(type_name, calculator)


def add_custom_composite_property(type_name, composite_calculator=None):
    """
    Adds a custom composite-property type to the program. In the input file, it
    can be referenced with the type ``type_name``.

    :param type_name: Name of the type of the custom composite property.
    :type type_name: str
    :param composite_calculator: (optional) The calculator function (see
        :py:func:`~gmak.api_signatures.composite_calculator`). If it is
        not supplied, the program implicitly assumes that the property has only
        one component and identifies the values and errors of the composite
        property with those of the component property.
    :type composite_calculator: callable
    """
    return CustomPropertyFactory.add_custom_composite_property(type_name,
                                                               composite_calculator)


def add_custom_protocol(type_name, simulator, calc_initial_len=None,
                        calc_extend=None, get_last_frame=None):
    """
    Adds a custom protocol type to the program. In the input file, it can be
    referenced with the type ``type_name``.

    :param type_name: Name of the type of the custom protocol.
    :type type_name: str
    :param simulator: The simulator function (see
        :py:func:`~gmak.api_signatures.simulator`)
    :type simulator: callable
    :param calc_initial_len: (optional) The function to calculate the initial
        length of the production run (see
        :py:func:`~gmak.api_signatures.calc_initial_len`). This is used only
        when support for extensions is desired. Defaults to :py:obj:`None`,
        indicating that extensions are not desired.
    :type calc_initial_len: callable
    :param calc_extend: (optional) The function to calculate the new length of
        the production run (see :py:func:`~gmak.api_signatures.calc_extend`).
        This is used only when support for extensions is desired. Defaults to
        :py:obj:`None`, indicating that extensions are not desired.
    :type calc_extend: callable
    :param get_last_frame: (optional) The function to extract the configuration
        file to be used as a starting point for the following protocol (see
        :py:func:`~gmak.api_signatures.get_last_frame`).  This is used only when a
        followable protocol is desired.  Defaults to :py:obj:`None`.
    :type get_last_frame: callable
    """
    CustomProtocolFactory.add_custom_protocol(type_name, simulator,
                                              calc_initial_len, calc_extend,
                                              get_last_frame)


def add_custom_surrogate_model(type_name, compute, corners=False):
    """
    Adds a custom surrogate-model type to the program. In the input file, it
    can be referenced with the type ``type_name``.

    :param type_name: Name of the type of the custom surrogate model
    :type type_name: str
    :param compute: The surrogate-model function (see
        :py:func:`~gmak.api_signatures.compute`)
    :type compute: callable
    :param corners: Indicates whether the surrogate model requires the corners
        of the grid to be simulated (e.g., interpolation does). Defaults to ``False``.
    :type corners: bool
    """
    CustomSurrogateModelFactory.add_custom_surrogate_model(type_name,
                                                           compute,
                                                           corners)


def add_custom_system(type_name, topo_out_creator, topo_out_writer):
    """
    Adds a custom system type to the program. In the input file, it
    can be referenced with the type ``type_name``.

    :param type_name: Name of the type of the custom system.
    :type type_name: str
    :param topo_out_creator: A function that creates a
        :py:class:`~gmak.systems.TopologyOutput` object for a given
        system, grid-shift iteration and grid point (see
        :ref:`overview/general_workflow:general workflow`
        and :py:func:`~gmak.api_signatures.topo_out_creator`).
    :type topo_out_creator: callable
    :param topo_out_writer: A function that applies the interaction-parameter
        values to the topology-output object (see
        :py:func:`~gmak.api_signatures.topo_out_writer`)
    :type topo_out_writer: callable
    """
    CustomSystemFactory.add_custom_system(type_name, topo_out_creator,
                                          topo_out_writer)


