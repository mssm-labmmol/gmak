from abc import ABC, abstractmethod
import gmak.traj_ana as traj_ana
from gmak.custom_attributes import CustomizableAttributesMixin
import numpy as np

class PropertyNotInitialized(Exception):
    pass


class GmxCalculatorMixin:
    def calculator(self, input_traj: dict, out_path: str):
        # member dependencies:
        #
        # self.gmx_name
        traj_ana.analyzeWrapper(input_traj, self.gmx_name, out_path)


class BaseAtomicProperty(ABC,
                         CustomizableAttributesMixin):
    # basic attributes:
    #
    # name: str
    #
    # is_timeseries: bool
    #
    # calculator: callable(
    #   input_traj: dict({'xtc': ..., 'edr': ..., ...}),
    #   out_path: str)
    def __eq__(self, other):
        if isinstance(other, str):
            return self.name == other
        elif isinstance(other, BaseAtomicProperty):
            return self.name == other.name
        else:
            raise NotImplementedError

    def calc(self, input_traj, out_path, state=None):
        return self.calculator(input_traj, out_path)


class GmxEnergyProperty(BaseAtomicProperty, GmxCalculatorMixin):
    def __init__(self, name, gmx_name):
        """
        Parameters:

        name (str): name of the property as identified in the gmak program.

        gmx_name (str): name of the property as listed in the `gmx energy`
                        command.
        """
        self.name = name
        self.gmx_name = gmx_name
        self.is_timeseries = True


# GmxPolcorr is an exception---it is not really a GmxEnergyProperty.
class GmxPolcorr(BaseAtomicProperty):

    def __init__(self, mu, alpha):
        self.name = "polcorr"
        self.is_timeseries = True
        if (mu is None) or (alpha is None):
            raise ValueError(f"Polarization correction requires float values of "
                             f"\"mu\" and \"alpha\".")
        self.mu = mu
        self.alpha = alpha

    def calc(self, input_traj, out_path, state=None):
        traj_ana.obtain_polcorr(input_traj['xtc'],
                                input_traj['edr'],
                                input_traj['gro'],
                                input_traj['tpr'],
                                self.mu,
                                self.alpha,
                                out_path)


# exception---it is not really a GmxEnergyProperty.
class GmxDG(BaseAtomicProperty):

    def __init__(self, temperature):
        self.name = "dg"
        self.temperature = temperature
        self.is_timeseries = False

    def calc(self, input_traj, out_path, state=None):
        # in this particular case, we expect each member of input_traj to be
        # a list of files for each core protocol
        # ------------------------------------------------------------------------
        #  Alchemlyb imports
        # ------------------------------------------------------------------------
        from   alchemlyb import preprocessing
        from   alchemlyb.parsing import gmx
        from   alchemlyb import estimators
        import pandas as pd
        import pymbar

        dhdlFiles = input_traj['dhdl']

        # calculate kbT value
        kbT = self.temperature * 0.83144626
        # intialize MBAR object with default settings
        mbar = estimators.MBAR()
        # extract u_nk from the dHdl files
        u_nk = [gmx.extract_u_nk(dhdl_file, self.temperature)
                for dhdl_file in dhdlFiles]
        # subsample within each u_nk
        u_nk = [preprocessing.statistical_inefficiency(df, df.iloc[:,i])
                for i, df in enumerate(u_nk)]
        # concatenate u_nk for all sampled states
        u_nk = pd.concat(u_nk)
        # run MBAR on the u_nk data
        mbar.fit(u_nk)
        # get free-energy-difference estimate and estimate error
        df  = kbT * mbar.delta_f_.iloc[0,-1]
        ddf = kbT * mbar.d_delta_f_.iloc[0,-1]
        # write these values to disk
        _fp = open(out_path, 'w')
        _fp.write('{:>10.3f}\n{:>10.4f}\n'.format(df, ddf))
        _fp.close()


class GmxPropertyFactory:

    @classmethod
    def create_gmx_component_property(cls, name, *args, **kwargs):
        if (name == 'density'):
            return GmxEnergyProperty("density", "Density")
        elif (name == 'potential'):
            return GmxEnergyProperty("potential", "Potential")
        elif (name == 'volume'):
            return GmxEnergyProperty("volume", "Volume")
        elif (name == 'gamma'):
            return GmxEnergyProperty("gamma", "#Surf*SurfTen")
        elif (name == 'polcorr'):
            return GmxPolcorr(*args, **kwargs)
        elif (name == 'dg'):
            return GmxDG(*args, **kwargs)
        elif (name.startswith("gmx_")):
            return GmxEnergyProperty(name, name[4:])
        else:
            raise PropertyNotInitialized


class CustomAtomicProperty(BaseAtomicProperty):

    def __init__(self, name, calculator, is_timeseries, *args, **kwargs):
        self.name = name
        self._calculator = calculator
        self.is_timeseries = is_timeseries
        self.system = kwargs['system']

    def calc(self, input_traj, out_path, state=None):
        return self.calculator(input_traj, out_path, state)

    def calculator(self, input_traj, out_path, state):
        topology = self.system.getPathsForStatepath(state)
        prop_data = self._calculator(topology,
                                     input_traj,
                                     self.get_custom_attributes())
        np.savetxt(out_path, prop_data)


class CustomAtomicPropertyFactory:

    # ptable is a dict { name: {calculator, is_timeseries} }
    ptable = {}

    @classmethod
    def add_custom_component_property(cls, name, calculator, is_timeseries):
        cls.ptable[name] = {'calculator': calculator,
                            'is_timeseries': is_timeseries}

    @classmethod
    def create_custom_component_property(cls, name, *args, **kwargs):
        return CustomAtomicProperty(name, cls.ptable[name]['calculator'],
                                    cls.ptable[name]['is_timeseries'],
                                    *args, **kwargs)

def add_custom_component_property(type_name, component_calculator, is_timeseries):
    """
    Adds a custom component property to the program. In the input file, it can
    be referenced with the type ``type_name``.

    :param type_name: Name of the type of the custom component property.
    :type type_name: str
    :param component_calculator: The function used to calculate the custom
        component property (see
        :py:func:`~gmak.custom_properties.component_calculator`)
    :type component_calculator: callable
    :param is_timeseries: ``True`` indicates that the component property is
        obtained as a timeseries; ``False``, as a tuple ``(EA, dEA)`` with the
        expected value and statistical uncertainty.
    :type is_timeseries: bool
    """
    CustomAtomicPropertyFactory.add_custom_component_property(type_name,
                                                           component_calculator,
                                                           is_timeseries)


def create_component_property(name, *args, **kwargs):
    ptable = [
        GmxPropertyFactory.create_gmx_component_property,
        CustomAtomicPropertyFactory.create_custom_component_property,
    ]
    while True:
        factory = ptable.pop(0)
        try:
            return factory(name, *args, **kwargs)
        except PropertyNotInitialized:
            continue
    raise PropertyNotInitialized(f"Couldn't create atomic property `{name}`.")
