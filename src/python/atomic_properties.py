from abc import ABC, abstractmethod
import traj_ana

class PropertyNotInitialized(Exception):
    pass


class GmxCalculatorMixin:
    def calculator(self, input_traj: dict, out_path: str):
        # member dependencies:
        #
        # self.gmx_name
        traj_ana.analyzeWrapper(input_traj, self.gmx_name, out_path)


class BaseAtomicProperty(ABC):
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

    def calc(self, input_traj, out_path):
        return self.calculator(input_traj, out_path)


class GmxEnergyProperty(BaseAtomicProperty, GmxCalculatorMixin):
    def __init__(self, name, gmx_ana):
        """
        Parameters:

        name (str): name of the property as identified in the gmak program.

        gmx_name (str): name of the property as listed in the `gmx energy`
                        command.
        """
        self.name = name
        self.gmx_name = gmx_name
        self.is_timeseries = True

    @staticmethod
    def create_gmx_atomic_property(name, *args, **kwargs):
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
        else:
            raise PropertyNotInitialized


# GmxPolcorr is an exception---it is not really a GmxEnergyProperty.
class GmxPolcorr(BaseAtomicProperty):
    def __init__(self, mu, alpha):
        self.name = "polcorr"
        self.is_timeseries = True
        self.mu = mu
        self.alpha = alpha

    def calc(self, input_traj, out_path):
        traj_ana.obtain_polcorr(input_traj['xtc'],
                                input_traj['edr'],
                                input_traj['gro'],
                                input_traj['tpr'],
                                self.mu,
                                self.alpha,
                                out_path)


class CustomAtomicProperty(BaseAtomicProperty):
    def __init__(self, name, calculator, is_timeseries):
        self.name = name
        self.calculator = calculator
        self.is_timeseries = is_timeseries


class CustomAtomicPropertyFactory:

    # ptable is a dict { name: {calculator, is_timeseries} }
    ptable = {}

    @classmethod
    def add_custom_atomic_property(cls, name, calculator, is_timeseries):
        cls.ptable[name] = {'calculator': calculator,
                            'is_timeseries': is_timeseries}

    @classmethod
    def create_custom_atomic_property(cls, name):
        return CustomAtomicProperty(name, cls.ptable[name]['calculator'],
                                    cls.ptable[name]['is_timeseries'])


def create_atomic_property(name):
    ptable = [
        GmxEnergyProperty.create_gmx_atomic_property,
        CustomAtomicPropertyFactory.create_custom_atomic_property,
    ]
    while True:
        factory = ptable.pop()
        try:
            return factory(name)
        except PropertyNotInitialized:
            continue
    raise PropertyNotInitialized(f"Couldn't create atomic property `{name}`.")
