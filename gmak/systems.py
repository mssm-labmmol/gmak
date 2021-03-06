from gmak.interaction_parameter import InteractionParameter
from gmak.custom_attributes import CustomizableAttributesMixin
from typing import Optional, Callable, List


class TopologyOutput:
    """
    This class encodes a topology. It is intended to be used as a
    generalization of the path of the topology file, when more flexibility is
    needed in constructing it. It is essentially a container for data that the
    user wants to transfer around and that is identified by the program as
    representing the topology of a system. The way to store data in objects
    of this type is via instance attributes, which can be freely set by the
    user. These objects are created individually for each system, grid-shift
    iteration and grid point (:ref:`overview/general_workflow:general
    workflow`).
    """
    def __init__(self):
        pass

    def get_system(self):
        return self._system_ref


class System(CustomizableAttributesMixin):

    def __init__(self,
                 name: str,
                 type: str,
                 topo_out_creator: Callable[
                     [str, str, int, int,
                         CustomizableAttributesMixin.InputParameters],
                     TopologyOutput],
                 topo_out_writer: Callable[
                     [List[InteractionParameter], TopologyOutput,
                      CustomizableAttributesMixin.InputParameters], None]):
        self.name = name
        self.type = type
        self.topo_out_creator = topo_out_creator
        self.topo_out_writer = topo_out_writer

    def write_topology(self,
                       workdir: str,
                       grid: int,
                       state: int,
                       params: List[InteractionParameter]):
        topo_out = self.topo_out_creator(workdir,
                                         self.name,
                                         grid,
                                         state,
                                         self.get_custom_attributes())
        self.topo_out_writer(params, topo_out, self.get_custom_attributes())
        topo_out._system_ref = self
        return topo_out


class CustomSystemFactory:

    ptable = {}

    @classmethod
    def add_custom_system(cls,
                          type_name: str,
                          topo_out_creator: Callable[
                              [str, str, int, int,
                                  CustomizableAttributesMixin.InputParameters],
                              TopologyOutput],
                          topo_out_writer: Callable[
                              [List[InteractionParameter], TopologyOutput,
                               CustomizableAttributesMixin.InputParameters],None]):
        cls.ptable[type_name] = topo_out_creator, topo_out_writer

    @classmethod
    def create(cls,
               type_name: str,
               name: str):
        try:
            return System(name, type_name, cls.ptable[type_name][0],
                          cls.ptable[type_name][1])
        except KeyError:
            raise ValueError(f"Unknown custom-system type \"{type_name}\".")


def create_system(type_name: str,
                  name: str,
                  blockdict: dict,
                  systems: dict,
                  coords: dict,
                  prots: dict,
                  exclude: Optional[List[str]] = None):
    from gmak.gmx_system import GmxSystem
    if type_name == "gmx":
        system = GmxSystem(name)
    else:
        system = CustomSystemFactory.create(type_name, name)
    system.set_custom_attributes_from_blockdict(
        blockdict, systems, coords, prots, exclude=exclude)
    return system
