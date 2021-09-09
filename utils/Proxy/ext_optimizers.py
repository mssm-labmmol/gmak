from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable, List, Tuple
from scipy.optimize import (OptimizeResult,
                            minimize,
                            differential_evolution,
                            basinhopping)
import numpy as np
import objective_function

class BaseOptimizer(ABC):

    @abstractmethod
    def run(self, *args, **kwargs):
        pass


@dataclass
class SimplexOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    x0: np.ndarray

    def run(self, *args, **kwargs) -> OptimizeResult:
        return minimize(self.obj_func.eval, self.x0, method='Nelder-Mead')


@dataclass
class CGOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    x0: np.ndarray
    #jac: Callable[[np.ndarray], float]

    def run(self, *args, **kwargs) -> OptimizeResult:
        return minimize(self.obj_func.eval, self.x0, method='CG')


@dataclass
class DEOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    popsize: int
    bounds: Tuple[Tuple[float, float],...]

    def run(self, *args, **kwargs) -> OptimizeResult:
        return differential_evolution(self.obj_func.eval, self.bounds)


@dataclass
class BasinHopOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    x0: np.ndarray

    def run(self, *args, **kwargs) -> OptimizeResult:
        return basinhopping(self.obj_func.eval, self.x0)


# -----------------------------------
# yMHG  Sat Aug 21 21:01:09 -03 2021 
# These are more complicated and will be implemented later on.
# -----------------------------------

@dataclass
class GPRUnrestrictedOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    x0s: List[np.ndarray]


@dataclass
class CMAESOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    mu_0: float
    s_0: float
    popsize: int


