from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable, List, Tuple
from scipy.optimize import minimize, OptimizeResult
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
class GPRUnrestrictedOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    x0s: List[np.ndarray]


@dataclass
class CMAESOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    mu_0: float
    s_0: float
    popsize: int


@dataclass
class DEOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    bounds: Tuple[Tuple[float, float],...]
    popsize: int
    

@dataclass
class BasinHopOptimizer(BaseOptimizer):

    obj_func: objective_function.BaseObjectiveFunction
    x0: np.ndarray

