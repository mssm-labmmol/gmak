from abc import ABC, abstractmethod

class BaseObjectiveFunction(ABC):
    @abstractmethod
    def eval(self, xs):
        pass

    @abstractmethod
    def get_estimators(self):
        pass


class WRMSD_ObjectiveFunction(BaseObjectiveFunction):
    def __init__(self, rs, ws, ss):
        """
        Parameters:

        rs: list of reference values for each property
        ws: list of weights for each property
        ss: list of estimators for each property
        """
        self.rs = rs
        self.ws = ws
        self.ss = ss
        if len(rs) != len(ws):
            raise ValueError

    def get_estimators(self):
        return self.ss

    def eval(self, xs):
        """Parameters:

        xs: x value (tuple or list)
        """
        sum = 0
        sum_w = 0
        for s, r, w in zip(self.ss, self.rs, self.ws):
            x = s.get_estimate(xs)
            sum += w *  (x - r) ** 2
            sum_w += w
        sum /= sum_w
        sum = (sum) ** (0.5)
        return sum


class RealObjectiveFunction(BaseObjectiveFunction):
    def __init__(self, estimator):
        """
        Parameters: estimator of the real function.
        """
        self.ss = estimator

    def eval(self, xs):
        return self.ss.get_estimate(xs)

    def get_estimators(self):
        return [self.ss,]
