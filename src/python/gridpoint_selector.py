
from abc import ABC, abstractmethod
import random

class GridpointSelectorInterface(ABC):

    @abstractmethod
    def selectPoints(self): pass

class RandomSelector(GridpointSelectorInterface):

    def __init__(self, linearSize, howMany):
        self.linearSize = linearSize
        self.howMany = howMany

    def selectPoints(self):
        return [random.randint(0, self.linearSize) for i in range(self.howMany)]
        
class BestSelector(GridpointSelectorInterface):

    def __init__(self, gridOptimizer, howMany):
        self.gridOptimizer = gridOptimizer
        self.howMany = howMany

    def selectPoints(self):
        return self.gridOptimizer.getRankedBest(self.howMany)

def createSelector(typeString, *args, **kwargs):
    options_dict = {
        'random': RandomSelector,
        'best': BestSelector,
    }
    return options_dict[typeString](*args, **kwargs)
