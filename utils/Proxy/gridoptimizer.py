from abc import ABC, abstractmethod
import numpy as np
import grid
import estimator
import objective_function


# TODO: Init estimators from inside function, and use samples as
# parameters
def gmak_from_parameters(grid_origin, grid_spacing, grid_shape,
                         estimators, 
                         obj_refs, obj_weis,
                         shift_cut, 
                         conv_margins,
                         opt_type='min',
                         fix_samples=True, keep_samples=False):

    initial_grid = grid.Grid(origin=grid_origin,
                        spacing=grid_spacing,
                        shape=grid_shape)
    obj_function = objective_function.WRMSD_ObjectiveFunction(
        rs=obj_refs,
        ws=obj_weis,
        ss=estimators)
    return StandardGridMaker(initial_grid,
                             obj_function,
                             shift_cut,
                             conv_margins,
                             opt_type=opt_type, fix_samples=fix_samples,
                             keep_samples=keep_samples)


# ========================= Shifter Hierarchy ========================

# Shift-exclusive operations.

class BaseGridShifterMixin(ABC):
    @abstractmethod
    def shift_grid(self, new_center):
        pass
    


class StandardGridShifterMixin(BaseGridShifterMixin):
    """Shifts the grid possibly fixing or keeping samples."""
    def init_gridshifter(self, fix_samples=True, keep_samples=False):
        self.fix = fix_samples
        self.keep = keep_samples

    def shift_grid(self, new_center):
        if not self.fix:
            raise NotImplementedError
        else:
            # store current sample indexes
            sample_idxs = [self.grid.locate_linear(s)
                           for s in self.get_samples()]
            # shift grid
            self.grid.recenter(new_center)
            # set samples 
            new_gridpoints = self.grid.to_x_ndarray()
            new_samples = [new_gridpoints[s] for s in sample_idxs]
            self.update_samples(new_samples)
        if self.keep:
            raise NotImplementedError
        
    
# Determination of new center.

class BaseCenterUpdateMixin(ABC):
    @abstractmethod
    def calc_new_center(self):
        pass

class StandardCenterUpdateMixin(BaseCenterUpdateMixin):
    def init_center_update(self, shift_cut, opt_type='min', snap_to_grid=True):
        self.shift_cut = shift_cut
        self.opt_type = opt_type
        self.snap_to_grid = snap_to_grid

    def calc_new_center(self):
        """New center is CG of shift_cut * number_of_points best
        points."""
        func_vals = [self.calc_obj_func(xs)
                     for xs in self.grid.to_x_ndarray()]
        sorted_idxs = np.argsort(func_vals)
        if self.opt_type == 'max':
            sorted_idxs = np.flip(sorted_idxs)
        actual_points = [self.grid.to_x_ndarray()[s] for s in sorted_idxs]
        npoints = int(self.shift_cut * self.grid.get_size())
        cg = np.mean(actual_points[:npoints], axis=0)
        if self.snap_to_grid:
            cg = self.grid.snap_to_grid(cg)
        return cg
        
# ===================== End of Shifter Hierarchy =====================

class BaseGridMaker(ABC):
    def calc_obj_func(self, xs):
        return self.objective_function.eval(xs)
        
    def get_samples(self):
        return self.objective_function.get_estimators()[0].samples

    def update_samples(self, new_samples):
        for es in self.objective_function.get_estimators():
            es.replace_samples(new_samples) # this should
                                              # automatically update
                                              # the objective function
    def get_min(self):
        func_vals = [self.calc_obj_func(xs)
                    for xs in self.grid.to_x_ndarray()]

        sorted_idxs = np.argsort(func_vals)
        actual_points = [self.grid.to_x_ndarray()[s] for s in sorted_idxs]
        func_vals = np.sort(func_vals)
        return actual_points[0], func_vals[0]


    @abstractmethod
    def has_converged(self):
        pass


    # This assumes the Mixins exist, but does not necessarily specify
    # them.
    def run(self):
        self.current_grid = 0
        print("---- Grid %d ----" % self.current_grid)
        np.savetxt("grid_%d.dat" % self.current_grid,
                    [self.calc_obj_func(xs)
                    for xs in self.grid.to_x_ndarray()])
        while not self.has_converged():
            np.savetxt("grid_%d.dat" % self.current_grid,
                       [self.calc_obj_func(xs)
                        for xs in self.grid.to_x_ndarray()])
            new_center = self.calc_new_center()
            self.shift_grid(new_center)
            self.current_grid += 1
            print("---- Grid %d ----" % self.current_grid)
            print("Origin = ", self.grid.origin)


class StandardGridMaker(BaseGridMaker,
                        StandardGridShifterMixin,
                        StandardCenterUpdateMixin):
    
    def __init__(self,
                 grid,
                 objective_function,
                 shift_cut,
                 conv_margins, # list of (min,max) tuples for each dimension
                 snap_to_grid=True,
                 opt_type='min',
                 fix_samples=True,
                 keep_samples=False,
                 maxshifts=10):

        # Init Mixins
        self.init_gridshifter(fix_samples=fix_samples,
                              keep_samples=keep_samples)
        self.init_center_update(shift_cut=shift_cut,
                                opt_type=opt_type,
                                snap_to_grid=snap_to_grid)
        self.objective_function = objective_function
        self.grid = grid
        self.maxshifts = maxshifts
        self.conv_margins = conv_margins

    def has_converged(self):
        near_center = True
        if self.snap_to_grid:
            new_center = self.grid.locate_frac(self.grid.snap_to_grid(self.calc_new_center()))
        else:
            new_center = self.grid.locate_frac(self.calc_new_center())
        print("Opt. CG located at (frac):", new_center)
        for i in range(self.grid.get_dim()):
            if (new_center[i] < self.conv_margins[i][0]) or (new_center[i] > self.conv_margins[i][1]):
                near_center = False
                break
        if (self.current_grid <= self.maxshifts) and not (near_center):
            print("Has not converged")
            return False
        else:
            print("Converged")
            return True
