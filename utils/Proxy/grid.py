from copy import deepcopy
import numpy as np
import itertools

class Grid:

    def __init__(self, origin, spacing, shape):
        self.origin = deepcopy(origin)
        self.spacing = deepcopy(spacing)
        self.shape = deepcopy(shape)
        self.update()

    @staticmethod
    def from_limits_spacing(bottom_limit,
                            upper_limit,
                            spacing):

        bottom_limit_arr = np.array(bottom_limit)
        upper_limit_arr  = np.array(upper_limit)
        spacing_arr      = np.array(spacing)

        shape_arr = (upper_limit_arr - bottom_limit_arr)/spacing_arr
        shape     = [int(x) for x in shape_arr]

        if (shape_arr - shape > 0.05).any():
            raise ValueError("Recheck your limits or your spacing.")

        return Grid(bottom_limit_arr, spacing_arr, tuple(shape))

    def extract_pattern(self, snap_to_grid=True):
        # only applies to 2D arrays...
        if self.get_dim() != 2:
            raise ValueError
        iter_arg = [(0, n-1) for n in self.shape]
        iter_obj = itertools.product(*iter_arg)
        points   = [p for p in iter_obj]
        corners  = [self.origin + p * self.spacing for p in points]
        # center is halfway between first and last corners
        center   = .5 * (corners[0] + corners[-1])
        # now add points between center and corners
        first_layer  = [.5 * (center + c) for c in corners]
        second_layer = [
            .5 * (first_layer[0] + first_layer[1]),
            .5 * (first_layer[2] + first_layer[3]),
            .5 * (first_layer[0] + first_layer[2]),
            .5 * (first_layer[1] + first_layer[3]),
        ]
        unsnapped = corners + [center,] + first_layer + second_layer
        snapped   = [self.snap_to_grid(x) for x in unsnapped]
        return snapped
        
    def update(self):
        self.idxs, self.gridpoints = self._calc_ij_xs_ndarray()

    def shift(self, transl):
        import numpy as np
        self.origin += np.array(transl)
        self.update()

    def recenter(self, center):
        curr_center = self.get_center()
        shift_vector = (center - curr_center)
        self.shift(shift_vector)

    def get_center(self):
        center_loc = [int((s-1)/2) for s in self.shape]
        for xi, ii in zip(self.gridpoints, self.idxs):
            if (ii == center_loc).all():
                return xi

    def get_dim(self):
        return len(self.shape)

    def get_size(self):
        return np.prod(self.shape)

    def locate(self, xs):
        for xi, ii in zip(self.gridpoints, self.idxs):
            if ((xi == xs).all()):
                return ii
        raise IndexError

    def locate_linear(self, xs):
        for i, xi in enumerate(self.gridpoints):
            if ((xi == xs).all()):
                return i
        raise IndexError

    def locate_frac(self, xs):
        """Locate point and return its position in fractional coordinates."""
        loc = self.locate(xs)
        last_point = self.gridpoints[-1]
        loc_frac = tuple([(xs[i] - self.origin[i])/(last_point[i] - self.origin[i])
                          for i in range(self.get_dim())])
        return loc_frac
                          
    def _calc_ij_xs_ndarray(self):
        import numpy as np
        from itertools import product
        def generator_wrapper(idxs):
            return np.array(self.origin) + \
                np.array(self.spacing) * idxs
        ranges = map(np.arange, self.shape)
        idxs = list(product(*ranges))
        gridpoints = list(map(generator_wrapper, idxs))
        return np.array(idxs), np.array(gridpoints)

    def snap_to_grid(self, xs):
        xsarr = np.array(xs)
        loc_as_float = (xsarr - self.origin)/(self.spacing)
        loc_as_int = list(map(int, loc_as_float))
        out = (self.origin + np.array(loc_as_int) * self.spacing)
        return out

    def to_ij_ndarray(self):
        return self.idxs

    def to_x_ndarray(self):
        return self.gridpoints

    def to_griddata(self):
        import numpy as np
        return tuple(np.meshgrid(
            *[[self.origin[d] + self.spacing[d]*j
               for j in range(i)]
              for d,i in enumerate(self.shape)]))
