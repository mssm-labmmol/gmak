from surrogate_model import add_custom_surrogate_model
from cartesiangrid import flat2tuple

"""
    To add your own surrogate model to the program, there are two
steps.  First, you must implement the core function that calculates
the property estimates and errors for every point of the
parameter-space grid based on the estimates and errors of the sampled
parameter-space points. This function (let's say, `compute`) must have
the following signature:

    compute(EA_s, dEA_s, I_s, gridshape, X_ki) -> EA_k, dEA_k

where

Input:
    EA_s: flat 1D np.ndarray with the average values of the properties
          for each sampled point.

    dEA_s: flat 1D np.ndarray with the uncertainties of the properties
           for each sampled point.

    I_s: list with the flat indexes of the sampled points.

    gridshape: tuple with the shape of the grid

    X_ki: np.ndarray with parameter-space values for each grid point.
          The first index refers to the flat position in the grid, and
          the second index refers to the coordinate of the
          parameter-space point.  For example, X_ki[0,0] is the first
          parameter-space coordinate of the first point of the grid.

Returns:
    EA_k: flat 1D np.ndarray with the estimated values of the
          properties for all points in the grid.

    dEA_k: flat 1D np.ndarray with the estimate uncertainties of the
           properties for all points in the grid.

If convenient, you can convert flat indexes into tuple indexes using
the function `flat2tuple` imported from the `cartesiangrid` module.
For example:

>>> flat2tuple((10,10), 5)
(0,5)

>>> flat2tuple((10,10), 10)
(1,0)

In this case, it may also be convenient to convert the flat arrays
`EA_s` and `dEA_s` into multidimensional arrays using the
`numpy.reshape` function, with the parameter `gridshape` as the
desired shape. Just remember that you must return flat arrays in the
end, so, if you use this approach, also use the
`numpy.ndarray.flatten` method to flatten the estimated arrays EA_k,
dEA_k.

"""

def example_compute(EA_s, dEA_s, I_s, gridshape, X_ki):
    import numpy as np
    flattened_size = np.prod(gridshape)
    return np.zeros((flattened_size,)), np.zeros((flattened_size,))


add_custom_surrogate_model("smexample",
                           compute=example_compute,
                           corners=True)
