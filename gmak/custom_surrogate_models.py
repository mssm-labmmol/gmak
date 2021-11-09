from gmak.surrogate_model import add_custom_surrogate_model
from gmak.cartesiangrid import flat2tuple

def compute(EA_s, dEA_s, I_s, gridshape, X_ki):
    """
    Computes the expected values and uncertainties of a given property for the
    entire grid.

    :param EA_s: A 1D array of shape ``(NSAMP,)`` with the average values of the
        property for each sampled point.
    :type EA_s: np.ndarray
    :param dEA_s: A 1D array of shape ``(NSAMP,)`` with the uncertainties of the
        property for each sampled point.
    :type dEA_s: np.ndarray
    :param I_s: A list of length ``NSAMP``with the linear indexes of the
        sampled grid points.
    :type I_s: list
    :param gridshape: A tuple with the grid dimensions.
    :type gridshape: tuple
    :param X_ki: A 2D array with the parameter-space values for each grid point
        (more specifically, of the :ref:`main variation`). The first index is
        the linear index, and the second index is the coordinate the
        parameter-space point.  For example, ``X_ki[1,0]`` is the value of
        the first parameter of the :ref:`main variation` for the grid point
        with linear index equal to 1.
    :type X_ki: np.ndarray
    :return: A tuple ``(EA_k, dEA_k)`` containing the estimated values and
        uncertainties, respectively, of the property for each grid point.
        ``EA_k`` and ``dEA_k`` are 1D arrays indexed by the linear index of
        the grid point.
    :rtype: tuple


    .. note:: The function :py:func:`gmak.cartesiangrid.flat2tuple` can be
       used to convert a linear to a tuple index if desired.
    """
    pass

