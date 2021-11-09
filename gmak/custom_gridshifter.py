from gmak.gridshifter import add_custom_gridshifter

def calculator(tuple_indexes,
               scores,
               propnames,
               averages,
               uncertainties):
    """
    The grid-shifting procedure.

    :param tuple_indexes: The tuple indexes of all grid points, ordered by
        linear index.
    :type tuple_indexes: list
    :param scores: A list with the scores of the grid points, ordered by linear
        index.
    :type scores: list
    :param propnames: The names of the calculated properties.
    :type propnames: list
    :param averages: A 2D array with the estimated values of the properties for
        all grid points.  The first index is the linear index of the grid
        point.  The second index corresponds to the property index in
        ``propnames``.
    :type averages: np.ndarray
    :param uncertainties: A 2D array with the estimated uncertainties of the
        properties for all grid points.  The first index is the linear index of
        the grid point.  The second index corresponds to the property index in
        ``propnames``.
    :type uncertainties: np.ndarray
    :return: The tuple index of the new origin.
    :rtype: tuple
    """
    pass


