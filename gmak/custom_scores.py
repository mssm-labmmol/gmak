from gmak.gridoptimizer import add_custom_score

def calc_score(estimates, errors, weis, refs):
    """
    Calculates the score for a grid point.

    :param estimates: The dictionary of property estimates. The keys are the
        property names. Each key maps to the expected value of the property at
        a given grid point.
    :type estimates: dict
    :param errors: The dictionary of property uncertainties. The keys are the
        property names. Each key maps to the value of the uncertainty of the
        property at a given grid point.
    :type errors: dict
    :param weis: The dictionary of property weights. The keys are the
            property names. Each key maps to the weight of the property
            in the optimization.
    :type weis: dict
    :param refs: The dictionary of reference values of the properties. The keys
        are the property names. Each key maps to the reference values of the
        property in the optimization.
    :type refs: dict
    :return: The value of the score
    :rtype: float
    """
    pass


def calc_score_err(estimates, errors, weis, refs):
    """
    Calculates the score uncertainty for a grid point.

    :param estimates: The dictionary of property estimates. The keys are the
        property names. Each key maps to the expected value of the property at
        a given grid point.
    :type estimates: dict
    :param errors: The dictionary of property uncertainties. The keys are the
        property names. Each key maps to the value of the uncertainty of the
        property at a given grid point.
    :type errors: dict
    :param weis: The dictionary of property weights. The keys are the
            property names. Each key maps to the weight of the property
            in the optimization.
    :type weis: dict
    :param refs: The dictionary of reference values of the properties. The keys
        are the property names. Each key maps to the reference values of the
        property in the optimization.
    :type refs: dict
    :return: A tuple ``(MIN, MAX)`` that represents a confidence interval for
        the score
    :rtype: tuple
    """
    pass
