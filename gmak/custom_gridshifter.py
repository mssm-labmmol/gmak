from gmak.gridshifter import add_custom_gridshifter

# test - a shifter that always shifts to the first point
def shifter_always_last(tuple_indexes,
                        scores,
                        propnames,
                        averages,
                        uncertainties):
    return tuple_indexes[-1]


add_custom_gridshifter("first", shifter_always_last)
