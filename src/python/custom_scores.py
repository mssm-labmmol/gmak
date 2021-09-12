from gridoptimizer import add_custom_score

# example score -- always returns 0.0
def score_example_1(estimates, errors, weis, refs):
    return 0.0

add_custom_score("my-score-1",
                 calc_score=score_example_1)

# example score -- always returns 0.0 with 'error interval' (-1.0, +0.5)
score_example_2 = score_example_1

def score_err_example_2(estimates, errors, weis, refs):
    return (-1.0, +0.5)

add_custom_score("my-score-2",
                 calc_score=score_example_2,
                 calc_score_err=score_err_example_2)


