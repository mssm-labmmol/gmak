from protocols import add_custom_protocol

# Example
# -------
# 
# Simulator is a function that receives the current topology, a dictionary
# with the args specified in the input file (first field is key; remaining 
# fields are in a list which are the values) and the workdir of the simulation.
#
# It returns a dictionary with the output files, where keys are extensions
# such as xtc, edr, gro, etc. You may define custom keys to interact with
# custom properties also.
def example_simulator(top, args_dict, workdir):
    import os
    os.system(f"echo 1 20 3 4 {top}> {workdir}/Test.dat")
    return {'dat': f"{workdir}/Test.dat"}

add_custom_protocol("example-protocol", example_simulator)

# TODO deal with custom protocols that do not require mdps... then how will we
# treat extensions of simulations?
#
# TODO better document the addition of custom protocols and custom properties.
