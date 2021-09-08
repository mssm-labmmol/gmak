from protocols import add_custom_protocol

# Example
# -------
#
# Simulator is a function that receives the current topology, a
# dictionary with the args specified in the input file (first field is
# key; remaining fields are in a list which are the values) and the
# workdir of the simulation.
#
# It returns a dictionary with the output files, where keys are
# extensions such as xtc, edr, gro, etc. You may define custom keys to
# interact with custom properties also.
def example_simulator(top, args_dict, length, workdir):
    import os
    os.system(f"echo {length} {top}> {workdir}/Test.dat")
    return {'dat': f"{workdir}/Test.dat"}

def example_calc_inital_len(args_dict):
    return 10

def example_calc_extend(errs_tols, length):
    if length > 30:
        print("Reached limit.")
        return None
    else:
        return length + 5

add_custom_protocol("example-protocol",
                    example_simulator,
                    calc_initial_len=example_calc_inital_len,
                    calc_extend=example_calc_extend)
