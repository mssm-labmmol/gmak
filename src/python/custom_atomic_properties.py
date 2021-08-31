from atomic_properties import CustomAtomicPropertyFactory

# Example 1
# ---------
#
# Adding an atomic property that always returns the value 1.0 +/- 0.1 and that
# is named "example_1".
#
# This is supposed to illustrate properties that are not represented by a time
# series, and the user is responsible for calculating the error.
def calculator_example_1(input_dict, out_path):
    import numpy as np
    out = np.array([1.0, 0.1])
    np.savetxt(out_path, out)


CustomAtomicPropertyFactory.add_custom_atomic_property("example_1",
                                                       calculator_example_1,
                                                       is_timeseries=False)

# Example 2
# ---------
#
# Adding an atomic property that emulates the behavior of 'Density'.
def calculator_example_2(input_dict, out_path):
    import numpy as np
    import os
    os.system(f"echo \"Density\" | gmx energy -f {input_dict['edr']} "
              f"-s {input_dict['tpr']} "
              f"-o {out_path}_pre.xvg")
    data = np.loadtxt(f"{out_path}_pre.xvg", comments=['@','#'])
    os.system(f"rm {out_path}_pre.xvg")
    np.savetxt(out_path, data)

CustomAtomicPropertyFactory.add_custom_atomic_property("example_2",
                                                       calculator_example_2,
                                                       is_timeseries=True)
