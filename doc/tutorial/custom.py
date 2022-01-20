# import the customization API
from gmak.api import *
# gmak.config.ConfigVariables contains the path of the gmx binary
from gmak.config import ConfigVariables
# these are other useful modules
import os
import tempfile
import re

def msd_calc(topology, protocol_output, property_pars):
    """
    This function plays the role of the
    gmak.api_signatures.component_calculator() function
    (please consult the customization API documentation).
    
    It receives a topology file, the protocol-output dictionary with
    the simulation data and the input parameters (property_pars)
    passed to the program in the $compute block where the property is
    used.

    It returns a tuple with the expected value of the self-diffusion
    coefficient and the error in the estimate.
    """
    # path of the gmx binary
    gmx = ConfigVariables.gmx
    # tpr file of the simulation
    tpr = protocol_output['tpr']
    # xtc file of the simulation
    xtc = protocol_output['xtc']
    # create two temporary files to store the output of the ~gmx msd~ command
    d_out = tempfile.NamedTemporaryFile()
    msd_out = tempfile.NamedTemporaryFile(suffix='.xvg')
    # extract the path of the temporary files
    msd_out_path = msd_out.name
    d_out_path = d_out.name
    try:
        # try to retrieve the path of the index file and the name
        # of the group used as reference for calculating the
        # self-diffusion coefficient
        #
        # for the input file in this tutorial, this should work
        # since we supply the input parameters ~index_file~ and
        # ~group_name~ in the input file
        ndx = property_pars.index_file
        group = property_pars.group_name
        # finally, call the ~gmx msd~ command passing the group name
        # and the index file and storing the output in the temporary
        # files created above
        os.system(f"echo {group} | {gmx}  msd -n {ndx} -f {xtc} -s {tpr} -o {msd_out_path} > {d_out_path}")
    except AttributeError:
        # if retrieving the group name or the index file fails, use
        # all atoms as reference
        os.system(f"echo 0 | {gmx}  msd -f {xtc} -s {tpr} -o {msd_out_path} > {d_out_path}")
    # open the ~gmx msd~ output file search for the line containing
    # the value of the self-diffusion coefficient and the fitting
    # error
    with open(d_out_path, 'r') as fp:
        for line in fp:
            m = re.match('^D\[.*\]\s+(\S+)\s+\(\+/-\s+(\S+)\)', line)
            # the line that matches the regex above has the value of
            # the self-diffusion coefficient in the first group and
            # the fitting error in the second group
            if m:
                value = float(m.group(1))
                err = float(m.group(2))
    # close (thus removing) the temporary files
    d_out.close()
    msd_out.close()
    # return the value and the fitting error
    return (value, err)

# see the documentation for gmak.api.add_custom_component_property()
#
# add ~msd~ as a custom component property associated with the 
# ~msd_calc~ function defined above
#
# note that it does not correspond to a timeseries
add_custom_component_property("msd",
                              msd_calc,
                              is_timeseries=False)

# see the documentation for gmak.api.add_custom_composite_property()
#
# add ~diffusion_coeff~ as a custom composite property that is
# recognized by the program
#
# since a calculator function is not supplied, it is implicitly
# assumed that this property has only one component and that the value
# and error correspond to those of the component property
#
# the association between ~diffusion_coeff~ and the ~msd~ component
# property is established in the input file
add_custom_composite_property("diffusion_coeff")

