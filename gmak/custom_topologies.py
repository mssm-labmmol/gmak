from gmak.topology import add_custom_topology

# *TODO*
# 
# *Document usage of add_custom_topology*
# *Document correspondence in input_file*

# Example 1
def getfiles_example_1(prefix):
    #return [prefix + ".mytop1", prefix + ".mytop2"]
    return prefix + ".mytop1"

def writer_example_1(pars_dict, args_dict, output_files):
    #for fn in output_files:
    fn = output_files
    fp = open(fn, "w")
    fp.write("This is parameters:\n")
    fp.write(str(pars_dict))
    fp.write("This is args:\n")
    fp.write(str(args_dict))
    fp.write("Macro test (work): @custom\n")
    fp.write("Macro test (not work): @acustom\n")
    fp.close()

add_custom_topology("mytop", writer_example_1, getfiles_example_1)
