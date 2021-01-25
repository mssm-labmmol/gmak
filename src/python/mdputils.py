import runcmd

class mdpUtils:
    """A simple class to read info and possibly manipulate GROMACS .mdp files.
    Features are implemented as they are needed, so many functionalities may still
    be absent."""

    def __init__(self):
        """Default initializer."""
        self.filename = ""
        self.optionsDict = {}
        return

    def clean(self):
        """Clean contents."""
        self.filename = ""
        self.optionsDict = {}
        return

    def parse_file(self, fn, commchar=';'):
        """Parses a .mdp file, capturing the values for each option."""
        self.clean()
        self.filename = fn
        fp = open(fn, 'r')
        for line in fp:
            clean_line = line.split(commchar)[0].lstrip().rstrip()
            if (clean_line != ''):
                splitted_line = clean_line.split()
                option = splitted_line[0]
                values = splitted_line[2:]
                self.optionsDict[option] = values
        fp.close()

    def get_option(self, option):
        """Returns a list with the values of the option."""
        try:
            return self.optionsDict[option]
        except KeyError:
            raise KeyError("Option {} was not specified in file {}.".format(option, self.filename))

    def get_nsteps(self):
        """Return nsteps."""
        return int(self.get_option('nsteps')[0])

    def get_writefreq_energy(self):
        """Returns frequency (in number of steps) of writing energies.
         This also checks if the frequency is a multiple of the frequency of calculating energy."""
        energy_freq = int(self.get_option('nstenergy')[0])
        calc_freq = int(self.get_option('nstcalcenergy')[0])
        if ((energy_freq % calc_freq) == 0):
            return energy_freq
        else:
            raise ValueError("Frequency of writing energies is not a multiple of the frequency of calculating energies.")

    def get_writefreq_pos(self):
        """Returns frequency (in number of steps) of writing positions to .trr file."""
        return int(self.get_option('nstxout')[0])

    def get_writefreq_compressedpos(self):
        """Returns frequency (in number of steps) of writing positions to .xtc file."""
        return int(self.get_option('nstxtcout')[0])

    def get_reference_temperatures(self):
        """Returns a list of reference temperatures for each temperature bath."""
        try:
            string_list = self.get_option('ref_t')
        except KeyError:
            string_list = self.get_option('ref-t')
        float_list  = [float(x) for x in string_list]
        return float_list

    def get_nlambdas(self):
        """Returns the number of lambda values for a free energy mdp file."""
        # make a list of all *-lambdas options
        lambdaOptions = [k for k in self.optionsDict.keys() if k.endswith('-lambdas')]
        if len(lambdaOptions) > 0:
            # get number of lambda for each
            nlambdas = [len(self.get_option(k)) for k in lambdaOptions]
            nlambdasCte = nlambdas[0]
            for l in nlambdas:
               if (l != nlambdasCte):
                   raise Exception('Two *-lambdas options have a different number of lambda values.')
            return nlambdasCte
        else:
            raise Exception('Requesting nlambdas for a mdp file that does not contain *-lambdas options.')

    def set_option(self, option, value_list):
        self.optionsDict[option] = value_list
        
    def set_lambda_state(self, value):
        option = 'init-lambda-state'
        self.optionsDict[option] = [value]

    def write_to_file(self, fn):
        fp = open(fn, 'w')
        for k in self.optionsDict.keys():
            fp.write("%-16s = " % k)
            for v in self.optionsDict[k]:
                fp.write("%-10s" % v)
            fp.write("\n")
        fp.close()
