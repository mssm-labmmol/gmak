from math import sqrt

class Atomtype:

    def __init__ (self):
        self.type = ""
        self.parameters = {'c6': 0.00, 'c12': 0.00, 'cs6': 0.00, 'cs12': 0.00}

    def __init__ (self, typestring):
        self.type = typestring
        self.parameters = {'c6': 0.00, 'c12': 0.00, 'cs6': 0.00, 'cs12': 0.00}

    def set_parameters (self, c6, c12, cs6, cs12):
        self.parameters['c6'] = c6
        self.parameters['c12'] = c12
        self.parameters['cs6'] = cs6
        self.parameters['cs12'] = cs12

    def alter_parameter (self, parameterName, newValue):
        self.parameters[parameterName] = newValue

    def decouple (self, parameterName):
        self.alter_parameter (parameterName, 0.0)

class LoopData:
    
    def __init__ (self):
        self.atomtype = Atomtype ()
        self.start['c6']  = 0.00
        self.step['c6']   = 0.00
        self.size['c6']   = 0
        self.start['c12'] = 0.00
        self.step['c12']  = 0.00
        self.size['c12']  = 0

    def set_members (self, atomtype, c6start, c6step, c6size, c12start, c12step, c12size):
        # this is an Atomtype instance
        self.atomtype = atomtype
        self.start['c6'] = c6start
        self.step['c6']  = c6step
        self.size['c6']  = c6size
        self.start['c12'] = c12start
        self.step['c12']  = c12step
        self.size['c12']  = c12size

    def set_at_linear_position (self, position):
        c6_position = int(position / self.size['c12'])
        c12_position = position % self.size['c12']
        this_c6 = self.start['c6'] + c6_position * self.step['c6']
        this_c12 = self.start['c12'] + c12_position * self.step['c12']
        self.atomtype.alter_parameter ('c6', this_c6)
        self.atomtype.alter_parameter ('c12', this_c12)

class ParameterLoop:

    def __init__ (self):
        self.basic_itp = "/dev/null"
        # list of Atomtypes
        self.atomtypes = []
        # loop data
        self.loop = LoopData ()

    # Writes FF parameters to stream.
    # Decouple flag:
    #
    #       0: no decoupling
    #       1: decouple c12 of self-interactions
    #       2: 1 + decouple sqrt(c12) in non-self-interactions
    #       3: 2 + decouple c6 of self-interactions
    #       4: 3 + decouple sqrt(c6) in non-self-interactions
    #
    # The modified atomtype is inferred from self.loop.atomtype.
    def write_ff_header_to_stream (self, stream, decoupleFlag):
        stream.write("[ defaults ]\n")
        stream.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
        stream.write("       1                1              no           1.0     1.0\n")
        stream.write("[ atomtypes ]\n")
        for at in self.atomtypes:
            selfc6  = at.parameters['c6']
            selfc12 = at.parameters['c12']
            if (at == self.loop.atomtype):
                if (decoupleFlag > 0):
                    selfc12 = 0
                if (decoupleFlag > 2):
                    selfc6   = 0
            stream.write("%-20s%10d%10.2f%10.2f%18.7e%18.7e\n", at.type, 6, 0.0, 0.0, selfc6, selfc12)
        stream.write("\n[ nonbond_params ]\n")
        for at in self.atomtypes:
            selfc6  = at.parameters['c6']
            selfc12 = at.parameters['c12']
            if (at == self.loop.atomtype):
                if (decoupleFlag > 1):
                    selfc12 = 0
                if (decoupleFlag > 3):
                    selfc6   = 0
            for ato in self.atomtypes:
                if (ato != at):
                    selfc6_o  = ato.parameters['c6']
                    selfc12_o = ato.parameters['c12']
                    c6 = sqrt(selfc6 * selfc6_o)
                    c12 = sqrt(selfc12 * selfc12_o)
                    stream.write("%-20s%-20s%10d%18.7e%18.7e\n", at.type, ato.type, 1, c6, c12)
        stream.write("\n[ pairtypes ]\n")
        for at in self.atomtypes:
            selfc6  = at.parameters['cs6']
            selfc12 = at.parameters['cs12']
            for ato in self.atomtypes:
                if (ato != at):
                    selfc6_o  = ato.parameters['cs6']
                    selfc12_o = ato.parameters['cs12']
                    c6 = sqrt(selfc6 * selfc6_o)
                    c12 = sqrt(selfc12 * selfc12_o)
                    stream.write("%-20s%-20s%10d%18.7e%18.7e\n", at.type, ato.type, 1, c6, c12)
        stream.write("\n")

    # create a itp with the given decoupleFlag
    def create_itp (self, output, decoupleFlag):
        fp = open(output, 'w')
        self.write_ff_header_to_stream (fp, decoupleFlag)
        # also copy itp file into it
        with open(self.basic_itp, 'r') as fpi:
            fp.write(fpi.read())
        fp.close()

    # Move loop to position and create an itp file based on this position.
    def create_full_itp_for_position (self, linear_position, output_preffix):
        self.loop.set_at_linear_position (linear_position)
        self.create_itp (output_fn)

    # This assumes a stream is given.
    #
    # It will read from stream until line with terminating string '$end' is 
    # found.
    #
    def read_from_stream (self, stream):
        starts = []
        steps  = []
        sizes  = []
        pars   = []
        for line in stream:
            if line[0] == '#':
                continue
            if (re.match(r"^\$end.*",line)):
                # before ending, set loop members
                idx_c6 = pars.index('c6')
                idx_c12 = pars.index('c12')
                self.loop.set_members (self.atomtypes[index_of_atomtype], \
                        starts[idx_c6], steps[idx_c6], sizes[idx_c6],\
                        starts[idx_c12], steps[idx_c12], sizes[idx_c12])
                return
            splitted = line.split()
            if (line.split()[0] == 'itp'):
                self.basic_itp = splitted[1]
            if (line.split()[0] == 'make_top'):
                self.make_top_path = os.path.abspath(splitted[1])
            if (line.split()[0] == 'atomtypes'):
                for at in splitted[1:]:
                    new_at = Atomtype (at)
                    self.atomtypes.append(new_at)
            if (splitted[0] in ['c6','c12','cs6','cs12']):
                for i,par in enumerate(splitted[1:]):
                    self.atomtypes[i].alter_parameter(splitted[0], float(par))
            if (splitted[0] == 'loop'):
                # first, set the reference to the atomtype in the loop
                dummy_atomtype     = Atomtype (splitted[1])
                index_of_atomtype  = self.atomtypes.index(dummy_atomtype)
                # now set pars 
                pars = splitted[2:]
            if (splitted[0] == 'start'):
                starts = float(splitted[1:])
            if (splitted[0] == 'step'):
                steps = float(splitted[1:])
            if (splitted[0] == 'size'):
                size = int(splitted[1:])

