import os

class InputModifier:

    def __init__(self, src):
        self.src = src
        # Things that may be altered
        # If None, they won't be altered; otherwise, they will.
        self.workdir = None
        self.main_variation = None
        self.samples = None

    def set_samples(self, linear_ids):
        self.samples = linear_ids

    def set_workdir(self, wd):
        self.workdir = wd

    def set_main_variation(self, par_space_gen):
        self.main_variation = par_space_gen

    def _checkEndOfBlock(self, line):
        return (line.rstrip() == "$end")
        
    def write_to_file(self, dest):
        fpi = open(self.src, 'r')
        fpo = open(dest, 'w')
        line = fpi.readline()
        while line:
            if (len(line.split()) < 1):
                fpo.write(line)
                line = fpi.readline()
                continue
            # test workdir
            if (line.split()[0] == 'workdir'):
                if (self.workdir is not None):
                    fpo.write("workdir " + self.workdir + "\n")
                else:
                    fpo.write(line)

                line = fpi.readline()
                continue
            # test variation
            elif (line.split()[0] == '$variation'):
                if (self.main_variation is not None):
                    # store current position
                    start_pos = fpi.tell()
                    # check if variation is main
                    for line in fpi:
                        flds = line.split()
                        if (flds[0] == 'name'):
                            if (flds[1] == 'main'):
                                self.main_variation.writeMainVariationBlock(fpo)
                                for line in fpi:
                                    if self._checkEndOfBlock(line):
                                        break
                            else:
                                fpi.seek(start_pos, os.SEEK_SET)
                                for line in fpi:
                                    fpo.write(line)
                                    if self._checkEndOfBlock(line):
                                        break
                        if self._checkEndOfBlock(line):
                            break
                else:
                    for line in fpi:
                        fpo.write(line)
                        if self._checkEndOfBlock(line):
                            break
                line = fpi.readline()
                continue
            # test samples
            elif (line.split()[0] == 'samples'):
                if (self.samples is not None):
                    fpo.write("samples ")
                    fpo.write(" ".join(map(str, self.samples)))
                    fpo.write("\n")
                else:
                    fpo.write(line)
                line = fpi.readline()
                continue
            # other cases
            else:
                fpo.write(line)
                line = fpi.readline()

        fpo.close()
        fpi.close()
