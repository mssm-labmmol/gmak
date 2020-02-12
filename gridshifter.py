import os
from parameters import * 
from RunFromInput import GridPoint

# creates (soft) symbolic link from target to link_name
def create_symbolic_link (target, link_name):
    if not (os.path.isdir(link_name)):
        print ("Linking: ln -s %s %s" % (target, link_name))
        os.system ("ln -s %s %s" % (target, link_name))
    else:
        print ("Warning: Link %s already exists.\n" % link_name)

class GridShifter:

    def __init__ (self):

        # Margins delimiting the border region.  self.margins = {}
        self.maxshifts = 0
        self.nshifts   = 0
        self.margins = {}
        self.margins['left'] = 0.10
        self.margins['right'] = 0.90
        self.margins['bottom'] = 0.10
        self.margins['top'] = 0.90

        # Fraction of the best points which, if inside the border region,
        # trigger the shifting of the grid.
        self.ncut = 0.10

        # Stores the CG.
        self.cg = (0,0)

    # Calculates CG of the best ncut percent points.
    # Returns it as a tuple.
    # Store in internal variables.
    def calcCG (self, optimizer, grid):
        
        if (len(grid.size) != 2):
            print "Error: GridShifter expects a bi-dimensional grid."
            exit()

        nx = grid.size[0]
        ny = grid.size[1]

        cg_x = 0
        cg_y = 0
        n    = 0
        thr  = int (self.ncut * grid.linear_size)

        for k in range(thr):
            idx = optimizer.stateScores[k][0]
            i   = int(idx / ny)
            j   = idx % ny
            cg_x += i
            cg_y += j

        cg_x /= thr
        cg_y /= thr

        self.cg = (cg_x, cg_y)
        print "Note: CG is %d, %d" % (cg_x, cg_y)

        return (cg_x, cg_y)

    # With CG calculated, verify if we need to shift.
    def doWeShift (self, grid):
        if (self.nshifts == self.maxshifts):
            return False
        x = self.cg[0]
        y = self.cg[1]
        real_left = int(self.margins['left'] * grid.size[0])
        real_right = int(self.margins['right'] * grid.size[0])
        real_top = int(self.margins['top'] * grid.size[1])
        real_bottom = int(self.margins['bottom'] * grid.size[1])
        if (x <= real_left) or (x >= real_right):
            return True
        if (y <= real_bottom) or (y >= real_top):
            return True
        return False

    # Performs shifting operations.
    def shift (self, grid, paramLoop, old_workdir, new_workdir):

        # Create directory of new grid.
        os.system("mkdir -p " + new_workdir)

        # Move paramLoop.
        paramLoop.loop.set_new_center (self.cg[0]*grid.size[1] + self.cg[1])

        # Create a grid file containing the parameters.
        paramLoop.loop.write_parameter_grid_to_file (new_workdir + "/grid.dat")

        # Clean MBAR results.
        grid.hashOfMBAR = {}

        # New gridpoints.
        new_gridpoints = [0] * grid.linear_size

        for gp in grid.grid_points:
            # Clean properties and reweight outputs.
            gp.rw_outputs = {}
            gp.estimated_properties = {}
            gp.atomic_properties = {}
            i_center = int((grid.size[0] - 1)/2)
            j_center = int((grid.size[1] - 1)/2)
            i = int(gp.id / grid.size[1])
            j = gp.id % grid.size[1]
            new_i = i - (self.cg[0] - i_center)
            new_j = j - (self.cg[1] - j_center)
            new_id = new_i * grid.size[1] + new_j 
            # Link from the directories of the old simulation results to 
            # new simulation results. There is no need to change the 
            # itp path, I think.
            if (new_i >= 0) and (new_i < grid.size[0]) and (new_j >= 0) and (new_j < grid.size[1]):
                # loop over protocols
                for prot in gp.protocol_outputs:
                    # create a link for each protocol
                    old_dir = old_workdir + ("/%s/simu/%d" % (prot,gp.id))
                    new_dir = new_workdir + ("/%s/simu/%d" % (prot,new_id))
                    os.system("mkdir -p " + new_workdir + ("/%s/simu/" % prot))
                    create_symbolic_link(old_dir, new_dir)
                # after linking, I can now change the id 
                gp.id = new_id
                # and put in the new list
                new_gridpoints[gp.id] = gp
        
        # For the rest of the list, create gridpoints.
        for i in range(grid.linear_size):
            if (new_gridpoints[i] == 0):
                itp_path = new_workdir + "grid.dat"
                new_gridpoints[i] = GridPoint(itp_path, i)
        
        # Substitute the grid's gridpoints.
        grid.grid_points = new_gridpoints

        # Now everything should be working.
        return


    # Apply shifter, which means: shift if necessary, do nothing if not.
    def apply (self, optimizer, paramLoop, grid, old_workdir, new_workdir):
        self.calcCG(optimizer, grid)
        if not (self.doWeShift(grid)):
            print ("Note: No need to shift the grid.")
            return False
        print ("Note: We are shifting the grid.")
        self.shift(grid, paramLoop, old_workdir, new_workdir)
        self.nshifts += 1
        return True


    def readFromStream (self, stream):
        for line in stream:
            if line[0] == '#':
                continue
            if (re.match(r"^\$end.*",line)):
                return
            splitted = line.split()
            if (splitted[0] == 'margins'):
                self.margins['left'] = float(splitted[1])
                self.margins['right'] = float(splitted[2])
                self.margins['bottom'] = float(splitted[3])
                self.margins['top'] = float(splitted[4])
            if (splitted[0] == 'ncut'):
                self.ncut = float(splitted[1])
            if (splitted[0] == 'maxshifts'):
                self.maxshifts = int(splitted[1])


