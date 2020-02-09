import os

def create_symbolic_link (target, link_name):
    os.system ("ln -s %s %s" % (target, link_name))

class GridShifter:

    def __init__ (self):

        # Margins delimiting the border region.
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
            idx = self.stateScores[k][0]
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
    def shift (self, grid, old_workdir, new_workdir):

        # Create directory of new grid.
        os.mkdir(new_workdir)

        # Create a grid file containing the parameters.
        # TODO

        # Clean MBAR results.
        grid.hashOfMBAR = {}

        # New gridpoints.
        new_gridpoints = [0] * grid.linear_size

        for gp in grid.grid_points:
            # Clean properties and reweight outputs.
            gp.rw_outputs = {}
            gp.estimated_properties = {}
            gp.atomic_properties = {}
            i = int(gp.id / grid.size[1])
            j = gp.id % grid.size[1]
            new_i = i - self.cg[0]
            new_j = j - self.cg[1]
            new_id = new_i * grid.size[1] + new_j 
            # Link from the directories of the old simulation results to 
            # new simulation results. There is no need to change the 
            # itp path, I think.
            if (new_i >= 0) and (new_j >= 0):
                # loop over protocols
                for prot in gp.protocol_outputs:
                    # create a link for each protocol
                    old_dir = old_workdir + ("/%s/simu/%d" % (prot,gp.id))
                    new_dir = new_workdir + ("/%s/simu/%d" % (prot,new_id))
                    os.mkdir(new_workdir +  ("/%s/simu/" % prot))
                    create_symbolic_link(old_dir, new_dir)
                # after linking, I can now change the id 
                gp.id = new_id
                # and put in the new list
                new_gridpoints[gp.id] = gp

        
        # For the rest of the list, create gridpoints.
        for i in range(grid.linear_size):
            if (new_gridpoints[i] == 0):
                # TODO Here I should create the itp file based on the force
                # field loop
                new_gridpoints[i] = GridPoint(itp_path, i)
        
        # Substitute the grid's gridpoints.
        grid.grid_points = new_gridpoints


    # Apply shifter, which means: shift if necessary, do nothing if not.
    def apply (self, optimizer, grid, old_workdir, new_workdir):
        self.calcCG(optimizer, grid)
        if not (doWeShift(grid)):
            print ("Note: No need to shift the grid.")
            return
        self.shift(grid, old_workdir, new_workdir)


