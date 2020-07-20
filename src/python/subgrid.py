from  copy             import  deepcopy
from  grid_ana         import  *
from  property         import  *
from  cartesiangrid    import  *
from  gridbase         import  *
from  surrogate_model  import  *

def linearMaskSubgridToGrid(subgrid, grid):
    mask = subgrid.getVolume()
    lens = subgrid.getLens()
    lc   = 0
    for l,t in zip(CartesianGridLinearIterator(subgrid), CartesianGridIterator(subgrid)):
        success = True
        originalPosition = []
        for xi, li in zip(t, lens):
            if (xi % li == 0):
                originalPosition.append(int(xi/li))
            else:
                success = False
                break
        if (sucess):
            originalPosition = tuple(originalPosition)
            originalLinearPosition = grid.tuple2linear(originalPosition)
            mask[lc] = originalLinearPosition
        else:
            mask[lc] = -1
        lc += 1

    return mask

class SubgridPoint:

    def __init__(originalGridPoint, index):
        self.id = index
        self.is_sample = False
        self.estimated_properties = {}
        if (originalGridPoint) is not None:
            self.is_sample            = originalGridPoint.is_sample
            self.estimated_properties = deepcopy(originalGridPoint.estimated_properties)

    def get_property_estimate(self, prop):
        return self.estimated_properties[prop].value

    def get_property_err(self, prop):
        return self.estimated_properties[prop].err

    # obj is a PropertyBase::X object (X=Density,Gamma,dHvap)
    def add_property_estimate (self, prop_id, prop_name, obj):
        if (prop_name != obj.name):
            raise ValueError ("ERROR: expected {} but got {}.\n".format(prop_name, obj.name))
        self.estimated_properties[prop_id] = obj

class ParameterSubgrid:

    def __init__(self, parSpaceGen, originalParameterGrid, propid2type, interpModelString, boolLegacy):
        self.parSpaceGen  = parSpaceGen
        self.indexGrid    = CartesianGrid(self.get_size())
        self.originalGrid = originalParameterGrid
        self.model        = init_model_from_string(interpModelString, boolLegacy)
        self.xlabel       = originalParameterGrid.xlabel
        self.ylabel       = originalParameterGrid.ylabel

        _indexMask = linearMaskSubgridToGrid(self.indexGrid, self.originalGrid)
        # Initialize grid_points
        for i in range(self.get_linear_size()):
            originalIndex = _indexMask[i]
            if (originalIndex == -1):
                self.grid_points.append(SubgridPoint(None, i))
            else:
                self.grid_points.append(SubgridPoint(self.originalGrid[originalIndex], i))                 
        # Important warning
        if (self.get_dim() != 2):
            raise ValueError("With the current implementation, only bi-dimensional subgrids are allowed.")

        # Estimates for new gridpoints
        # Re-build property matrix for estimation
        I_e                    = self.get_samples_id()
        list_of_property_names = list(self.grid_points[I_e[0]].estimated_properties.keys())
        A_pe                   = [None] * len(list_of_property_names)
        dA_pe                  = [None] * len(list_of_property_names)        
        for p, prop_id in enumerate(list_of_property_names):
            A_pe[p]  = [None] * originalGrid.get_linear_size()
            dA_pe[p] = [None] * originalGrid.get_linear_size()
            for e in range(self.linear_size):
                A_pe[p][e]  = self.grid_points[I_e[e]].get_property_estimate(prop_id)
                dA_pe[p][e] = self.grid_points[I_e[e]].get_property_err(prop_id)
        A_pe            = np.array(A_pe)
        dA_pe           = np.array(dA_pe)
        (EA_pk, dEA_pk) = model.computeExpectationsFromAvgStd(A_pe, dA_pe, I_e, tuple(self.get_size()))
        for p, prop_id in enumerate(list_of_property_names):
            for gp in self.grid_points:
                propertyObj = init_property_from_string(propid2type[prop_id], EA_pk[p, gp.id], dEA_pk[p, gp.id])
                gp.add_property_estimate(prop_id, propid2type[prop_id], propertyObj)

    def get_size(self):
        return self.parSpaceGen.getSizes(self._mainString)

    def get_linear_size(self):
        return self.parSpaceGen.getNumberOfStates()

    def get_samples_id (self):
        sp = []
        for x in self.grid_points:
            if x.is_sample:
                sp.append(x.id)
        return sp
        
    # Output functions 
    def save_property_values_to_file (self, prop, filename):
        # make data
        data = [x.estimated_properties[prop].value for x in self.grid_points]
        data = np.array(data)
        np.savetxt(filename, data)

    def read_property_values_err_from_file (self, prop, propType, filename_value, filename_err):
        values = np.loadtxt(filename_value)
        err = np.loadtxt(filename_err)
        for i,x in enumerate(self.grid_points):
            x.estimated_properties[prop] = init_property_from_string(propType, values[i], err[i])

    def save_property_err_to_file (self, prop, filename):
        # make data
        data = [x.estimated_properties[prop].err for x in self.grid_points]
        data = np.array(data)
        np.savetxt(filename, data)

    def save_property_diff_to_file (self, prop, ref, filename):
        # make data
        data = [x.estimated_properties[prop].value - ref for x in self.grid_points]
        data = np.array(data)
        np.savetxt(filename, data)

    def plot_property_to_file (self, prop, filename):
        if (self.get_dim() != 2):
            warnings.warn("Can only plot 2-D grids.")
            return
        # make data
        data = [x.estimated_properties[prop].value for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.get_size())
        cbox_label = ""
        cbox_limits_colors = ()
        cbox_limits = ()
        title = self.grid_points[0].estimated_properties[prop].get_label()
        # plot
        # from grid_ana ...
        plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())

    def plot_property_err_to_file (self, prop, filename):
        if (self.get_dim() != 2):
            warnings.warn("Can only plot 2-D grids.")
            return
        # make data
        data = [x.estimated_properties[prop].err for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.get_size())
        cbox_label = ""
        cbox_limits_colors = ('white', 'blue')
        cbox_limits = (0.0, np.max(data))
        title = self.grid_points[0].estimated_properties[prop].get_label_err()
        # plot
        # from grid_ana ...
        plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())

    def plot_property_diff_to_file (self, prop, propname, ref, filename):
        if (self.get_dim() != 2):
            warnings.warn("Can only plot 2-D grids.")
            return
        # make data
        data = [x.estimated_properties[prop].value - ref for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.get_size())
        cbox_label = ""
        # hard-coded preferences
        cbox_limits_colors = ('red', 'white', 'blue')
        # hard-coded preferences
        limits = {\
                'density': 20,\
                'dhvap': 4,\
                  'gamma': 10, \
                  'ced': 10 * 1e-4, \
                  'gced': 10}
        cbox_limits = (-limits[propname], limits[propname])
        title = "$\\Delta^\\mathrm{exp}$" + self.grid_points[0].estimated_properties[prop].get_label()
        # plot
        # from grid_ana ...
        plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())
