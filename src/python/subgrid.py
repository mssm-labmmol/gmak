from  copy             import  deepcopy
import runcmd
from  grid_ana         import  *
from  property         import  *
from  cartesiangrid    import  *
from  gridbase         import  *
from  surrogate_model  import  *
from  collections      import  OrderedDict

def linearMaskSubgridToGrid(subgrid, grid):
    mask = [None for i in range(subgrid.getVolume())]
    lens = subgrid.getLens()
    orlens = grid.getLens()
    factors = [int((s-1)/(o-1)) for s,o in zip(lens, orlens)]
    lc   = 0
    for l,t in zip(CartesianGridLinearIterator(subgrid), CartesianGridIterator(subgrid)):
        success = True
        originalPosition = []
        for xi, li in zip(t, factors):
            if (xi % li == 0):
                originalPosition.append(int(xi/li))
            else:
                success = False
                break
        if (success):
            originalPosition = tuple(originalPosition)
            originalLinearPosition = grid.tuple2linear(originalPosition)
            mask[lc] = originalLinearPosition
        else:
            mask[lc] = -1
        lc += 1
    return mask

class SubgridPoint:

    def __init__(self, originalGridPoint, index):
        self.id = index
        self.is_sample = False
        self.is_original_sample = False
        self.estimated_properties = OrderedDict()
        if (originalGridPoint) is not None:
            self.is_sample            = True
            self.is_original_sample   = originalGridPoint.is_sample
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
        self.originalGrid = originalParameterGrid
        self._mainString  = originalParameterGrid._mainString
        self.indexGrid    = CartesianGrid(self.get_size())
        self.model        = init_surrogate_model_from_string(interpModelString, boolLegacy)
        self.xlabel       = originalParameterGrid.xlabel
        self.ylabel       = originalParameterGrid.ylabel
        self.gridWorkdir  = originalParameterGrid.workdir

        _indexMask = linearMaskSubgridToGrid(self.indexGrid, self.originalGrid.getCartesianGrid())
        # Initialize grid_points
        self.grid_points = []
        for i in range(self.get_linear_size()):
            originalIndex = _indexMask[i]
            if (originalIndex == -1):
                self.grid_points.append(SubgridPoint(None, i))
            else:
                self.grid_points.append(SubgridPoint(self.originalGrid[originalIndex], i))                 

        # Estimates for new gridpoints
        # Re-build property matrix for estimation
        I_e                    = self.get_estimated_id()
        list_of_property_names = list(self.grid_points[I_e[0]].estimated_properties.keys())
        A_pe                   = [None] * len(list_of_property_names)
        dA_pe                  = [None] * len(list_of_property_names)
        for p, prop_id in enumerate(list_of_property_names):
            A_pe[p]  = []
            dA_pe[p] = []
            for e,s in enumerate(I_e):
                A_pe[p].append(self.grid_points[s].get_property_estimate(prop_id))
                dA_pe[p].append(self.grid_points[s].get_property_err(prop_id))
        A_pe            = np.array(A_pe)
        dA_pe           = np.array(dA_pe)
        (EA_pk, dEA_pk) = self.model.computeExpectationsFromAvgStd(A_pe, dA_pe, I_e, tuple(self.get_size()))
        for p, prop_id in enumerate(list_of_property_names):
            for gp in self.grid_points:
                propertyObj = init_property_from_string(propid2type[prop_id], EA_pk[p, gp.id], dEA_pk[p, gp.id])
                gp.add_property_estimate(prop_id, propid2type[prop_id], propertyObj)

    def __getitem__(self, i):
        return self.grid_points[i]

    def save_to_binary(self, optimizer):
        # save grid
        fn = self.makeSubgridDir() + "/grid.bin"
        fp = open(fn, 'wb')
        pickle.dump(self, fp, pickle.HIGHEST_PROTOCOL)
        fp.close()
        # save optimizer
        fn = self.makeSubgridDir() + "/optimizer.bin"
        fp = open(fn, 'wb')
        pickle.dump(optimizer, fp, pickle.HIGHEST_PROTOCOL)
        fp.close()
                
    def get_dim(self):
        return self.parSpaceGen.getDimension(self._mainString)
                
    def get_size(self):
        return self.parSpaceGen.getSizes(self._mainString)

    def get_linear_size(self):
        return self.parSpaceGen.getNumberOfStates()

    def get_samples_id (self):
        sp = []
        for x in self.grid_points:
            if x.is_original_sample:
                sp.append(x.id)
        return sp

    def get_estimated_id(self):
        sp = []
        for x in self.grid_points:
            if x.is_sample:
                sp.append(x.id)
        return sp

    def makeSubgridDir(self):
        subgridDir = os.path.abspath(self.gridWorkdir + "/subgrid")
        runcmd.run("mkdir -p " + subgridDir)
        return subgridDir
    
    def save_property_values_to_file (self, prop):
        # make data
        data = [x.estimated_properties[prop].value for x in self.grid_points]
        data = np.array(data)
        filename = "{}/{}_EA_k.dat".format(self.makeSubgridDir(), prop)
        np.savetxt(filename, data)

    def save_property_err_to_file (self, prop):
        # make data
        data = [x.estimated_properties[prop].err for x in self.grid_points]
        data = np.array(data)
        filename = "{}/{}_dEA_k.dat".format(self.makeSubgridDir(), prop)
        np.savetxt(filename, data)

    def save_property_diff_to_file (self, prop, ref):
        # make data
        data = [x.estimated_properties[prop].value - ref for x in self.grid_points]
        data = np.array(data)
        filename = "{}/{}_diff.dat".format(self.makeSubgridDir(), prop)
        np.savetxt(filename, data)

    def plot_property_to_file (self, prop):
        if (self.get_dim() > 2):
            warnings.warn("Can only plot 1-D or 2-D grids.")
            return
        # make data
        data = [x.estimated_properties[prop].value for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.get_size())
        cbox_label = ""
        cbox_limits_colors = ()
        cbox_limits = ()
        title = self.grid_points[0].estimated_properties[prop].get_label()
        filename = "{}/{}_EA_k.pdf".format(self.makeSubgridDir(), prop)
        # plot
        # from grid_ana ...
        if (self.get_dim() == 2):
            plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())
        elif (self.get_dim() == 1):
            plot_1d_to_file(filename, title, self.xlabel, data, self.get_samples_id())
        else:
            return

    def plot_property_err_to_file (self, prop):
        if (self.get_dim() > 2):
            warnings.warn("Can only plot 1-D or 2-D grids.")
            return
        # make data
        data = [x.estimated_properties[prop].err for x in self.grid_points]
        data = np.array(data)
        data = data.reshape (self.get_size())
        cbox_label = ""
        cbox_limits_colors = ('white', 'blue')
        cbox_limits = (0.0, np.max(data))
        title = self.grid_points[0].estimated_properties[prop].get_label_err()
        filename = "{}/{}_dEA_k.pdf".format(self.makeSubgridDir(), prop)
        # plot
        # from grid_ana ...
        if (self.get_dim() == 2):
            plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())
        elif (self.get_dim() == 1):
            plot_1d_to_file(filename, title, self.xlabel, data, self.get_samples_id())
        else:
            return

    def plot_property_diff_to_file (self, prop, propname, ref):
        if (self.get_dim() > 2):
            warnings.warn("Can only plot 1-D or 2-D grids.")
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
        filename = "{}/{}_diff.pdf".format(self.makeSubgridDir(), prop)
        # plot
        # from grid_ana ...
        if (self.get_dim() == 2):
            plot_grid_to_file (filename, title, self.xlabel, self.ylabel, cbox_label,\
                    cbox_limits, cbox_limits_colors, data, self.get_samples_id())
        elif (self.get_dim() == 1):
            plot_1d_to_file(filename, title, self.xlabel, data, self.get_samples_id())
        else:
            return

    def printAndPlotScores(self, optimizer, plotFlag):
        optimizer.fillWithScores (self)
        optimizer.printToFile (self, self.makeSubgridDir() + "/optimizer_data.dat")
        optimizer.printToFile (self, self.makeSubgridDir() + "/full_data.dat", sorted=False)
        if (plotFlag):
            optimizer.plotToPdf (self, self.makeSubgridDir() + "/optimizer_score.pdf")

    def writeParameters(self):
        self.parSpaceGen.writeParameters(self.makeSubgridDir() + "/parameters")


