#!/usr/bin/env python3


class ParetoFront:
    import warnings

    """
    ------------------------------------------------------------------
    Auxiliary methods.
    ------------------------------------------------------------------
    """

    @staticmethod
    def _x_dominates_y(fx, fy):
        """
        Auxiliary function. Returns True if x dominates y, False
        otherwise.
        """
        eqflag = True
        for fx_i, fy_i in zip(fx, fy):
            if (fx_i > fy_i):
                return False
            elif (eqflag) and (fx_i != fy_i):
                eqflag = False
        if eqflag:
            return False
        return True


    def _prepare_front_standard(self):
        self.pf = []
        self.fp = []
        for xi, fxi in zip(self.x, self.fx):
            in_pareto = True
            for xj, fxj in zip(self.x, self.fx):
                if ParetoFront._x_dominates_y(fxj, fxi):
                    in_pareto = False
                    break
            if (in_pareto):
                self.pf.append(xi)
                self.fp.append(fxi)


    def _check(self):
        if self.get_dim_f() < 2:
            warnings.warn("Pareto front only makes sense for "
                          "multi-objective problems.")
            return False
        return True


    """
    ------------------------------------------------------------------
    Factories.
    ------------------------------------------------------------------
    """

    @staticmethod
    def read_from_csv(csv_fn, xdim, comments='#'):
        import numpy as np
        data = np.loadtxt(csv_fn, comments=[comments,])
        x = data[:, :xdim]
        f = data[:, xdim:]
        return ParetoFront(x, func_data=f)


    @staticmethod
    def read_from_gmak_step(x_fn, diff_k_fn_list):
        """
        Parameters:
          x_fn: file containing parameters
          diff_k_fn_list: list of diff_k files. Only absolute values
                          will be taken into account.
        """
        import numpy as np
        x = np.loadtxt(x_fn)
        nprops = len(diff_k_fn_list)
        npoints = len(x)
        f = np.zeros((npoints, nprops))
        for p, diff_k in enumerate(diff_k_fn_list):
            f[:,p] = np.abs(np.loadtxt(diff_k))
        return ParetoFront(x.tolist(), func_data=f.tolist())

    @staticmethod
    def read_from_gmak(rootdir):
        from glob import glob
        import numpy as np
        pfs = []
        for grid in glob(rootdir + "/grid_*/"):
            steps = glob(grid + "/step_*/")
            steps.sort()
            last_step = steps[-1]
            pf = ParetoFront.read_from_gmak_step(
                grid + "/parameters_main.dat",
                glob(last_step + "/*diff*"))
            pfs.append(pf)
        return ParetoFront.join(pfs)

    """
    ------------------------------------------------------------------
    Public methods.
    ------------------------------------------------------------------
    """

    def __init__(self, x, func_data=None, func=None):
        self.x = x
        if func is not None:
            self.fx = [func(xx) for xx in x]
        elif func_data is not None:
            self.fx = func_data.copy()
        else:
            raise ValueError("ParetoFront: either func_data or func "
                             "must not be None.")

        if not self._check():
            return None

        # Populates self.pf and self.fp.
        self._prepare_front_standard()


    @staticmethod
    def join(pf_list):
        x = []
        fx = []
        for pf in pf_list:
            x += pf.x
            fx += pf.fx
        return ParetoFront(x, func_data=fx)


    def get_dim_x(self):
        try:
            return len(self.x[0])
        except TypeError:
            return 1


    def get_dim_f(self):
        try:
            return len(self.fx[0])
        except TypeError:
            return 1


    def size(self):
        return len(self.pf)


    def plot(self, fn=None, dump=None):
        """
        Plots a Pareto front together with other samples in
        objective-function space.

        The plot is saved in file <fn>. If <fn> is none, the plot is
        shown for interaction.

        The figure instance is saved to file <dump>, if not None.

        pf is the Pareto front, a list of parameter-space points.

        x is an iterable of evaluated parameter-space points.
        """
        import matplotlib.pyplot as plt
        import pickle

        dim = self.get_dim_f()
        if dim > 3:
            raise ValueError("Can't plot more than three dimensions.")
        else:
            if dim == 2:
                # 2D plot
                fx_1 = [f[0] for f in self.fx]
                fx_2 = [f[1] for f in self.fx]
                fp_1 = [f[0] for f in self.fp]
                fp_2 = [f[1] for f in self.fp]
                fig, axs = plt.subplots()
                ax.scatter(fx_1, fx_2)
                ax.scatter(fp_1, fp_2)
            elif dim == 3:
                # 3D plot
                fx_1 = [f[0] for f in self.fx]
                fx_2 = [f[1] for f in self.fx]
                fx_3 = [f[2] for f in self.fx]
                fp_1 = [f[0] for f in self.fp]
                fp_2 = [f[1] for f in self.fp]
                fp_3 = [f[2] for f in self.fp]
                fig = plt.figure()
                ax  = plt.axes(projection='3d')
                ax.scatter(fx_1, fx_2, fx_3)
                ax.scatter(fp_1, fp_2, fp_3)
        if fn is None:
            plt.show()
        else:
            plt.savefig(fn)
        if dump is not None:
            pickle.dump(fig, open(dump,'wb'))


    def write(self, fn):
        """Saves the Pareto front to file <fn>."""
        import numpy as np
        xdim = self.get_dim_x()
        fdim = self.get_dim_f()
        paretolen = self.size()
        all = np.zeros((paretolen, xdim + fdim))
        for i, (xx, fp) in enumerate(zip(self.pf, self.fp)):
            if xdim > 1:
                all[i,:len(xx)] = xx
                all[i,len(xx):] = fp
            else:
                all[i,:1] = xx
                all[i,1:] = fp
        np.savetxt(fn, all)


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        print("Usage: {} <project_directory> <output_prefix>".format(sys.argv[0]))
        sys.exit()

    paretoFront = ParetoFront.read_from_gmak(sys.argv[1])
    paretoFront.write(sys.argv[2] + ".dat")
    paretoFront.plot()
    paretoFront.plot(fn=sys.argv[2] + ".pdf",
                     dump=sys.argv[2] + ".fig")
