#!/usr/bin/env python3

"""
TODO:

* Allow passing a <*.csv> file instead of the project directory as the
  first argument.
"""

class ParetoFront:

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
            raise Exception("Pareto front only makes sense for "
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
    def join(pf_list, remove_duplicates=True):
        import numpy as np
        import warnings
        x = []
        fx = []
        dim_f = pf_list[0].get_dim_f()
        for pf in pf_list:
            x += pf.x
            fx += pf.fx
        if (remove_duplicates):
            for ix, xi in enumerate(x):
                # get all indexes corresponding to xi
                ii = [i for i,_ in enumerate(x) if (x[i] == xi) and (i != ix)]
                # get the corresponding function values
                fs = [fx[ix],] + [fx[i] for i in ii]
                # do statistics
                avg = np.mean(fs, axis=0)
                if len(ii) > 0:
                    warnings.warn(f"{len(ii)+1} occurrences of {xi} were replaced "
                                  f"by their average function value: {avg}.")
                # replace all occurrences by a single occurence with
                # the average values
                for index in sorted(ii, reverse=True):
                    del x[index]
                    del fx[index]
                fx[ix] = avg
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


    def _plot_core(self, data_1, data_2, dim,
                   labels, fn=None, dump=None):
        import matplotlib.pyplot as plt
        import pickle

        if dim > 3:
            raise ValueError("Can't plot more than three dimensions.")
        else:
            if dim == 2:
                # 2D plot
                fx_1 = [f[0] for f in data_1]
                fx_2 = [f[1] for f in data_1]
                fp_1 = [f[0] for f in data_2]
                fp_2 = [f[1] for f in data_2]
                fig, axs = plt.subplots()
                ax = axs
                ax.set_xlabel(labels[0])
                ax.set_ylabel(labels[1])
                ax.scatter(fx_1, fx_2)
                ax.scatter(fp_1, fp_2)
            elif dim == 3:
                # 3D plot
                fx_1 = [f[0] for f in data_1]
                fx_2 = [f[1] for f in data_1]
                fx_3 = [f[2] for f in data_1]
                fp_1 = [f[0] for f in data_2]
                fp_2 = [f[1] for f in data_2]
                fp_3 = [f[2] for f in data_2]
                fig = plt.figure()
                ax  = plt.axes(projection='3d')
                ax.set_xlabel(labels[0])
                ax.set_ylabel(labels[1])
                ax.set_zlabel(labels[2])
                ax.scatter(fx_1, fx_2, fx_3)
                ax.scatter(fp_1, fp_2, fp_3)
        if fn is None:
            plt.show()
        else:
            plt.savefig(fn)
        if dump is not None:
            pickle.dump(fig, open(dump,'wb'))



    def plot_objective_space(self, fn=None, dump=None):
        """
        Plots a Pareto front together with other samples in
        objective-function space.

        The plot is saved in file <fn>. If <fn> is none, the plot is
        shown for interaction.

        The figure instance is saved to file <dump>, if not None.
        """
        dim = self.get_dim_f()
        if dim == 2:
            labels = ("F_1", "F_2")
        elif dim == 3:
            labels = ("F_1", "F_2", "F_3")
        self._plot_core(self.fx, self.fp, dim, labels, fn, dump)


    def plot_domain_space(self, fn=None, dump=None):
        """
        Same as plot_objective_space but in domain space.
        """
        dim = self.get_dim_x()
        if dim == 2:
            labels = ("X_1", "X_2")
        elif dim == 3:
            labels = ("X_1", "X_2", "X_3")
        self._plot_core(self.x, self.pf, dim, labels, fn, dump)


    def write(self, fn):
        """Saves the Pareto front to file <fn>."""
        import numpy as np
        xdim = self.get_dim_x()
        fdim = self.get_dim_f()
        paretolen = self.size()
        # prepare header
        header = ""
        for i in range(xdim):
            header += f"    X_{i+1}"
        for i in range(fdim):
            header += f"    F_{i+1}"
        all = np.zeros((paretolen, xdim + fdim))
        for i, (xx, fp) in enumerate(zip(self.pf, self.fp)):
            if xdim > 1:
                all[i,:len(xx)] = xx
                all[i,len(xx):] = fp
            else:
                all[i,:1] = xx
                all[i,1:] = fp
        np.savetxt(fn, all, header=header)


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        print("Usage: {} <project_directory> <output_prefix>".format(sys.argv[0]))
        sys.exit()

    paretoFront = ParetoFront.read_from_gmak(sys.argv[1])
    paretoFront.write(sys.argv[2] + "_pf.dat")
    paretoFront.plot_objective_space()
    paretoFront.plot_domain_space()
    paretoFront.plot_objective_space(fn=sys.argv[2] + "_pf.pdf",
                                     dump=sys.argv[2] + "_pf.fig")
    paretoFront.plot_domain_space(fn=sys.argv[2] + ".pdf",
                                  dump=sys.argv[2] + ".fig")

