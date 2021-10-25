import pandas as pd
import re
import os
import numpy as np
from glob import glob
from collections import OrderedDict

class GmakOutput:
    """
    Class to collect and manipulate the output of Gmak jobs.

    Attributes
    ----------
    data : pandas.DataFrame object
        The output data of a Gmak job, in memory. The index is a
        pandas.MultiIndex object with levels ('grid', 'gridpoint') and
        columns:

        * (name_1, 'mu'), ... , (name_N, 'mu'): Expected value of the property
          name_i (float).
        * (name_1, 'sigma'), ... , (name_N, 'sigma'): Statistical uncertainty
          value of the property name_i (float).
        * (name_1, 'diff'), ... , (name_N, 'diff'): Difference with
          respect to reference value of the property name_i (float).
        * (score, 'mu'): The expected value of the score (float).

    X : pandas.DataFrame object
        The index is a pandas.MultiIndex object with levels ('grid',
        'gridpoint') and columns:

        * ('X', '1'), ...: The main variation parameters (float).

    """
    def __init__(self, df_X, df_data):
        """Initialize GmakOutput with data. Should not be used; use
        `from_gmak_job` instead.

        Parameters
        ----------
        df_X : pandas.DataFrame object
        df_data: pandas.DataFrame object

        """
        self.X = df_X
        self.data = df_data

    @classmethod
    def from_gmak_bin(cls, binpath):
        raise NotImplementedError

    @classmethod
    def from_gmak_job(cls, jobdirpath, decimals=13):
        """Create a GmakOutput from the directory of a Gmak job.

        Parameters
        ----------
        jobdirpath : os.PathLike object or string

        decimals: int, default: 13
            Number of decimals to which the values of the parameters are
            rounded. This is *very* important to set appropriately, because
            precision issues may affect the grouping of the data by parameter
            values.

        Returns
        -------
        GmakOutput
            A GmakOutput object with the data of the Gmak job in the
            directory.

        """
        # Convert to absolute path.
        a_jobdirpath = os.path.abspath(jobdirpath)

        # Check if input path is a directory.
        if not os.path.isdir(a_jobdirpath):
            raise IOError(f"{jobdirpath} is not a directory.")

        # Grid paths.
        grids = glob(os.path.join(a_jobdirpath, "grid_*"))
        grids.sort()

        # These will hold the columns.
        grid_numbers  = []
        point_numbers = []
        parameters    = []
        scores        = []
        # A bit more involved--dictionary where keys are tuples with
        # two levels (property name, mu or sigma) and values are lists
        # of property estimates or uncertainties.
        properties   = {}

        # Loop over grid_* inside directory.
        for grid_number, grid in enumerate(grids):
            # Catch, sort and loop over steps.
            steps = glob(os.path.join(grid, "step_*"))
            steps.sort()

            # Use only last step!
            step = steps[-1]

            scorefile = os.path.join(step, "full_data.dat")
            parametersfile = os.path.join(grid, "parameters_main.dat")

            # Infer property names and files.
            propfiles_EA_k = glob(os.path.join(step, "*_EA_k.dat"))
            propfiles_dEA_k = glob(os.path.join(step, "*_dEA_k.dat"))
            propfiles_diff = glob(os.path.join(step, "*_diff.dat"))
            propfiles_EA_k.sort()
            propfiles_dEA_k.sort()
            propfiles_diff.sort()
            propnames = [re.search('(.*)_EA_k.dat', os.path.basename(fn)).group(1)
                         for fn in propfiles_EA_k]
            # Init dict of properties if necessary.
            if (grid_number == 0):
                for prop in propnames:
                    properties[(prop, 'mu')] = []
                    properties[(prop, 'sigma')] = []
                    properties[(prop, 'diff')] = []
            # Fill data.
            parameters += list(np.around(np.loadtxt(parametersfile),
                                         decimals=decimals))
            new_scores = list(np.loadtxt(scorefile, comments=['#',])[:,-1])
            scores += new_scores
            number_of_points = len(new_scores)
            for i in range(number_of_points):
                grid_numbers.append(grid_number)
                point_numbers.append(i)
            # parameters?
            for prop, propmu, propsigma, propdiff in zip(
                    propnames, propfiles_EA_k, propfiles_dEA_k, propfiles_diff):
                # read property values as a np.ndarray then convert to list
                muvalues = np.loadtxt(propmu)
                sigmavalues = np.loadtxt(propsigma)
                diffvalues = np.loadtxt(propdiff)
                properties[(prop, 'mu')] += list(muvalues)
                properties[(prop, 'sigma')] += list(sigmavalues)
                properties[(prop, 'diff')] += list(diffvalues)

        df_props_score = pd.DataFrame.from_dict({
            'grid': grid_numbers,
            'gridpoint': point_numbers,
            ('score', 'mu'): scores,
            **properties,
        })

        # Process parameters into a dictionary where keys are the components
        # of each X, ('X', '1'), ('X', '2'), etc...
        # Unless it is 1D, in which case key is just 'X'
        parameter_dict = dict()
        try:
            ncomps = len(parameters[0])
            for c in range(ncomps):
                parameter_dict[('X', str(c+1))] = [x[c] for x in parameters]

        except:
            ncomps = 1
            parameter_dict['X'] = parameters

        df_X = pd.DataFrame.from_dict({
            'grid': grid_numbers,
            'gridpoint': point_numbers,
            **parameter_dict,
        })
        df_props_score.set_index(['grid', 'gridpoint'], inplace=True)
        df_props_score.columns = pd.MultiIndex.from_tuples(df_props_score.columns)
        df_X.set_index(['grid', 'gridpoint'], inplace=True)
        try:
            df_X.columns = pd.MultiIndex.from_tuples(df_X.columns)
        except:
            pass
        return cls(df_X, df_props_score)

    def groupby_X(self, agg=True, agg_arg=None, mu_agg=None, sigma_agg=None):
        """Collects data that corresponds to the same force-field parameters and
        optionally aggregate these values using functions.

        Parameters
        ----------
        agg : {True, False}, default: True
            Indicates whether aggregation is desired. This has priority over
            agg_arg, mu_agg and sigma_agg.

        agg_arg : function, str, list, dict or None, default: None
            If not None, it is the `func` argument of the
             `pandas.DataFrameGroupBy.aggregate` function used to aggregate.
            This has priority over the parameters mu_agg and sigma_agg.

        mu_agg : function, str or None, default: None
            If not None, use it to aggregate expected values and their
            differences with respect to reference data. By default, use mean.

        sigma_agg : function, str or None, default: None
            If not None, use it to aggregate uncertainties. By default, use 1/n
            * sqrt(k_1**2 + ... + k_n**2), where the k_i's are the n values to
            be aggregated.

        Returns
        -------
        pandas.DataFrame
            A DataFrame object with the force-field parameters as MultiIndex and
            the score and properties as columns (like `self.data`).

            If aggregation is not requested, each entry in a column is either a
            float or a list of floats, depending on whether it was estimated for
            the force-field parameter once or more than once.

            Otherwise, each entry results from aggregating the values for the
            same force-field parameter based on the parameters `agg_arg`,
            `mu_agg` and `sigma_agg`.
        """
        # build full dataframe, with X and data
        full = pd.merge(self.X, self.data, how='inner',
                        left_index=True, right_index=True)
        Xcols = list(self.X.columns.values)
        propcols = self.data.columns.values
        # groupby X
        gb_df = full.groupby(by=Xcols)
        # aggregate
        if not(agg):
            # collect the data for duplicate indexes in a list
            agg_func = lambda x: list(x) if len(x) > 1 else x.iloc[0]
        else:
            # build aggregate function
            if (agg_arg is not None):
                agg_func = agg_arg
            else:
                # build my own functions and replace by those given
                # by the user if it is the case
                agg_func = OrderedDict()
                for proptuple in propcols:
                    if (proptuple[1] == 'mu') or (proptuple[1] == 'diff'):
                        if mu_agg is None:
                            agg_func[proptuple] = "mean"
                        else:
                            agg_func[proptuple] = mu_agg
                    elif (proptuple[1] == 'sigma'):
                        if sigma_agg is None:
                            lfunc = lambda r: np.sqrt(np.sum(r*r))/len(r)
                            agg_func[proptuple] = lfunc
                        else:
                            agg_func[proptuple] = sigma_agg
        return gb_df.agg(agg_func)


    def get_dataframe(self):
        """
        Returns a pandas.DataFrame object containing both the main parameters
        and the properties and scores.

        Returns
        -------
        pandas.DataFrame
            A DataFrame object indexed by grid-shift iteration and gridpoint
            (linear index). The columns are the main parameters and the scores
            and property estimates. Naming convention follows that for `self.X`
            and `self.data`.
        """
        return pd.merge(self.X, self.data, how='inner', left_index=True,
                        right_index=True)


    def compute_pareto(self, properties=None, groupby_X=None, **kwargs):
        """
        Returns the indexes (parameters) that are Pareto optimal. It first
        groups the data by the parameters using `self.groupby_X`, unless the
        groupby_X DataFrame is passed, in which case it is used instead.  The
        optimization problem considered is minimizing the absolute value of the
        difference between estimated properties and their reference values.

        Parameters
        ----------
        properties : list of str, default: None
            Properties considered in the optimization problem. If None, all
            properties are used.

        groupby_X : pandas.DataFrame object, default: None
            DataFrame obtained from a previous call of `self.groupby_X`.
            If None, `self.groupby_X` will be called with the keyword
            arguments kwargs.

        **kwargs
            Additional keyword arguments passed to the aggregation method
            (`self.groupby_X`) called before computing the Pareto front.

        Returns
        -------
        List of force-field parameters that are Pareto optimal.
        """
        from pandas_pareto.pareto import compute_pareto as pd_pareto
        # get diff columns
        if properties is None:
            diffcols = [c for c in self.data.columns if c[1] == 'diff']
        else:
            diffcols = [c for c in self.data.columns if c[1] == 'diff' and c[0] in properties]
        # compute pareto with absolute values of diff columns
        if groupby_X is None:
            input = self.groupby_X(**kwargs)
        else:
            input = groupby_X
        return pd_pareto(input[diffcols].apply(np.abs))

