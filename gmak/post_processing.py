import pandas as pd
import re
import os
import numpy as np
from glob import glob
from collections import OrderedDict
from abc import ABC, abstractmethod


def default_agg(x):
    if len(x) > 1:
        return list(x)
    else:
        return x.iloc[0]


class GroupbyStrategy(ABC):

    @abstractmethod
    def groupby_X(self, gmak_output, agg=True, agg_arg=None, mu_agg=None, sigma_agg=None):
        pass


class GroupbyParameterValues(GroupbyStrategy):

    def __init__(self):
        pass

    def _groupby(self, gmak_output):
        # build full dataframe, with X and data
        full = pd.merge(gmak_output.X, gmak_output.data, how='inner',
                        left_index=True, right_index=True)
        Xcols = list(gmak_output.X.columns.values)
        return full.groupby(by=Xcols)

    def groupby_X(self, gmak_output, agg=True, agg_arg=None, mu_agg=None, sigma_agg=None):
        gb_df = self._groupby(gmak_output)
        propcols = gmak_output.data.columns.values
        parcols = gmak_output.X.columns.values
        # aggregate
        if not(agg):
            # collect the data for duplicate indexes in a list
            agg_func = OrderedDict()
            for proptuple in propcols:
                agg_func[proptuple] = default_agg
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


class GroupbyAbsoluteIndex(GroupbyStrategy):

    def __init__(self, aidx):
        """
        aidx is a DataFrame with MultiIndex([grid, gridpoint]) and column 'aidx'.
        """
        self.aidx = aidx

    def _groupby(self, gmak_output):
        # merge X with aidx
        full_X = pd.merge(gmak_output.X, self.aidx, how='inner',
                        left_index=True, right_index=True)
        # merge again
        full = pd.merge(full_X, gmak_output.data, how='inner',
                        left_index=True, right_index=True)
        # grouping column
        Xcols = list(self.aidx.columns.values)
        return full.groupby(by=Xcols)

    def groupby_X(self, gmak_output, agg=True, agg_arg=None, mu_agg=None, sigma_agg=None):
        gb_df = self._groupby(gmak_output)
        propcols = gmak_output.data.columns.values
        parcols = gmak_output.X.columns.values
        # aggregate
        if not(agg):
            # collect the data for duplicate indexes in a list
            agg_func = OrderedDict()
            for partuple in parcols:
                # partuple *must* aggregate, otherwise the parameters cannot be used 
                # as index
                agg_func[partuple] = "mean"
            for proptuple in propcols:
                agg_func[proptuple] = default_agg
        else:
            # build aggregate function
            if (agg_arg is not None):
                agg_func = agg_arg
            else:
                # build my own functions and replace by those given
                # by the user if it is the case
                agg_func = OrderedDict()
                for partuple in parcols:
                    agg_func[partuple] = "mean"
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
        agg_df = gb_df.agg(agg_func)
        # set X as index
        return agg_df.set_index(list(parcols))


class GmakOutput:
    """
    Class that collects the output of ``gmak`` jobs.

    Attributes
    ----------
    :var data: The output data of a ``gmak`` job. This contains the estimates
        (second column-level ``mu``), uncertainties (second column-level ``sigma``)
        and difference with respect to the reference value (second column-level
        ``diff``), for all explored grid points, of the composite properties
        that make up the score function, as well as the score itself. The data
        frame has two index levels: the grid-shift iteration and the linear
        index of the grid point.
    :vartype data: pandas.DataFrame
    :var X: The main-variation parameters. The index is the same as that of
        :py:data:`~gmak.post_processing.GmakOutput.data`.
    :vartype X: pandas.DataFrame
    """
    def __init__(self, df_X, df_data, groupby_strategy):
        """Initialize GmakOutput with data.
        Should not be used; use `from_gmak_*` instead.

        Parameters
        ----------
        df_X : pandas.DataFrame object
        df_data: pandas.DataFrame object groupby_strategy: GroupbyStrategy object
        """
        self.X = df_X
        self.data = df_data
        self.groupby_strategy = groupby_strategy

    @classmethod
    def from_gmak_bin(cls, binpath):
        """
        Creates a :py:class:`~gmak.post_processing.GmakOutput` instance from a
        binary state file.

        :param binpath: The path of the binary file
        :type binpath: str
        :return: A :py:class:`~gmak.post_processing.GmakOutput` instance containing the data in the binary file.
        :rtype: ~gmak.post_processing.GmakOutput
        """
        from gmak.state import State as State
        state = State()
        fp = open(binpath, "rb")
        state.setFromBinary(fp)
        fp.close()
        rec = state.record_book
        idx = rec.to_dataframe(tag="idx")
        aidx = rec.to_dataframe(tag="aidx")
        properties = rec.to_dataframe(tag="property")
        score = rec.to_dataframe(tag="score")
        X = rec.to_dataframe(tag="X")
        # create object
        df_X = pd.merge(idx, X, left_index=True, right_index=True)
        data = pd.merge(idx, properties, left_index=True, right_index=True)
        data = pd.merge(data, score, left_index=True, right_index=True)
        df_aidx = pd.merge(idx, aidx, left_index=True, right_index=True)
        df_X = df_X.set_index(["grid", "gridpoint"])
        data = data.set_index(["grid", "gridpoint"])
        df_aidx = df_aidx.set_index(["grid", "gridpoint"])
        groupby_strategy = GroupbyAbsoluteIndex(df_aidx)
        return cls(df_X, data, groupby_strategy)

    @classmethod
    def _from_gmak_job_noaidx(cls, jobdirpath, decimals=13):
        """Create a GmakOutput from the directory of a Gmak job.

        Parameters
        ----------
        jobdirpath : os.PathLike object or string

        decimals: int, default: 13
            Number of decimals to which the values of the parameters are
            rounded. This is *very* important to set appropriately, because
            precision issues may affect the grouping of the data by parameter
            values.

            Negative values -> don't round

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
            if decimals < 0:
                parameters += list(np.loadtxt(parametersfile))
            else:
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
        return cls(df_X, df_props_score, GroupbyParameterValues())

    def groupby_X(self, agg=True, agg_arg=None, mu_agg=None, sigma_agg=None):
        """Collects data that corresponds to the same force-field parameters and
        optionally aggregate these values using functions.

        :param agg: Indicates whether aggregation is desired. This has priority
            over ``agg_arg``, ``mu_agg`` and ``sigma_agg`` (default is
            :py:obj:`True`).
        :type agg: bool
        :param agg_arg: The `func` argument of the
            :py:meth:`~pandas.core.groupby.DataFrameGroupBy.aggregate` function
            used to aggregate values of a common group. This has priority over
            the parameters ``mu_agg`` and ``sigma_agg``.
        :type agg_arg: function, str, list, dict or None
        :param mu_agg: The function used to aggregate expected values and their
            differences with respect to reference data. By default, uses mean.
        :type mu_agg: function, str or None
        :param sigma_agg: The function used to aggregate uncertainties. By
            default, uses :math:`\sqrt{k_1^2 + \cdots + k_n^2}/n`, where the
            :math:`k_i` are the values to be aggregated.
        :type sigma_agg: function, str or None

        :return: A data frame with the force-field parameters as index and the
            score and properties as columns.  If aggregation is not requested, each
            entry in a column is either a float or a list of floats, depending on
            whether it was estimated for the force-field parameter once or more
            than once.  Otherwise, each entry results from aggregating the values
            for the same force-field parameter based on the parameters `agg_arg`,
            `mu_agg` and `sigma_agg`.
        :rtype: pandas.DataFrame
        """
        return self.groupby_strategy.groupby_X(self, agg, agg_arg, mu_agg,
                                               sigma_agg)


    def get_dataframe(self):
        """
        Returns a data frame containing both the main parameters and the
        properties and scores.

        :return: A data frame indexed by grid-shift iteration and gridpoint
            (linear index). The columns are the main parameters and the scores
            and property estimates.
        :rtype: pandas.DataFrame
        """
        return pd.merge(self.X, self.data, how='inner', left_index=True,
                        right_index=True)


    def compute_pareto(self, properties=None, groupby_X=None, **kwargs):
        """
        Returns the parameters that are Pareto optimal. It first groups the
        data using the ``groupby_X`` method, unless the parameter `groupby_X`
        is passed, in which case it is used instead.  The optimization problem
        considered is minimizing the absolute value of the difference between
        estimated properties and their reference values.

        :type properties: list of str
        :param properties: Properties considered in the optimization problem.
            By default, all properties are used.
        :type groupby_X: pandas.DataFrame
        :param groupby_X: data frame obtained from a previous call of
            the ``groupby_X`` method. If not given, the ``groupby_X`` method
            will be called with the keyword arguments ``kwargs``.
        :param kwargs:
            Additional keyword arguments passed to the aggregation method
            (``groupby_X``) called before computing the Pareto front.
        :return: The list of force-field parameters (each one a tuple of float)
            that are Pareto optimal.
        :rtype: list
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

