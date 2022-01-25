#####
Score
#####

The quality of each grid point is evaluated based on a score function,
which plays the role of the objective function to be minimized in the
parameter-optimization problem.

The default score function is a weighted root mean square difference
(RMSD) between the estimated property values and their corresponding
reference values,

.. math::
   s_{\text{RMSD}}(i, \mathbf{a}^{(k)}, \mathbf{w}, \mathbf{r}) = \sqrt{\frac{1}{\sum\limits_{p=1}^{N_P}w_p}\,\sum\limits_{p=1}^{N_P} w_{p} \left( a^{(k)}_{pi} - r_p\right)^2}
   \quad ,

where

-  the index :math:`i` is the :ref:`linear index <overview/grid:grid
   indexing>` of the grid point under evaluation;

-  the index :math:`p` runs over the composite properties;

-  :math:`N_P` is the number of composite properties;

-  :math:`\mathbf{w} = (w_1,\cdots,w_{N_P})` are the weights of the
   properties (given in the :doc:`input file </usage/input_file>`);

-  :math:`\mathbf{r} = (r_1,\cdots,r_{N_P})` are the reference values of
   the properties (given in the :doc:`input file </usage/input_file>`);

-  :math:`\mathbf{a}^{{(k)}} = (a_{pi}^{(k)})_{p}` is the vector of
   surrogate-model-estimated values of the properties for the grid
   point :math:`i`.

For flexibility, one can define customized score functions by means of
the :ref:`customization API <usage/customization_api:score function>`.
