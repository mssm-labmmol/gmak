###############
Surrogate Model
###############

After the expected values and errors of the component properties have
been estimated for the sampled grid points, these quantities are
predicted for all points. The is done based on the *surrogate model*.
In ``gmak``, one can individually specify surrogate models for each
property. The options available are

-  cubic interpolation;
-  linear interpolation;
-  (homoscedastic) Gaussian Process regression.

For flexibility, one can also create customized surrogate models and
use them in the program by means of the :ref:`customization
API<usage/customization_api:surrogate model>`.

.. warning:: For the linear and cubic interpolation models, ``gmak``
   enforces the simulation of the *corners* of the grid---they are
   *automatically* added to the list of sampled points.

.. note:: The motivation for applying the surrogate model to predict
   the component properties (instead of the composite ones) is to
   allow for the future use of reweighting-based surrogate models.
   Since the component properties associated with a given composite
   property may originate from different simulation trajectories, it
   is not possible in general to use a reweighting-based surrogate
   model for composite properties.

Linear/Cubic Interpolation
==========================

The linear and cubic interpolations are carried out based on the
:py:func:`scipy.interpolate.griddata` function by specifying the
``method`` parameter as ``"linear"`` or ``"cubic"``. The arguments
``points`` and ``xi`` are set as the :ref:`tuple
indexes<overview/grid:grid indexing>` of the sampled grid points and
of the entire grid, respectively. The expected values and the
statistical errors are independently interpolated for each component
property.


Gaussian Process Regression (GPR)
=================================

The GPR surrogate model is based on the
:py:class:`sklearn.gaussian_process.GaussianProcessRegressor`.

The features of the model are the normalized axes of the :ref:`main
variation <overview/grid:main variation>`.  The normalization is done
independently for each axis so that the smallest interval containing
the normalized values of the corresponding parameter is :math:`[0,1]`.

Each regression round individually takes into account a component
property, using as target the expected values for the sampled grid
points.

The kernel is chosen from a list of commonly used ones

.. code-block:: python

   kernels = [
       WhiteKernel(noise_level=noiseLevel, noise_level_bounds='fixed') + ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * RBF(), 
       WhiteKernel(noise_level=noiseLevel, noise_level_bounds='fixed') + ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * Matern(),
       WhiteKernel(noise_level=noiseLevel, noise_level_bounds='fixed') + ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * DotProduct(),
       WhiteKernel(noise_level=noiseLevel, noise_level_bounds='fixed') + ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * ExpSineSquared(),
       WhiteKernel(noise_level=noiseLevel, noise_level_bounds='fixed') + ConstantKernel(constant_value_bounds=(1e-10,1e+10)) * RationalQuadratic()]


by using each one to indepedently fit a surrogate model and selecting
the one associated with the best log marginal likelihood.  The
training set of each fitting comprises all grid points, thus implying
an empty validation set. The ``noiseLevel`` is the mean statistical
error of the property along the sampled cells.

.. note:: Contrary to the linear and cubic interpolations, which
   require the expected values and statistical errors to be treated as
   two independent targets, the GPR model estimates these
   two quantities simultaneously.

