#############
Grid Shifting
#############

After the values of the score function have been calculated, ``gmak``
allows one to update the origin of the grid based on those values,
effectively displacing the search-space boundaries in the parameter
space.  This is provided as a way to automatically explore nearby
parameter values whenever the results suggest that this is likely to
improve the optimal-parameter estimate. The main assumption is that,
within the limitations of the statistical uncertainties, one should be
able to roughly determine whether there is a significant chance that
the optimal parameter does not lie, in fact, within the lastly
considered search-space boundaries. The process will be referred to as
*grid shifting*.

The default grid-shifting strategy used in ``gmak`` is rather simple
and intuitive. It starts by ranking the grid cells from best to worst
score.  After that, the center-of-geometry (CG) of the
:math:`n_{\text{cut}} \times N_{1} \times N_{2} \times \cdots \times
N_{D}` best ones is calculated, where :math:`n_{\text{cut}}` is a
parameter defined in the :doc:`input file </usage/input_file>` and
:math:`N_{1}, N_{2}, \cdots ,N_{D}` are the dimensions of the
:ref:`main variation <overview/grid:main variation>`.  The positions
considered in this calculation are the :ref:`tuple indexes
<overview/grid:grid indexing>` of the points. If the location
:math:`\mathbf{x}_{\text{CG}}^{(k)} =
(x_{\text{CG},i}^{(k)})_{i=1\cdots D}` of the CG lies within some
pre-defined (in the :doc:`input file </usage/input_file>`) margins
:math:`(\delta_{i})_{i=1\cdots D}`, :math:`(\Delta_{i})_{i=1\cdots
D}`, *i.e.*

.. math::
   \delta_{i} N_{i} \leq x_{\text{CG},i}^{(k)} \leq \Delta_{i} N_{i} \quad , \quad i=1 \cdots D \, ,

or the number of grid shifts, :math:`k - 1`, reaches a maximum value
(also defined in the :doc:`input file </usage/input_file>`), then the
grid is *not* shifted. This signals the end of the run for ``gmak``.
Otherwise, the grid is shifted (see :ref:`overview/grid:variation
shifting` and :ref:`overview/grid:the parameter-search grid`) so that
the cell closest to the CG is the center of the new grid, and the
counter :math:`k` is incremented.

For flexibility, one can create customized grid-shifting procedures
and use them by means of the :ref:`customization API
<usage/customization_api:grid shifting>`.

.. note::
    The indexes defining the sampled grid points are not altered
    during grid shifting. This means that the parameter values are
    shifted together with the grid. 

.. note::
    Grid-shifting iterations are indepedent from each other.  The data
    and output files of each iteration are kept in disk, but are not
    accessed by the program during other iterations. This includes
    simulation files, property estimates, etc.
    
