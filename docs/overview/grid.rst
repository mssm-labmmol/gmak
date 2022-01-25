#####################
Parameter-search Grid
#####################

A :math:`d`-dimensional *parameter-search grid* is a
:math:`N_1 \times N_2 \times \cdots \times N_d` multidimensional array
where each element, referred to as a *cell* or *point*, corresponds to
some combination of values of the force-field parameters in
optimization.

.. _Indexing:

Grid Indexing
=============

Any given cell of the grid can be unambiguously identified, or
indexed, in two ways:

Tuple Index
    using a *tuple index*
    :math:`(n_1, n_2, \cdots, n_d),\, n_\alpha \in [0, N_{\alpha} - 1]`,
    following zero-based array notation;

Linear Index
    using a *linear index* (or *flattened index*), whereby the
    :math:`i`-th member of the array is identified by the index
    :math:`i - 1`, following the `row-major
    order <https://en.wikipedia.org/wiki/Row-_and_column-major_order>`__.

.. admonition:: Example
   
   For example, in a :math:`10 \times 20` bidimensional grid, the
   element with tuple index :math:`(4,12)` has linear index
   :math:`4\times 20 + 12 = 92`, and the element with linear index
   :math:`136` has tuple index :math:`(6, 16)`, since :math:`136 = 6
   \times 20 + 16`. For grids of higher dimension, an analogous
   reasoning can be used to switch between these two indexing systems.

Variations
==========

The values of different force-field parameters may be coupled for many
different reasons. For example, a set of partial charges in a molecule
must usually satisfy electroneutrality. As another example, one may
impose that the changes on some parameters are coupled together by a
common scaling factor, to simplify the parameter-search space. The
common feature in all cases is that the parameters can be partitioned
into a group that varies *independently* and one that varies as a
*function* of the values of the first. These two groups will be
referred to as *variation domains*, with the first as a *main* and the
second as a *coupled* one. The array containing all the investigated
parameter values for a variation domain at a given moment of a
``gmak`` run is called a *variation*.


.. admonition:: Example

   For example, consider a neutral molecule with :math:`N` partial charges
   :math:`q_1, \cdots, q_N`. In the most general case, the charges can be
   varied systematically based on a *main* variation domain encompassing
   :math:`q_1,\cdots,q_{N-1}` and on a *coupled* variation domain with
   :math:`q_N = -(q_1+\cdots + q_{N-1})`.

Notation for Variations
-----------------------

More specifically, a :math:`d`-dimensional variation for
:math:`d'`-parameters is a :math:`d`-dimensional array where each
element is a :math:`d'`-dimensional tuple of real numbers, each the
value of some adjustable force-field parameter. The set of all
:math:`d`-dimensional variations with some given dimensions
:math:`\mathbf{D} \in \mathbb{N}^{d}` for :math:`d'` parameters is
represented as :math:`\mathcal{V}_d^{d'}(\mathbf{D})`.

.. admonition:: Example

   For instance, in the example of partial charges *above*, the main
   variation is a member of :math:`\mathcal{V}_{N-1}^{N-1}(\mathbf{D}_{q})`
   for some unspecified dimensions
   :math:`\mathbf{D}_{q} \in \mathbb{N}^{N-1}`. The coupled variation is a
   member of :math:`\mathcal{V}_{N-1}^{1}(\mathbf{D}_{q})`.

Function-induced Variations
---------------------------

Given some dimensions :math:`\mathbf{D}\in \mathbb{N}^{d}`, a function
:math:`\mathbf{h}: \mathbb{Z}^d \to \mathbb{R}^{d'}` naturally induces a
variation :math:`\mathbf{X} \in \mathcal{V}_d^{d'}(\mathbf{D})` such
that :math:`\mathbf{X}(\mathbf{n}) = \mathbf{h}(\mathbf{n})`
for all tuple indexes :math:`\mathbf{n}`. This means that, on its
finite domain, :math:`\mathbf{X}` agrees with :math:`\mathbf{h}`
everywhere. For clarity, the notation
:math:`\mathbf{h} \leadsto \mathbf{X}` is used to indicate that
:math:`\mathbf{X}` is constructed from :math:`\mathbf{h}` in this
fashion. In this case, :math:`\mathbf{X}` is called a *function-induced*
variation.

Variation Shifting
------------------

If a variation :math:`\mathbf{X} \in \mathcal{V}_d^{d'}(\mathbf{D})` is
induced by a function
:math:`\mathbf{h}:\mathbb{R}^{d} \to \mathbb{R}^{d'}` and
:math:`\mathbf{s}\in \mathbb{Z}^{d}` is a vector of integers, the
operation of *shifting the variation by* :math:`\mathbf{s}` results in a
new variation
:math:`\mathbf{X}|_{\mathbf{s}} \in \mathcal{V}_{d}^{d'}(\mathbf{D})`
with members

.. math:: \mathbf{X}|_{\mathbf{s}}(\mathbf{n}) = \mathbf{h}(\mathbf{n}+\mathbf{s})\, ,

where :math:`\mathbf{n}=(n_1, \cdots, n_d)` is a :ref:`tuple
index<overview/grid:grid indexing>`.  Alternatively,
:math:`\mathbf{X}|_{\mathbf{s}}` is the variation induced by the
function :math:`\mathbf{h}|_{\mathbf{s}}: \mathbb{R}^{d} \to
\mathbb{R}^{d'}` that results from shifting :math:`\mathbf{h}` by
:math:`\mathbf{s}`, with

.. math:: \mathbf{h}|_{\mathbf{s}}(\mathbf{n}) = \mathbf{h}(\mathbf{\mathbf{n}+\mathbf{s}}).

.. _CartesianVariation:

Cartesian Variation
-------------------

A cartesian variation is a special type of function-induced variation
:math:`\mathbf{h} \leadsto \mathbf{X} \in \mathcal{V}_{d}^{d}(\mathbf{D})`
for which each component of :math:`\mathbf{h}` is a discrete affine
function,

.. math::
   \mathbf{h}(n_{1},\cdots,n_{d}; \mathcal{O}, \Delta) & = (h_{1}(n_{1}), \cdots , h_{d}(n_d)) \\
   h_{i}(n_{i}; \mathcal{O}_i, \Delta_{i}) & = \mathcal{O}_{i} + \Delta_i n_i \, ,

where :math:`\mathcal{O}:=(\mathcal{O}_1, \cdots , \mathcal{O}_d)` is
the *origin* and :math:`\Delta:=(\Delta_1, \cdots , \Delta_d)` is the
*step*. This represents, perhaps, the simplest way of systematically and
independently varying a set of :math:`d` parameters—mapping them to a
cartesian grid with a given origin and step.

For flexibility, ``gmak`` provides also an *enhanced* cartesian
variation :math:`\mathbf{h}' \circ \mathbf{h} \leadsto \mathbf{X}`,
where
:math:`\circ` stands for function composition,
:math:`\mathbf{h}` is defined as above and
:math:`\mathbf{h}': \mathbb{R}^{d} \to \mathbb{R}^{d}` is a second,
user-specified, function (which will be referred to as a *scaling
function*) that is mapped over the members of the original cartesian
variation.

.. admonition:: Example

   For instance, in the example of partial charges *above*, the
   independent charges can be varied systematically and independently
   in the range :math:`\{0.0^{2}, 0.05^{2}, 0.10^{2}, \cdots,
   1.0^{2}\}` by choosing an enhanced cartesian variation with origin
   :math:`\mathcal{O} = \mathbf{0}\in \mathbb{R}^{N-1}`, step
   :math:`\Delta = \mathbf{0.05} \in \mathbb{R}^{N-1}` and a scaling
   function :math:`\mathbf{h}'(x_{1}, \cdots, x_{N-1}) = (x_1^2,
   \cdots, x_{N-1}^{2})`.


.. _MainVariation:

Main Variation
--------------

In the context of ``gmak``, the *main* variation
:math:`\mathbf{X}_\text{main} \in \mathcal{V}_{d}^{d}(\mathbf{D})` is
the variation holding the explored values of the :math:`d` force-field
parameters that are *independently* adjusted, at a given moment of
program execution. For completeness, it can also be represented as
:math:`\mathbf{X}_{\text{main}}^{(k)}`, where the superscript :math:`k`
indicates the number of grid-shifting iterations (see
:doc:`/overview/general_workflow`).

The (enhanced) cartesian variation is the only type of main variation
supported so far in ``gmak``.

.. _VariationsTie:

Coupled Variation
-----------------

In the context of ``gmak``, the *coupled* variation
:math:`\mathbf{X}_{\text{coupled}} \in \mathcal{V}_{d}^{d_{c}}(\mathbf{D})`
is the variation holding the explored values of the :math:`d_c`
force-field parameters that are *coupled* to the independent parameters.
The members are the results of mapping some function
:math:`\mathbf{f}: \mathbb{R}^{d} \to \mathbb{R}^{d_{c}}` (which will be
referred to as the *coupling function*) over the main variation,

.. math:: \mathbf{X}_{\mathrm{coupled}}(n_1, \cdots, n_d) = \mathbf{f}(\mathbf{X}_{\mathrm{main}}(n_1, \cdots, n_d))

where :math:`(n_1, \cdots , n_d)` is a :ref:`tuple
index<overview/grid:grid indexing>`. Note that :math:`\mathbf{h}
\leadsto \mathbf{X}_{\mathrm{main}} \Rightarrow \mathbf{f} \circ
\mathbf{h} \leadsto \mathbf{X}_{\mathrm{coupled}}`.  In other words,
if the main variation is function-induced, then the coupled variation
also is.

For completeness, the coupled variation can also be represented as
:math:`\mathbf{X}_{\text{coupled}}^{(k)}`, where the superscript
:math:`k` indicates the number of grid-shifting iterations (see
:doc:`/overview/general_workflow`).

.. admonition:: Example

   For instance, in the example of partial charges *above*, the
   coupled variation is defined by the function
   :math:`f:\mathbb{R}^{N-1} \to \mathbb{R}` such that :math:`f(x_{1},
   \cdots, x_{N-1}) = -(x_{1} + \cdots + x_{N-1})`.


The Parameter-search Grid
-------------------------

The parameter-search grid
:math:`\mathbf{X}_{\mathrm{pars}} \in \mathcal{V}_{d}^{d+d_{c}}(\mathbf{D})`
combines the main variation
:math:`\mathbf{X}_{\mathrm{main}} \in \mathcal{V}_{d}^d(\mathbf{D})` and
the coupled variation
:math:`\mathbf{X}_{\mathrm{coupled}}\in \mathcal{V}_{d}^{d_c}(\mathbf{D})`
into a single entity. In ``gmak``, this is what is referred to as the
"grid", when no further specifications are given. It is a variation with
members

.. math::
   \mathbf{X}_{\mathrm{pars}}(n_1, \cdots , n_d) = (\mathbf{X}_{\mathrm{main}}(n_1, \cdots , n_d),  \mathbf{X}_{\mathrm{coupled}}(n_1, \cdots , n_d)) \, ,

where :math:`(n_1, \cdots , n_d)` is a :ref:`tuple
index<overview/grid:grid indexing>`. One can easily prove that a
necessary and sufficient condition for the grid to be function-induced
is that the main variation is function-induced.  Since only (enhanced)
cartesian variations are allowed as main variations in ``gmak``, this
means that the parameter-search grid is automatically function-induced
and admits shifting.

For completeness, the grid can also be represented as
:math:`\mathbf{X}_{\mathrm{pars}}^{(k)}`, where the superscript
:math:`k` indicates the number of grid-shifting iterations (see
:doc:`/overview/general_workflow`).

The List of Sampled Points
==========================

The parameter-search grid is associated with an immutable list of
indexes that determines the grid points that are simulated. This list,
referred to as the *sampled* grid points, is initialized in the
:doc:`input file </usage/input_file>` and applies to all
:doc:`protocols </overview/protocols>` and :doc:`grid-shift iterations
</overview/grid_shifting>`.

.. note::
    This does not mean that the simulated force-field-parameter values
    are also immutable. In fact, these values change everytime the
    parameter-search grid is shifted. The corresponding grid-point
    indexes, however, are always the same.
