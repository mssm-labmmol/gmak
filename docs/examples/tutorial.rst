========
Tutorial
========


This tutorial teaches the user how to:

- run a simple parameter-optimization job using ``gmak``

- analyze the results using the :doc:`post-processing module </usage/post_processing>`

- implement a custom property using the :doc:`customization API </usage/customization_api>`

- validate the results for some parameters using the :doc:`validation mode </usage/validation_mode>`

The system considered in this tutorial is the rigid three-point
water model OPC3
:footcite:`izadi16_accur_limit_rigid_water_model`. In particular,
the values of the :math:`C_6` and :math:`C_{12}` parameters of the oxygen atom
are optimized so as to make the model compatible with the GROMOS
treatment of long-range Lennard-Jonnes interactions (straight cutoff
truncation at a distance of 1.4 nm)
:footcite:`oostenbrink04_biomol_force_field_based_free`.  These two
parameters will be coupled by keeping their ratio constant, *i.e.*



.. math::

    \frac{C_{12}}{C_{6}} = 0.0010229709330817824 \text{ nm}^6

In this manner, the parameter search space effectively has only one
dimension.  Furthermore, this allows us to illustrate how to couple
the values of parameters within ``gmak``.

All the simulations and analyses will be carried out based on the
GROMACS software package, mostly relying on the GROMACS-compatible
:doc:`protocols </overview/protocols>` and :doc:`properties </overview/properties>`. The target properties of the
parametrization are the pure-liquid density and the enthalpy of
vaporization. In the validation stage, the self-diffusion
coefficient (based on the displacement of the oxygen atoms) is also
considered to illustrate the implementation of a :doc:`custom property </usage/customization_api>`.

.. note::

    Throughout this tutorial, ``ROOT`` denotes the root directory of the
    ``gmak`` `Git repository <http:github.com/mssm-labmmol/gmak>`_. All command snippets are executed from inside
    the directory ``$ROOT/docs/tutorial``, referred to as the tutorial
    directory. All relative paths are relative to the tutorial
    directory.

File structure
--------------

The tutorial directory is structured as follows:

::

    .
    ├── custom.py
    ├── input_files
    │   ├── gromos53a6.ff
    │   │   ├── aminoacids.c.tdb
    │   │   ├── aminoacids.hdb
    │   │   ├── aminoacids.n.tdb
    │   │   ├── aminoacids.r2b
    │   │   ├── aminoacids.rtp
    │   │   ├── aminoacids.vsd
    │   │   ├── atomtypes.atp
    │   │   ├── ffbonded.itp
    │   │   ├── ff_dum.itp
    │   │   ├── ffnonbonded.itp
    │   │   ├── forcefield.doc
    │   │   ├── forcefield.itp
    │   │   ├── ions.itp
    │   │   ├── spce.itp
    │   │   ├── spc.itp
    │   │   ├── tip3p.itp
    │   │   ├── tip4p.itp
    │   │   └── watermodels.dat
    │   ├── mdp
    │   │   ├── em_liq.mdp
    │   │   ├── md_liq.mdp
    │   │   ├── npt_liq.mdp
    │   │   └── nvt_liq.mdp
    │   ├── opc3_1024.top
    │   ├── opc3.gro
    │   ├── OW.ndx
    │   ├── tutorial_01.gmi
    │   └── tutorial_01_validation.gmi

These files are briefly explained below:

``custom.py``
    This is the file used to interact with the
    :doc:`customization API </usage/customization_api>`. This is where the calculation of the
    self-diffusion coefficient is implemented.

``gromos53a6.ff``
    This is the directory containing the GROMACS
    files for the GROMOS 53A6 force field. These files are only used
    to set the form of the Lennard-Jones potential to the
    :math:`C_6`-:math:`C_{12}` one (rather than the :math:`\epsilon`-:math:`\sigma`
    one), to set the mixing rule to the geometric-mean one (rather
    than the arithmetic-mean or the Lorentz-Berthelot one), and to
    define the atom types used in the simulations (in this case, H
    for the hydrogen atom and OW for the oxygen atom).

``mdp``
    This is the directory containing the input-parameter
    files (``.mdp`` files) for the GROMACS simulations.  These files
    set up a simulation protocol consisting of the following steps:

    1. Energy minimization (``em``)

    2. Equilibration in an NVT ensemble at 298.15K (``nvt``)

    3. Equilibration in an NPT ensemble at 298.15K and 1 bar (``npt``)

    4. Production run in the same conditions as the previous step (``md``)

    More details (*e.g.*, the length of the simulations) can be
    consulted in the files themselves.

``opc3_1024.top``
    This is the GROMACS topology file
    corresponding to a system of 1024 molecules of OPC3 water.  This
    is the topology file used for all simulations.

``opc3.gro``
    This is the geometry structure of a single OPC3
    water molecule in vacuum. This structure will be used by ``gmak``
    to construct the initial coordinates for all simulations.

``tutorial_01.gmi``
    This is the :doc:`input file </usage/input_file>` of the
    parameter-optimization job (described in details below).

``tutorial_01_validation.gmi``
    This is the :doc:`input file </usage/input_file>` of the
    validation job (described in more details below).

``OW.ndx``
    This is a GROMACS index file containing a single
    group named ``OW`` with the indexes of all oxygen atoms in the
    simulation box.

Parameter optimization
----------------------

The parameter-optimization part of this tutorial consists of three
parts. First, the content of the input file is explained, including
references to relevant parts of the documentation that relate to
each part of the file. After that, we show how to invoke ``gmak``
from the command line in order to carry out the
optimization. Finally, the results of the job are briefly analyzed
using the :doc:`post-processing module </usage/post_processing>`. In this analysis, the three
top-performing grid points in the estimated Pareto front are
selected for further validation.

Input file
~~~~~~~~~~

The input file of the parameter-optimization job is described part
by part below.

.. code:: gmi


    workdir tutorial_01




This sets the working directory to ``tutorial_01``.


.. seealso::

    :ref:`workdir <usage/input_file:global input parameters>`
        Global input-file parameter.


.. code:: gmi


    $variation
    name main
    pars V_OW
    type cartesian
    start 2.797579e-03
    step  0.02e-03
    size 33
    $end





This block configures the :ref:`main variation <overview/grid:main variation>`
to explore values from
2.797579e-03 to 3.437579e-03 using 33 points adjacently spaced by 0.02e-03.
The string ``V_OW`` indicates that these values replace the value of the
:math:`C_6` parameter of the ``OW`` atom type (``V`` is
interpreted as :math:`C_6` because the Lennard-Jones potential is
in the :math:`C_6`-:math:`C_{12}` form).


.. code:: gmi


    $variation
    name c12
    pars W_OW
    type coupled
    using main
    function x[0]*0.0010229709330817824
    $end





This block configures an additional variation named ``c12`` that is
:ref:`coupled to the main variation <overview/grid:coupled variation>` by a coupling function


.. math::

    f: \mathbb{R} \to \mathbb{R}




.. math::

    x_0 \mapsto 0.0010229709330817824 \cdot x_0


where :math:`x_0 \in \mathbb{R}` is an element of the main variation.
The string ``W_OW`` indicates that the elements of this variation
replace the value of the :math:`C_{12}` parameter of the ``OW`` atomtype
(``W`` is interpreted as :math:`C_{12}`
because the Lennard-Jones potential is
in the :math:`C_6`-:math:`C_{12}` form).

.. seealso::

    :doc:`$variation </usage/blocks/variation>`
        Input-file block.

    :ref:`overview/grid:variations`
        Section about variations in ``gmak``.

    :doc:`/overview/interaction_parameters`
        Section about interaction parameters in ``gmak``.


.. code:: gmi


    $gridshift
    maxshifts 10
    ncut 0.10
    margins 0.25 0.75
    $end





This block sets a default :doc:`grid-shifting procedure </overview/grid_shifting>` 
with a maximum number of iterations of 10, a :math:`n_\text{cut}` value of 0.10 and
margins :math:`\delta_1 = 0.25` and :math:`\Delta_1 = 0.75`.

.. seealso::

    :doc:`$gridshift </usage/blocks/gridshift>`
        Input-file block.

    :doc:`/overview/grid_shifting`
        Section about grid shifting.

.. code:: gmi


    $grid
    samples 0 16 32
    $end





This block sets the :ref:`sampled grid points <overview/grid:the list of sampled points>` to those with :ref:`linear indexes <overview/grid:grid indexing>`
0 (first point), 16 (middle point) and 32 (last point).

.. seealso::

    :doc:`$grid </usage/blocks/grid>`
        Input-file block.

    :doc:`/overview/grid`
        Section about the parameter-search grid in ``gmak``.

.. code:: gmi


    $coordinates
    name opc3_1024
    type gmx_liquid
    coords input_files/opc3.gro
    nmols 1024
    box cubic 3.135
    $end





This block sets a :ref:`configuration-construction routine <overview/coordinates:configuration-construction routines>` named ``opc3_1024`` to construct the initial configuration of the simulations.
It creates a :ref:`pure-liquid simulation box <overview/coordinates:pure-liquid configuration>` containing 1024 molecules of OPC3 water within a cubic box with edge length of 3.135 nm.
The basic molecular structure (one water molecule) replicated within the box is in the file ``input_files/opc3.gro``.

.. seealso::

    :doc:`$coordinates </usage/blocks/coordinates>`
        Input-file block.

    :doc:`/overview/coordinates`
        Section about coordinates in ``gmak``.

.. code:: gmi


    $system
    name opc3_1024
    type gmx
    template input_files/opc3_1024.top
    $end





This block sets a :ref:`GROMACS-compatible system <overview/systems_and_topologies:gromacs-compatible systems>`
named ``opc3_1024`` that uses the file ``input_files/opc3_1024.top`` as a template topology.

.. seealso::

    :doc:`$system </usage/blocks/system>`
        Input-file block.

    :doc:`/overview/systems_and_topologies`
        Section about systems and topologies in ``gmak``.

.. code:: gmi


    $protocol
    name opc3_1024
    type gmx
    system opc3_1024
    coords opc3_1024
    mdps input_files/mdp/em_liq.mdp input_files/mdp/nvt_liq.mdp input_files/mdp/npt_liq.mdp input_files/mdp/md_liq.mdp
    maxsteps 2500000
    $end





This block sets a :ref:`GROMACS-compatible general protocol <overview/protocols:gromacs-compatible general protocol>`
named ``opc3_1024`` that relies on the ``opc3_1024`` system (defined above)
and on the ``opc3_1024`` coordinates (defined above).
The input parameters of the simulation are given in the files
``input_files/mdp/em_liq.mdp`` to ``input_files/mdp/md_liq.mdp``.
The production run of the protocol is limited to a maximum duration of 2500000 steps (5 ns).

.. seealso::

    :doc:`$protocol </usage/blocks/protocol>`
        Input-file block.

    :doc:`/overview/protocols`
        Section about protocols in ``gmak``.

.. code:: gmi


    $compute
    name dens
    type density
    protocols opc3_1024
    surrogate_model linear
    $end





This block configures the program to compute the :ref:`density <overview/properties:density>` based on the production run of the ``opc3_1024`` protocol (defined above).
The surrogate model selected to compute the estimates of this property for the grid points that are not simulated is the :ref:`linear interpolation <overview/surrogate_model:linear/cubic interpolation>`.
The property is given the name ``dens``.

.. code:: gmi


    $compute
    name dhvap
    type dhvap
    protocols opc3_1024 none none
    surrogate_model linear
    C -7.186
    nmols from coordinates opc3_1024
    $end





This block configures the program to compute the :ref:`enthalpy of vaporization <overview/properties:enthalpy of vaporization>`
using the production run of the ``opc3_1024`` protocol (defined above) for the calculation of the liquid-phase potential energy.
The gas-phase potential energy and the polarization-energy correction are not calculated by the program based on simulations, which is indicated by associating them
with protocols named ``none``.
However, a constant corrrection of -7.186 kJ/mol is used, which encompasses the polarization-energy correction.
The number of molecules in the liquid phase is recycled from the ``opc3_1024`` coordinates (defined above).
The surrogate model selected to compute the estimates of this property for the grid points that are not simulated is the :ref:`linear interpolation <overview/surrogate_model:linear/cubic interpolation>`.
The property is given the name ``dhvap``.

.. seealso::

    :doc:`$compute </usage/blocks/compute>`
        Input-file block.

    :doc:`/overview/properties`
        Section about properties in ``gmak``.

    :doc:`/overview/surrogate_model`
        Section about surrogate models in ``gmak``.


.. code:: gmi


    $optimize
    properties   dens     dhvap
    references   997.00   43.989
    weights      1.0      1.0
    tolerances   0.30     0.10
    $end





This block sets a :doc:`default score function </overview/score>` based on the properties named ``dens`` and ``dhvap``.
The reference values of the properties are 997 kg/m\ :sup:`3`\ and 43.989 kJ/mol, respectively.
The weights of the properties are both 1.0.
This block also sets the tolerances for the statistical errors of these properties: 0.3 kg/m\ :sup:`3`\ and 0.1 kJ/mol, respectively.

.. seealso::

    :doc:`$optimize </usage/blocks/optimize>`
        Input-file block.

    :doc:`/overview/score`
        Section about the score function in ``gmak``.

    :ref:`Simulation Extensions <extensions1>`
        Section about the simulation extensions for GROMACS-compatible protocols in ``gmak``.

Running the job
~~~~~~~~~~~~~~~

Running the optimization job is very simple: in the command-line,
one can execute the ``gmak`` program as follows:

.. code:: bash

    gmak --gmx $GMXPATH --gnp $NPROCS input_files/tutorial_01.gmi

where ``$GMXPATH`` should be replaced by the path of the ``gmx``
binary and ``$NPROCS`` should be replaced by the desired number of
parallel threads (this number is passed along to the option ``-nt``
of ``gmx mdrun``). These two options are not mandatory---if they are
not supplied, the program will guess the path of the ``gmx`` binary
and delegate the choice of the number of threads to the ``mdrun``
program.

.. seealso::

    :doc:`/usage/command-line`
        Section about the ``gmak`` command.

Post-processing
~~~~~~~~~~~~~~~

After the job has completed, a new directory named ``tutorial_01``
(as specified in the input file) should have been created, storing
the :doc:`output files </usage/output_files>` of the job. Out of these files, only the binary
state file ``tutorial_01/state_%jobid.bin`` is used, where ``%jobid``
is the PID of the ``gmak`` job and is specific to your run. In our
case, this file is ``tutorial_01/state_8949.bin``, and will by
analyzed using the :doc:`post-processing module </usage/post_processing>`.

.. seealso::

    :doc:`/usage/output_files`
        Section about output files in ``gmak``.

    :doc:`/usage/post_processing`
        Section about the post-processing module in ``gmak``.

In a Python interpreter session, import the post-processing module
and read the state binary file:

.. code:: python

    import gmak.post_processing as pp

    jobdata = pp.GmakOutput.from_gmak_bin('%s/docs/tutorial/tutorial_01/state_8949.bin' % ROOT)

The variable ``jobdata`` is an instance of the
:py:class:`~gmak.post_processing.GmakOutput` class. It contains in
its attributes the main-variation elements, the estimates and
errors of the properties and the score for all grid points and all
grid-shift iterations. This data can be visualized more
effectively by converting this variable to a
:py:class:`pandas.DataFrame`, as shown below:

.. code:: python

    df = jobdata.get_dataframe()
    print(df)

::

                      (X, 1)  (dens, mu)  (dens, sigma)  (dens, diff)  (dhvap, mu)  (dhvap, sigma)  (dhvap, diff)  (score, mu)
    grid gridpoint                                                                                                            
    0    0          0.002798  992.678092       0.226104     -4.321908    44.365680        0.006418       0.376680     3.067636
         1          0.002818  991.915764       0.228639     -5.084236    44.259750        0.006425       0.270750     3.600191
         2          0.002838  991.153437       0.231175     -5.846563    44.153820        0.006432       0.164820     4.135787
         3          0.002858  990.391110       0.233710     -6.608890    44.047890        0.006439       0.058890     4.673377
         4          0.002878  989.628782       0.236245     -7.371218    43.941960        0.006446      -0.047040     5.212344
    ...                  ...         ...            ...           ...          ...             ...            ...          ...
    1    28         0.003058  982.513245       0.259369    -14.486755    43.003454        0.006766      -0.985546    10.267360
         29         0.003078  981.769364       0.260034    -15.230636    42.900889        0.006642      -1.088111    10.797135
         30         0.003098  981.025482       0.260699    -15.974518    42.798325        0.006518      -1.190675    11.327024
         31         0.003118  980.281601       0.261365    -16.718399    42.695761        0.006393      -1.293239    11.857009
         32         0.003138  979.537719       0.262030    -17.462281    42.593197        0.006269      -1.395803    12.387080

    [66 rows x 8 columns]

Another advantage of using the dataframe is the ease of interacting
with the underlying data. For example, the entries above can easily be
ordered by the value of the score function with the
:py:meth:`~pandas.DataFrame.sort_values` method:

.. code:: python

    df_sorted = df.sort_values(('score', 'mu'))
    print(df_sorted)

::

                      (X, 1)  (dens, mu)  (dens, sigma)  (dens, diff)  (dhvap, mu)  (dhvap, sigma)  (dhvap, diff)  (score, mu)
    grid gridpoint                                                                                                            
    1    11         0.002718  996.728604       0.264170     -0.271396    44.915398        0.008116       0.926398     0.682594
         10         0.002698  997.786360       0.266727      0.786360    45.051633        0.008087       1.062633     0.934759
         12         0.002738  995.670847       0.261613     -1.329153    44.779163        0.008144       0.790163     1.093390
         9          0.002678  998.844116       0.269284      1.844116    45.187867        0.008059       1.198867     1.555321
         13         0.002758  994.613091       0.259056     -2.386909    44.642929        0.008172       0.653929     1.749994
    ...                  ...         ...            ...           ...          ...             ...            ...          ...
    0    28         0.003358  973.295909       0.245010    -23.704091    41.730591        0.006343      -2.258409    16.837226
         29         0.003378  972.697163       0.243205    -24.302837    41.652240        0.006328      -2.336760    17.263955
         30         0.003398  972.098418       0.241400    -24.901582    41.573889        0.006313      -2.415111    17.690697
         31         0.003418  971.499673       0.239595    -25.500327    41.495539        0.006297      -2.493461    18.117451
         32         0.003438  970.900927       0.237791    -26.099073    41.417188        0.006282      -2.571812    18.544215

    [66 rows x 8 columns]

There is a lot of information to unpack from the dataframes above:

Indexing
    The index of the dataframe is a
    :py:class:`~pandas.MultiIndex` with the levels ``grid`` (the
    :doc:`grid-shift iteration </overview/grid_shifting>`) and ``gridpoint`` (the :ref:`linear index <overview/grid:grid indexing>`). This
    particular job involved only two grid-shift iterations.

Parameters
    The first :math:`d` columns ``(X,1)``, ``(X,2)``, ... ``(X,d)``
    contain the elements of the main variation (:math:`d` is the number of
    main-variation parameters). In this case, :math:`d=1`, and the column
    ``(X,1)`` corresponds to the :math:`C_6` coefficient of the OW atom
    type. The :math:`C_{12}` coeffient is not associated with main variation
    and is not shown.

Properties
    The next columns contain the estimated expected values
    (``mu``), statistical errors (``sigma``) and differences with respect to
    the reference value (``diff``) for all composite properties involved
    in the score function. In this case, the properties are only ``dens``
    and ``dhvap``.

Score
    The final columns show the estimated value (``mu``) and
    statistical error (``sigma``), when available, of the score function
    (``score``). For the default score function used in this tutorial, the
    error is not reported.

We proceed in the analysis by computing the approximate Pareto front
of the optimization problem:

.. code:: python

    pareto = jobdata.compute_pareto()

The method :py:meth:`~gmak.post_processing.GmakOutput.compute_pareto`
returns a list of the main-variation elements associated with the
Pareto front. In this case, the main variaton is one-dimensional, and
the variable ``pareto`` is a list of :math:`C_6` values:

::

    [0.0027175790000000003,
     0.0027375790000000004,
     0.0027575790000000005,
     0.002777579,
     0.002797579,
     0.002817579,
     0.0028375790000000002,
     0.0028575790000000003]


In order to verify the values of the properties and of the score
function for these points, it is first necessary to reindex the output
data based on the main variation. This can be done with
the :py:meth:`~gmak.post_processing.GmakOutput.groupby_X` method:

.. code:: python

    df_group = jobdata.groupby_X()
    print(df_group)

::

                     dens                           dhvap                          score
                       mu     sigma       diff         mu     sigma      diff         mu
    (X, 1)                                                                              
    0.002498  1008.363921  0.292297  11.363921  46.413979  0.007804  2.424979   8.216423
    0.002518  1007.306164  0.289740  10.306164  46.277744  0.007832  2.288744   7.465098
    0.002538  1006.248408  0.287183   9.248408  46.141509  0.007861  2.152509   6.714401
    0.002558  1005.190652  0.284626   8.190652  46.005275  0.007889  2.016275   5.964568
    0.002578  1004.132896  0.282069   7.132896  45.869040  0.007917  1.880040   5.215973
    0.002598  1003.075140  0.279512   6.075140  45.732806  0.007946  1.743806   4.469238
    0.002618  1002.017384  0.276955   5.017384  45.596571  0.007974  1.607571   3.725482
    0.002638  1000.959628  0.274398   3.959628  45.460336  0.008002  1.471336   2.986929
    0.002658   999.901872  0.271841   2.901872  45.324102  0.008031  1.335102   2.258690
    0.002678   998.844116  0.269284   1.844116  45.187867  0.008059  1.198867   1.555321
    0.002698   997.786360  0.266727   0.786360  45.051633  0.008087  1.062633   0.934759
    0.002718   996.728604  0.264170  -0.271396  44.915398  0.008116  0.926398   0.682594
    0.002738   995.670847  0.261613  -1.329153  44.779163  0.008144  0.790163   1.093390
    0.002758   994.613091  0.259056  -2.386909  44.642929  0.008172  0.653929   1.749994
    0.002778   993.555335  0.256499  -3.444665  44.506694  0.008201  0.517694   2.463100
    0.002798   992.587836  0.170007  -4.412164  44.368070  0.005218  0.379070   3.131367
    0.002818   991.677794  0.169905  -5.322206  44.246987  0.005231  0.257987   3.767826
    0.002838   990.924689  0.171005  -6.075311  44.142740  0.005184  0.153740   4.297285
    0.002858   990.171585  0.172109  -6.828415  44.038493  0.005138  0.049493   4.828552
    0.002878   989.418480  0.173215  -7.581520  43.934246  0.005092 -0.054754   5.361085
    0.002898   988.665376  0.174324  -8.334624  43.829999  0.005046 -0.159001   5.894542
    0.002918   987.912272  0.175436  -9.087728  43.725752  0.005001 -0.263248   6.428690
    0.002938   987.159167  0.176551  -9.840833  43.621505  0.004956 -0.367495   6.963371
    0.002958   986.406063  0.177668 -10.593937  43.517257  0.004911 -0.471743   7.498470
    0.002978   985.652958  0.178788 -11.347042  43.413010  0.004867 -0.575990   8.033903
    0.002998   984.899854  0.179910 -12.100146  43.308763  0.004823 -0.680237   8.569609
    0.003018   984.146749  0.181035 -12.853251  43.204516  0.004780 -0.784484   9.105538
    0.003038   983.393645  0.182163 -13.606355  43.100269  0.004737 -0.888731   9.641653
    0.003058   982.640540  0.183293 -14.359460  42.996022  0.004694 -0.992978  10.177926
    0.003078   981.887436  0.184425 -15.112564  42.891774  0.004652 -1.097226  10.714332
    0.003098   981.134332  0.185560 -15.865668  42.787527  0.004610 -1.201473  11.250852
    0.003118   980.381227  0.186697 -16.618773  42.683280  0.004569 -1.305720  11.787471
    0.003138   979.709914  0.186287 -17.290086  42.592823  0.004520 -1.396177  12.265737
    0.003158   979.283363  0.263057 -17.716637  42.514098  0.006498 -1.474902  12.570891
    0.003178   978.684617  0.261252 -18.315383  42.435747  0.006482 -1.553253  12.997420
    0.003198   978.085872  0.259447 -18.914128  42.357396  0.006467 -1.631604  13.423978
    0.003218   977.487126  0.257643 -19.512874  42.279046  0.006451 -1.709954  13.850563
    0.003238   976.888381  0.255838 -20.111619  42.200695  0.006436 -1.788305  14.277171
    0.003258   976.289636  0.254033 -20.710364  42.122344  0.006421 -1.866656  14.703802
    0.003278   975.690890  0.252228 -21.309110  42.043994  0.006405 -1.945006  15.130453
    0.003298   975.092145  0.250424 -21.907855  41.965643  0.006390 -2.023357  15.557122
    0.003318   974.493399  0.248619 -22.506601  41.887292  0.006374 -2.101708  15.983808
    0.003338   973.894654  0.246814 -23.105346  41.808942  0.006359 -2.180058  16.410510
    0.003358   973.295909  0.245010 -23.704091  41.730591  0.006343 -2.258409  16.837226
    0.003378   972.697163  0.243205 -24.302837  41.652240  0.006328 -2.336760  17.263955
    0.003398   972.098418  0.241400 -24.901582  41.573889  0.006313 -2.415111  17.690697
    0.003418   971.499673  0.239595 -25.500327  41.495539  0.006297 -2.493461  18.117451
    0.003438   970.900927  0.237791 -26.099073  41.417188  0.006282 -2.571812  18.544215

Finally, the rows that belong to the Pareto front can be filtered out from
the dataframe above using the :py:attr:`~pandas.DataFrame.loc` accessor:

.. code:: python
    :name: ParetoResults

    df_pareto = df_group.loc[pareto]
    print(df_pareto)

::

                    dens                          dhvap                         score
                      mu     sigma      diff         mu     sigma      diff        mu
    (X, 1)                                                                           
    0.002718  996.728604  0.264170 -0.271396  44.915398  0.008116  0.926398  0.682594
    0.002738  995.670847  0.261613 -1.329153  44.779163  0.008144  0.790163  1.093390
    0.002758  994.613091  0.259056 -2.386909  44.642929  0.008172  0.653929  1.749994
    0.002778  993.555335  0.256499 -3.444665  44.506694  0.008201  0.517694  2.463100
    0.002798  992.587836  0.170007 -4.412164  44.368070  0.005218  0.379070  3.131367
    0.002818  991.677794  0.169905 -5.322206  44.246987  0.005231  0.257987  3.767826
    0.002838  990.924689  0.171005 -6.075311  44.142740  0.005184  0.153740  4.297285
    0.002858  990.171585  0.172109 -6.828415  44.038493  0.005138  0.049493  4.828552

.. seealso::

    :py:meth:`~gmak.post_processing.GmakOutput.compute_pareto`, :py:meth:`~gmak.post_processing.GmakOutput.groupby_X`

Validation
----------

For the sake of illustration, the first three points in the Pareto
front above are chosen for further validation. For these points, the
values of the target properties are computed directly from simulation,
without relying on a surrogate model. Furthermore, the self-diffusion
coefficient (which is not implemented by default in the program) is
also calculated so as to illustrate the use of the customization API.

Input file
~~~~~~~~~~

The input file of the validation job
(``input_files/tutorial_01_validation.gmi``) is in many respects
very similar to the input file of the opimization job.  The most
important changes are the addition of a ``$compute`` block that sets
up the calculation of the self-diffusion coefficient and the
modification of some previously existing parts so as to take this
new property into account in the score function, explore only the
three parameter sets chosen for validation and change the directory
of the output files.

The first line

.. code:: gmi


    workdir tutorial_01_validation




sets the working directory to ``tutorial_01_validation``, so as to 
not overwrite the results of the optimization job.

The main variation is also set differently, and now reads

.. code:: gmi


    $variation
    name main
    pars V_OW
    type explicit
    dim 1
    values 2.71758e-03 2.73758e-03 2.75758e-03
    $end





This sets the main variation to explore the explicitly given values
of 2.71758e-03, 2.73758e-03 and 2.75758e-03. These values correspond
to the top-three performing points of the Pareto front selected above.

Note that only the main variation needs to be altered with respect
to the parameter-optimization job. The coupled variation (associated
with :math:`C_{12}`) automatically reflects the changes in the main
variation.

The parameter-search grid must also be updated to reflect
the new main variation. The block

.. code:: gmi


    $grid
    samples 0 1 2
    $end





sets the sampled grid points as those
with linear indexes 0, 1 and 2. This corresponds
to all three points of the parameter-search grid.

As explained above, the self-diffusion coefficient is also calculated
in the validation job. This requires setting up the following
``$compute`` block:

.. code:: gmi


    $compute
    name D
    type diffusion_coeff
    components msd
    protocols opc3_1024
    surrogate_model linear
    index_file input_files/OW.ndx
    group_name OW
    $end





This block sets a :ref:`composite property <overview/properties:composite properties>`
of type ``diffusion_coeff`` named ``D`` that has the ``msd`` property as its single :ref:`component property <overview/properties:component properties>`
and is calculated based on the ``opc3_1024`` protocol.
The estimates for the grid points that are not simulated are obtained based on a :ref:`linear-interpolation surrogate model <overview/surrogate_model:linear/cubic interpolation>`.
This property also requires an index file (given in ``index_file``) and a group (given in ``group_name``) to use
as a basis for calculating the self-diffusion coefficient.
The ``diffusion_coeff`` composite property and the ``msd`` component property
are defined in the customization file (see below).

Finally, the self-diffusion coefficient is also included
in the score function by altering the ``$optimize`` block
to

.. code:: gmi


    $optimize
    properties   dens     dhvap    D
    references   997.00   43.989   2.30
    weights      1.0      1.0      1.0
    tolerances   0.30     0.10     0.20
    $end

Customization
~~~~~~~~~~~~~

The file ``custom.py`` contains the implementation of the routine
for calculating the self-diffusion coefficient. We highly
recommend reading the :ref:`corresponding section <usage/customization_api:properties>` regarding custom
properties before continuing this tutorial in order to understand
the contents of the customization file. Also, consult the
:ref:`genindex` whenever a reference to an object needs to be
clarified.

The customization file is reproduced below with comments:

.. code:: python

    # import the customization API
    from gmak.api import *
    # gmak.config.ConfigVariables contains the path of the gmx binary
    from gmak.config import ConfigVariables
    # these are other useful modules
    import os
    import tempfile
    import re

    def msd_calc(topology, protocol_output, property_pars):
        """
        This function plays the role of the
        gmak.api_signatures.component_calculator() function
        (please consult the customization API documentation).

        It receives a topology file, the protocol-output dictionary with
        the simulation data and the input parameters (property_pars)
        passed to the program in the $compute block where the property is
        used.

        It returns a tuple with the expected value of the self-diffusion
        coefficient and the error in the estimate.
        """
        # path of the gmx binary
        gmx = ConfigVariables.gmx
        # tpr file of the simulation
        tpr = protocol_output['tpr']
        # xtc file of the simulation
        xtc = protocol_output['xtc']
        # create two temporary files to store the output of the ~gmx msd~ command
        d_out = tempfile.NamedTemporaryFile()
        msd_out = tempfile.NamedTemporaryFile(suffix='.xvg')
        # extract the path of the temporary files
        msd_out_path = msd_out.name
        d_out_path = d_out.name
        try:
            # try to retrieve the path of the index file and the name
            # of the group used as reference for calculating the
            # self-diffusion coefficient
            #
            # for the input file in this tutorial, this should work
            # since we supply the input parameters ~index_file~ and
            # ~group_name~ in the input file
            ndx = property_pars.index_file
            group = property_pars.group_name
            # finally, call the ~gmx msd~ command passing the group name
            # and the index file and storing the output in the temporary
            # files created above
            os.system(f"echo {group} | {gmx}  msd -n {ndx} -f {xtc} -s {tpr} -o {msd_out_path} > {d_out_path}")
        except AttributeError:
            # if retrieving the group name or the index file fails, use
            # all atoms as reference
            os.system(f"echo 0 | {gmx}  msd -f {xtc} -s {tpr} -o {msd_out_path} > {d_out_path}")
        # open the ~gmx msd~ output file search for the line containing
        # the value of the self-diffusion coefficient and the fitting
        # error
        with open(d_out_path, 'r') as fp:
            for line in fp:
                m = re.match('^D\[.*\]\s+(\S+)\s+\(\+/-\s+(\S+)\)', line)
                # the line that matches the regex above has the value of
                # the self-diffusion coefficient in the first group and
                # the fitting error in the second group
                if m:
                    value = float(m.group(1))
                    err = float(m.group(2))
        # close (thus removing) the temporary files
        d_out.close()
        msd_out.close()
        # return the value and the fitting error
        return (value, err)

    # see the documentation for gmak.api.add_custom_component_property()
    #
    # add ~msd~ as a custom component property associated with the 
    # ~msd_calc~ function defined above
    #
    # note that it does not correspond to a timeseries
    add_custom_component_property("msd",
                                  msd_calc,
                                  is_timeseries=False)

    # see the documentation for gmak.api.add_custom_composite_property()
    #
    # add ~diffusion_coeff~ as a custom composite property that is
    # recognized by the program
    #
    # since a calculator function is not supplied, it is implicitly
    # assumed that this property has only one component and that the value
    # and error correspond to those of the component property
    #
    # the association between ~diffusion_coeff~ and the ~msd~ component
    # property is established in the input file
    add_custom_composite_property("diffusion_coeff")

.. seealso::

    :doc:`/usage/customization_api`
        Section about the customization API in ``gmak``.

Running the job
~~~~~~~~~~~~~~~

The validation job is run very similarly to the
parameter-optimization job:

.. code:: bash

    gmak --gmx $GMXPATH --gnp $NPROCS --custom --validate input_files/tutorial_01_validation.gmi

where ``$GMXPATH`` is the path of your ``gmx`` binary and ``$NPROCS`` is
the number of parallel threads requested (option ``-nt`` of ``gmx mdrun``). The option ``--custom`` is used to read the file
``custom.py`` in the current directory. The option ``--validate``
is used to activate the validation mode.

.. seealso::

    :doc:`/usage/validation_mode`
        Section about the validation mode in ``gmak``.

Post-processing
~~~~~~~~~~~~~~~

After the job has completed, a new directory named
``tutorial_01_validation`` (as specified in the input file) should
have been created, storing the :doc:`output files </usage/output_files>` of the job. Out of
these files, only the binary state file
``tutorial_01_validation/state_%jobid.bin`` is used, where ``%jobid``
is the PID of the ``gmak`` job and is specific to your run. In our
case, this file is ``tutorial_01_validation/state_18767.bin``, and
will by analyzed using the :doc:`post-processing module </usage/post_processing>`.


In a Python interpreter session, import the post-processing module
and load the state binary file:

.. code:: python

    import gmak.post_processing as pp

    jobdata = pp.GmakOutput.from_gmak_bin('%s/docs/tutorial/tutorial_01_validation/state_18767.bin' % ROOT)

The variable ``jobdata`` is an instance of the
:py:class:`~gmak.post_processing.GmakOutput` class. It contains in
its attributes the main-variation elements, the estimates and
errors of the properties and the score for all grid points and all
grid-shift iterations. This data can be visualized more
effectively by converting this variable to a
:py:class:`pandas.DataFrame`, as shown below:

.. code:: python

    df = jobdata.get_dataframe()
    print(df)

::

                      (X, 1)  (dens, mu)  (dens, sigma)  (dens, diff)  (dhvap, mu)  (dhvap, sigma)  (dhvap, diff)  (D, mu)  (D, sigma)  (D, diff)  (score, mu)
    grid gridpoint                                                                                                                                            
    0    0          0.002718  996.446743       0.280857     -0.553257    44.841276        0.008420       0.852276   2.0700      0.0608    -0.2300     0.601489
         1          0.002738  995.550751       0.290102     -1.449249    44.719801        0.005507       0.730801   1.9924      0.0973    -0.3076     0.953766
         2          0.002758  994.159267       0.274727     -2.840733    44.610744        0.008299       0.621744   2.0848      0.1346    -0.2152     1.683512


From these results, one can see that there is a good match between
the surrogate-model based property estimates obtained previously
(the first three points below) and the estimates obtained directly
from simulation.

::

                    dens                          dhvap                         score
                      mu     sigma      diff         mu     sigma      diff        mu
    (X, 1)                                                                           
    0.002718  996.728604  0.264170 -0.271396  44.915398  0.008116  0.926398  0.682594
    0.002738  995.670847  0.261613 -1.329153  44.779163  0.008144  0.790163  1.093390
    0.002758  994.613091  0.259056 -2.386909  44.642929  0.008172  0.653929  1.749994
    0.002778  993.555335  0.256499 -3.444665  44.506694  0.008201  0.517694  2.463100
    0.002798  992.587836  0.170007 -4.412164  44.368070  0.005218  0.379070  3.131367
    0.002818  991.677794  0.169905 -5.322206  44.246987  0.005231  0.257987  3.767826
    0.002838  990.924689  0.171005 -6.075311  44.142740  0.005184  0.153740  4.297285
    0.002858  990.171585  0.172109 -6.828415  44.038493  0.005138  0.049493  4.828552

The first point, with :math:`C_6 \approx 0.002718` kJ mol\ :sup:`-1`\ nm\ :sup:`6`\, seems to be the
best choice among the three selected points.

References
----------

.. footbibliography::
