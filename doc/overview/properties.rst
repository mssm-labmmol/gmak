##########
Properties
##########

Component and Composite Properties
==================================

The *component* properties are those that are calculated directly from
the output of a single protocol.  One or more component properties can
be combined via some mathematical function to yield a *composite*
property.

In terms of implementation, the difference between these two types of
properties is that, for the former, the input of the calculation is
the :ref:`output of a protocol<overview/protocols:protocol output>`,
while, for the latter, it is the expected values and uncertainties of
the component properties.  Furthermore, component properties are
treated as intermediaries in the calculation of composite properties,
which in turn are treated as the final goal.  For this reason, the
estimated expected values and uncertainties of only composite
properties enter in the calculation of the :doc:`score function
</overview/score>` and are accessible via the :doc:`post-processing
module </usage/post_processing>`.

For flexibility, one can define custom component and composite
property types by means of the :ref:`customization
API<usage/customization_api:properties>`.

Component Properties
====================

In the context of ``gmak``, the calculated component properties are
determined by the requirements of the composite properties.  They are
calculated only for the sampled grid points and then estimated for all
grid points based on the :doc:`surrogate model
</overview/surrogate_model>`. Except when explicitly stated otherwise,
the expected value and uncertainty of a component property for a
sampled grid point are calculated following the ensemble-averaging
procedure :ref:`described below<overview/properties:ensemble averages
and uncertainties>`.

Component properties inherit the :py:class:`input
parameters<gmak.custom_attributes.CustomizableAttributesMixin.InputParameters>`
of their overlying composite property.  Different component-property
types require different input parameters to be computed, as well as
different data-access interfaces for the protocol-output object on
which the computation is based.

The following component-property types are provided by ``gmak`` and
are listed together with the corresponding requirements in terms of
input parameters and protocol-output interface. All properties rely
on the GROMACS software package.


Gmx-energy based Properties
---------------------------

Any property that is an energy component extractable from the GROMACS
energy file (``.edr`` file) by means of the ``gmx energy`` program can
be readily used as a component property. To do this, the
component-property type in the :doc:`input file</usage/input_file>`
must be specified as ``gmx_PROP``, where ``PROP`` is the name of the
property as represented in the ``gmx energy`` interactive prompt.

Input Parameters
    This component-property type requires no input parameters.


Protocol-output Interface
    This component property requires the object to be subscriptable
    and to expose the following items.

    'tpr'
        The absolute path of the ``.tpr`` file of the analyzed simulation.
    'edr'
        The absolute path of the ``.edr`` file of the analyzed simulation.


Polarization-energy Correction
------------------------------

The polarization-energy correction\ :footcite:`berendsen87_missin_term_effec_pair_poten` 
is available as a component property under the type ``polcorr``. It is
calculated as

.. math::
   E_{\text{pol}} = + \frac{1}{2} \frac{(\mu - \mu_{0})^{2}}{\alpha} \, ,

where

-  :math:`\mu` is the ensemble average of the simulated gas-phase molecular dipole moment;

-  :math:`\mu_{0}` is the experimental gas-phase molecular dipole moment;

-  :math:`\alpha` is the molecular polarizability.

The simulation samples of :math:`\mu` are extracted from the energy
file (``.edr`` file) using the ``gmx dipoles`` program within the
GROMACS package.

Input Parameters
    This component-property type requires the following input
    parameters.

    mu
        The experimental value :math:`\mu_0` of the molecular dipole
        moment in the gas phase (in :math:`\text{nm}^3`).

    alpha
        The value :math:`\alpha` of the molecular
        isotropic polarizability in the gas phase (in D).


Protocol-output Interface
    This component property requires the object to be subscriptable
    and to expose the following items.

    'tpr'
        The absolute path of the ``.tpr`` file of the analyzed gas-phase simulation.
    'edr'
        The absolute path of the ``.edr`` file of the analyzed gas-phase simulation.


.. _dg1:

Free-energy Difference
----------------------

The free-energy difference for an alchemical transformation is
available as a component property under the type ``dg``.

It is calculated using the MBAR technique\ :footcite:`shirts08_statis_optim_analy_sampl_from`,
as implemented in the `alchemlyb
<https://github.com/alchemistry/alchemlyb>`__ library.
This library provides easy-to-use methods for processing the energy
files of the simulations, preprocessing the data, and, finally,
estimating :math:`\Delta G` with MBAR.

Input Parameters
    This component-property type requires the following input
    parameters.

    temperature
        (optional) The reference temperature of the
        alchemical-transformation simulations. By default, it is
        inferred from the input-parameter files.


Protocol-output Interface
    This component property requires the object to be subscriptable
    and to expose the following items.

    'dhdl'
        A list of the absolute paths of the ``dhdl.xvg`` file of
        the production run for each alchemical state.


.. warning::
   This component property requires that the value of the option
   ``calc-lambda-neighbors`` is set to -1 in the input-parameter
   files.


Ensemble Averages and Uncertainties
-----------------------------------

In ``gmak``, ensemble averages and the corresponding statistical
errors are calculated as follows:

#. The values of the component property are obtained for each
   configuration of the ensemble and stored in an array;

#. The *ensemble average* is set as the mean of the array data;

#. The array is subsampled with a stride equal to its `autocorrelation
   statistical
   inefficiency <https://pymbar.readthedocs.io/en/master/timeseries.html#pymbar.timeseries.statisticalInefficiency>`__;

#. The *statistical error* is set as the standard error of the
   subsampled array, *i.e.* the standard deviation (using `Bessel's
   correction <https://en.wikipedia.org/wiki/Bessel%27s_correction>`__)
   divided by the square root of the number of uncorrelated samples.


Composite Properties
====================

The composite properties are calculated for all grid points based on
the values of the component properties estimated using the
:doc:`surrogate model </overview/surrogate_model>`.

Different composite-property types require different component
properties to be computed, and each of those requires specific input
parameters and protocol-output interfaces.  The input-parameter
requirements are satisfied by passing the input parameters of the
composite property to the underlying component properties. Also, each
component of the composite property is associated with an existing
protocol at the same time it is defined.

The following composite-property types are provided by ``gmak`` and
are listed together with the corresponding requirements in terms of
input parameters.

Density
-------

The expected value and uncertainty of the density (in kg m\
:superscript:`3`) are calculated as the expected value and uncertainty
of the ``gmx_Density`` component property (see
:ref:`overview/properties:gmx-energy based properties`).

It is identified with the type ``density`` in the :doc:`input file
</usage/input_file>`.

Enthalpy of Vaporization
------------------------

The expected value and uncertainty of the enthalpy of vaporization (in
kJ mol\ :superscript:`-1`) are calculated based on the
``gmx_Potential`` component property for a liquid-phase system (the
potential energy of the liquid), and, optionally, on the
``gmx_Potential`` (potential energy of the gas) and ``polcorr`` (the
polarization-energy correction) component properties for a gas-phase
system (see :ref:`overview/properties:gmx-energy based properties` and
:ref:`overview/properties:polarization-energy correction`).  They are
obtained as

.. math::
   \mu_{\Delta H_\text{vap}} & = \mu_{U_{\text{gas}}} - \frac{\mu_{U_{\text{liq}}}}{N_{\text{liq}}} - \mu_{E_{\text{pol}}} + RT +  C  \\
   \sigma_{\Delta H_\text{vap}} & = \sqrt{\sigma_{U_{\text{gas}}}^{2} + \frac{\sigma_{U_{\text{liq}}}^{2}}{N_{\text{liq}}^{2}} + \sigma_{E_{\text{pol}}}^{2} } \, ,

where

-  :math:`\mu_{U_{\text{gas}}}` (or :math:`\mu_{U_{\text{liq}}}`) is the
   expected value of the potential energy of the gas (or liquid)

-  :math:`\sigma_{U_{\text{gas}}}` (or
   :math:`\sigma_{U_{\text{liq}}}`) is the statistical error of the
   potential energy of the gas (or liquid)

-  :math:`\mu_{E_{\text{pol}}}` and :math:`\sigma_{E_{\text{pol}}}` are
   the expected value and statistical error, respectively, of the
   polarization-energy correction

-  :math:`N_{\text{liq}}` is the number of molecules in the liquid;

-  :math:`R` is the universal gas constant;

-  :math:`T` is the temperature;

-  :math:`C` contains any additional constant corrections (e.g. to
   account for quantum effects)


It is identified with the type ``dhvap`` in the :doc:`input file
</usage/input_file>`.

In practice, not all terms above are necessarily taken into
account—for instance, if there is no need to simulate the gas (*e.g.*,
for three-point rigid water molecules), one can force
:math:`\mu_{U_{\text{gas}}} = \sigma_{U_{\text{gas}}} = 0` by
configuring appropriate options in the :doc:`input file
</usage/input_file>`.  Likewise, polarization corrections can be
turned off, effectively setting :math:`\mu_{E_{\text{pol}}} =
\sigma_{E_{\text{pol}}} = 0`.

Input Paramaters
    This composite-property type requires the following input
    parameters.

    nmols
        The number :math:`N_\text{liq}` of molecules in the liquid
        phase.

    mu
        (required for polarization-energy correction only) The
        experimental value :math:`\mu_0` of the molecular dipole
        moment in the gas
        phase (in :math:`\text{nm}^3`).

    alpha
        (required for polarization-energy correction only) The value
        :math:`\alpha` of the molecular isotropic polarizability in
        the gas phase (in D).

    C
        (optional, default is 0.0) The value :math:`C` of the
        additional constant corrections.

    temperature
        (optional) The reference temperature :math:`T`. By default, it
        is inferred from the input-parameter of the liquid
        simulations.


Surface-tension Coefficient
---------------------------

The expected value and uncertainty of the surface-tension coefficient
(in mN m\ :superscript:`-1`) are calculated based on the expected
value and uncertainty of the ``gmx_#Surf*SurfTen`` component property
by halving these values and converting them to the appropriate units
(see :ref:`overview/properties:gmx-energy based properties`).

It is identified with the type ``gamma`` in the :doc:`input file
</usage/input_file>`. 

Free-energy Difference
----------------------

The expected value and uncertainty of the free-energy difference (in
kJ mol\ :superscript:`-1`) are calculated as the expected value and
uncertainty of the ``dg`` component property (see :ref:`Free-energy
Difference <dg1>`).

It is identified with the type ``dg`` in the :doc:`input file
</usage/input_file>`. 

References
==========

.. footbibliography::

