###############
Validation mode
###############

In many cases, one might be interested in estimating directly from
simulation the values of the target properties for hand-picked
parameter sets. In this scenario, one is interested not in the
parameter optimization *per se*, but instead in the evaluation,
without relying on surrogate models, of the target properties for
these different parameter sets.

The application described above can be achieved using the *validation
mode*.  In this mode, the :doc:`grid-shifting procedure
</overview/grid_shifting>` and the use of :doc:`surrogate models
</overview/surrogate_model>` are disabled.  The target properties are
evaluated only for the sampled grid points, and the expected values
and uncertainties are obtained directly from the simulations.

The validation mode is activated by using the option ``--validate`` in
the :doc:`command-line </usage/command-line>`.

.. tip:: To use the validation mode together with explicitly specified
   parameter sets, use the ``explicit`` :doc:`variation
   </usage/blocks/variation>`.

