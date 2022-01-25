
 .. list-table::
   :header-rows: 1
   :align: center

   * - Block name
     - Description
     - Remarks

   * - :doc:`variation </usage/blocks/variation>`
     - The ``$variation`` block sets up the variations (see :doc:`/overview/grid`).
     - 

   * - :doc:`gridshift </usage/blocks/gridshift>`
     - The ``$gridshift`` block sets up the routine for the grid-shifting procedure.
     - 

   * - :doc:`system </usage/blocks/system>`
     - The ``$system`` block sets up a routine for constructing the topology of a system.
     - 

   * - :doc:`coordinates </usage/blocks/coordinates>`
     - The ``$coordinates`` block sets up a routine for constructing initial configurations.
     - 

   * - :doc:`optimize </usage/blocks/optimize>`
     - The ``$optimize`` block sets up the composite properties that enter in the calculation of the score function, as well as their weights and uncertainty tolerances.
     - 

   * - :doc:`compute </usage/blocks/compute>`
     - The ``$compute`` block sets up a composite property and the associated component properties (see :ref:`overview/properties:component and composite properties`).
     - 

   * - :doc:`grid </usage/blocks/grid>`
     - The ``$grid`` block sets up the :ref:`sampled points <overview/grid:the list of sampled points>`.
     - Depends on :doc:`$variation <blocks/variation>`, :doc:`$gridshift <blocks/gridshift>`.

   * - :doc:`protocol </usage/blocks/protocol>`
     - The ``$protocol`` block sets up a :doc:`protocol </overview/protocols>`.
     - 
