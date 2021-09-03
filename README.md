# Requirements

- ```alchemlyb``` for MBAR-based free-energy calculations

# TODO

- More atoms in topology
- $molecules -> $systems; itp -> top;
- More general support for complex topologies
- What about dgsolv?
- Refactor read_input for god's sake
- replace is_timeseries by apply_filter
- grep all TODO and todo


# Documentation Points

- dgsolv protocol:
    - uses MBAR so calc-lambda-neighbors -1
    - assumes NPT production runs in MBAR
    - assumes desolvation, so calculates -DeltaG_{A \to B}
