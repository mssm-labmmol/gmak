# Requirements

- ```alchemlyb``` for MBAR-based free-energy calculations

# TODO

- More atoms in topology
- $molecules -> $systems; itp -> top;
- More general support for complex topologies
- What about dgsolv?
- Refactor read_input for god's sake


# Documentation Points

- dgsolv protocol:
    - uses MBAR so calc-lambda-neighbors -1
    - assumes NPT production runs in MBAR
    - assumes desolvation, so calculates -DeltaG_{A \to B}
