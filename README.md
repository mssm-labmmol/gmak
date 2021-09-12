# Requirements

- ```alchemlyb``` for MBAR-based free-energy calculations

# TODO

- ```setuptools``` installer with executable ```gmak```
- update example input file in share
- More general support for complex topologies (for instance, ALTER some atomtypes instead of defining all of them)
- Refactor read_input for god's sake
- grep all TODO and todo
- consistent naming of methods (e.g., CamelCase vs camel_case)

# Documentation Points

- dgsolv protocol:
    - uses MBAR so calc-lambda-neighbors -1
    - assumes NPT production runs in MBAR
    - assumes desolvation, so calculates -DeltaG_{A \to B}
