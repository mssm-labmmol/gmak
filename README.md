# TODO

## Improvements

- [DONE] python3 compatibility
- [TODO] allow specifying ```gmx``` path
- [TODO] ```make test``` for all tests
- [DONE] identify properties by id instead of nature, allowing using the same property more than once (/e.g./ optimize using density at two temperatures)
    - [DONE] gridOptimizer 
    - [DONE] rest of code...

## Tests

- [DONE] Surrogate models:
    - [DONE] MBAR 
    - [DONE] Interpolation 
        - [DONE] Nearest 
        - [DONE] Linear 
        - [DONE] Cubic 

- [DONE] GridOptimizer 

- [TODO] Full simple grid

# Possible Bugs/Things that don't make sense

- [DONE] Filtering seems to be filtering too much!
    - Actually, it is not - my test is actually at a low temperature, so auto-correlations are bigger in this case.
- [TODO] Using tolerance for Interpolation makes no sense!
