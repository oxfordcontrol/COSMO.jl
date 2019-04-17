## Upcoming changes
- Update of PSD decomposition related code to Julia v1.0

## Version 0.4.1 (16. April 2019)
The most important changes are listed below:

- (`023f517`) Fix printing of rho
- (`cb2a483`) Memory efficient projections in BLAS
- (`c582f83`) Remove lower-triangle legacy code in MOI wrapper
- (`6b89271`) Dont export Model and convex sets
- (`d98c7ea`) Bug fix in norm computation
- (`0ad53d6`) Add Aox constraint type
- (`360ca0f`) Handle edge cases for `P`, `q`  in `assemble!`


## Version 0.4 (22. February 2019)
- The optional arguments `settings`, `x0`, `y0`, `s0` are now keyword arguments (breaks backward compatibility), see `d0f22aab`
- Some smaller bug fixes and function renaming
- Documentation via Documenter.jl


## Version 0.3.1 (2. February 2019)
- Support warm starting of primal and dual variables in MOI / JuMP interface
- Fix a bug when printing the header
- Improve unit test coverage


## Version 0.3 (7. January 2019)
__First public release of COSMO solver__
- Add MOI / JuMP support

