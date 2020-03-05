## Version 0.7.0 (05. March 2020)
- This release adds support for indirect linear system solver from the `IterativeSolvers.jl` package. We also drop support of PsdConeSquare constraints ins MOI. By default MOI will bridge them to PSDTriangle constraints. We further fixed a number of bugs related to the compact transformation for decomposable SDPs. Another change is that when used with `verbose_timing = true`, COSMO will now measure the initial factorisation time and the time for refactorisations separately.
The major commits are listed here:

- `439f218` Return NearlyFeasible status code if max_iter reached
- `fda29e1` Fix problem with empty columns in compact_transform
- `4324560` Drop support for MOI.PSDConeSquare
- `f0a99fe` ⏲  Add more detailed timing to solver
- `92aa32c` Store rho-updates in result variable
- `eb0be3b` Assemble dual variable like in Grone's Theorem
- `14e862e` Make Cholmod (SuiteSparse) default KKT solver
- `e64cfc9` Fix bug with CGIndirectKKTSolver and MINRESIndirectKKTSolver
- `f514caa` Merge pull request #106 from oxfordcontrol/pg/pardiso_fix
- `6290cbe` fixes minres
- `b2a0e49` fix pardiso default options
- `5ed8a9c` [ci skip] Add info about custom linear solvers
- `cfa6f76` Merge pull request #92 from oxfordcontrol/cg


## Version 0.6.0 (26. November 2019)
We add smart clique merging and fix a number of bugs in COSMO and the MOIWrapper.

- `94babde` Add Documentation for Chordal Decomposition
- `3ba43dd` Renaming merge concepts, adding docstrings
- `96e9c1c` Implementation of a SparseCoLO like transformation
- `752117b` Add unit test for psd completion with clique merge
- `975f87a` 🔧 Fix bug about recovering the dual variable
- `e524b59` Working recomputation of clique tree
- `09dc024` Implement clique merging
- `f90cc62` <#96> Map var idxs when processing MOI constraints
- `9f1b5ff` <#97> Loop over all MOI-Attributes during copy_to
- `59f0480` Add unit test for 1d psd projection
- `ab821d8` Fix bug in 1D PSD cone projection
- `857bf18` Merge pull request #98 from rschwarz/rs/moi095
- `969cb2b` Fix MOI tests "result_count".
- `c3d0263` Exclude MOI unit test "number_threads".
- `ed25b4b` [#94] Fix bug related to DualExpCone problems
- `7f91e71` 💾 Implement memory-efficient infeasibility checks
- `dc75bc5` 🔧 Code refactoring
- `a754b0c` Implement memory-efficient residual computation
- `020ada4` [#94] Merge pull request #94 from oxfordcontrol/check_psd_clean
- `4b8396d` Use (faster) in-place versions for in_pos_def()
- `5565fa9` Remove sparse version of is_pos_def!()
- `dceea67` Avoids eigendecomposition in psd tests.
- `340ecde` Assemble entry selector matrix more efficient


## Version 0.5.0 (30. August 2019)

- `8986827` Update to support MOI v0.9

## Version 0.4.4 (23. August 2019)

- `b514287` Fix memory issue in MOI's interface (PR #87)
- `cc32a30` Remove parameter printing for sets
- `d31e501` Code refactoring in chordal decomposition
- `128b409` Remove num of elem in A and size of b from header
- `b267b1c` Remove eigendecomposition of P in MOI wrapper
- `21eb8da` Move PowerCone to typed set in MOI wrapper test

## Version 0.4.3 (3. June 2019)
We add the exponential cone, the power cone and their dual cones. Furthermore, we offer the option to use different linear
system solvers. The most important changes are listed here:

- `eb7cb68` Add power cone and dual exp and pow cone
- `f3b5c2f` Add a logistic regression example
- `44cc8b3` Implement exponential cone constraints
- `4c3428c` Merge pull request #65 from oxfordcontrol/pg/lin_solvers
- `863168a` <#65> Support custom linear system solver

## Version 0.4.2 (5. May 2019)
This is a major upgrad since it adds the chordal decomposition functionality. The most important changes are listed here:
- `8dd70eb` Make sure graphs are connected
- `b44364c` Add chordal decomp info to header printing
- `5acaf73` Move PSDTriangle memory allocation after decomp
- `25593b4` Remove graph computation
- `a115196` Add DensePsdCones to handle non-decomposable cones
- `dbe70bc` Add PSD completion
- `4262bfe` Add chordal decomposition algorithm and unit tests
- `bdebbc5` 🔧 Fix bug in scaling of MOI-PSDTriangle constraint
- `726d03b` <#55> Implement set merging in MOI Wrapper
- `2a3957f` 🤫 Add support for MOI.Silent
- `d0f0e44` Fix bugs in warm starting
- `f5bd31b` 32 bit BLAS support
- `aeab91c` Check positive semidefiniteness of P in MOIWrapper
- `9ad8f36` Fix bug in MOI termination codes
- `34ebb89` Add deepcopy step for convex sets in set!()


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
