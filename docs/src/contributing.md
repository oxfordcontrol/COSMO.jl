# Contributing

This package is currently in maintenance mode. 

We are trying to keep it compatible with new releases of JuMP/MOI. Contributions for bug fixes or new features are always welcome. Enhancement ideas are tagged as such in the [Issues](https://github.com/oxfordcontrol/ossdp/issues) section. If you are unsure how to develop, continue with the section below.


# How-to develop
Let's assume you found a bug that you want to fix or develop a new feature. You have already downloaded the latest julia version as well as git. Follow these steps to make changes, test the changes, and open your first PR. This assumes that you don't yet have push access to the repository and need to fork.

- Create a new fork in Github by clicking on `Fork` at the top right of the repository page and name it e.g. `migarstka/COSMO.jl`. (If you have repository push rights, you don't need to fork and can clone the repository directly.)

- Clone the forked repository to your machine with:
```bash
$ git clone git@github.com:migarstka/COSMO.jl.git
```

- Navigate into the new folder and create a new branch:
```bash
$ git checkout -b mg/bug_fix
```

- Make your code changes by modifying the files in the repository.

- As a first sanity check that nothing obvious is broken, you can run one of the examples in the `examples/` folder. Open julia and navigate into the git repository. Next, we want to activate and instantiate the repositories virtual environment. In the julia REPL type `]`, then to activate the current environment (i.e. current `Project.toml` file):

```julia
(@v1.8) pkg>  activate .
```
Next, instantiate the environment, i.e. install the solver's dependencies defined in `Project.toml`:

```julia
(COSMO) pkg> instantiate
```

- Since we instantiated the local envrionment, `COSMO`  now refers to the local version that includes your changes. You can solve an example problem using the modified code:
```julia
julia> include("examples/qp.jl")
```

Let's assume you have made some changes and want to check whether the bug is fixed or the new feature broke anything. In both cases you should add some new tests for the change to `/test/UnitTest` and `/test/run_cosmo_tests.jl`. If you want to run existing tests follow the steps below.

### How-to run tests
You can run all package-native tests with the main test file:
```julia
julia> include("test/run_cosmo_tests.jl")
```
The full test suite also includes several MathOptInterface test problems that you can run with:
```julia
julia> include("test/runtests.jl")
```
These take longer to complete, so I recommend using `test/run_cosmo_tests.jl` during iterative development and then switch to 
`test/runtests.jl` once the previous tests pass.

### Open a PR
Once all tests pass, you have added additional tests for your changes, and used docstrings to describe any new functionality in the code, you can open a PR. Commit your changes with a detailed commit message, and push your branch to the (fork) remote.
```bash
$ git add .
$ git commit 
[Write commit message]
$ git push
```
Then navigate to [https://github.com/oxfordcontrol/COSMO.jl/compare](https://github.com/oxfordcontrol/COSMO.jl/compare) and create a new PR with `base repository: oxfordcontrol/COSMO.jl` `base: master` and `head repository: [your-fork-repo]` `base:mg/bug_fix`.
Ask me or other collaborators to review it, before it can get merged.

## Making a new release
Changes to `master` are periodically bundled in a new release version of the package. If you have permissions, you can make a new release by following these steps:

1. Ensure all tests pass (locally and in Github Actions)
2. Locally run `bumpversion patch` (or `bumpversion minor` or `bumpversion major` depending on [semver](https://semver.org/) convention) to increment the version number across the repository (most importantly in `Project.toml`).
3. Check and commit the changes and use the commit message `Bump version to X`. 
4. Push changes to remote with `git push`.
5. Let the julia package registry know about the new release by commenting `@JuliaRegistrator register` on the bump commit (like [here](https://github.com/oxfordcontrol/COSMO.jl/commit/20764ba075e1f598ec17990fb1721dcb5a3b418b#commitcomment-141674207)).
6. This will open a PR in the `JuliaRegistries` repository and a new release will typically be approved within hours.
7. Update the [CHANGELOG](https://github.com/oxfordcontrol/COSMO.jl/blob/master/CHANGELOG.md) and describe the changes bundled in the new release.

## Updating COSMO's documentation
The code for COSMO's documentation resides in the `/docs` folder of the repository and is created with the [Documenter.jl package](https://documenter.juliadocs.org/stable/). There are two main versions of this documentation. **Stable** is based on the committed changes of the latest tagged release. **Dev** is based on the latest commit of `master`.

The documentation has its own environment defined inside `/docs` (see `docs/Project.toml` file). To edit the documentation, edit the files in `/docs`. It's advised to run the documentation generation process locally do ensure the layout is as expected and that `Literate` code examples are built without errors. To generate the documentation locally, follow these steps:

- Clone the repository and make your changes (as described in the section above).
- Navigate into the `/docs` folder and start julia.
- Activate and instantiate the (docs) environment, after `]`, type:

```julia 
(@v1.8) pkg> activate .
```

```julia 
(docs) pkg> instantiate
```

- Build the documentation locally into the `/build` folder with

```julia 
include("make.jl")
```

(The first build will take several minutes. Subsequent changes will take ~30s.)

- After the new documentation has been built, you can serve the files in `/build` in your browser, e.g. using a python webserver. Navigate into `/build` and run:
```bash
python3 -m http.server --bind localhost
```

- View the local documentation in your browser of choice at `http://127.0.0.1:8000/`

- If you are satisfied with the results, open a PR for the changes as described in the section above.



# Code Style Guide

The code in this repository follows the naming and style conventions of [Julia Base](https://docs.julialang.org/en/v1.0/manual/style-guide/#Style-Guide-1) with a few modifications. This style guide is heavily "inspired" by the guides of [John Myles White](https://github.com/johnmyleswhite/Style.jl) and [JuMP](https://jump.dev/JuMP.jl/stable/developers/style/).

### Formatting
* Use one tab when indenting a new block (except `module`)

* Use spaces between operators, except for `^`, `'`, and `:`
* Use single space after commas and semicolons
* Don't use spaces around parentheses, or braces

**Bad**: `f(x,y) = [5*sin(x+y);y']` **Good**: `f(x, y) = [5 * sin(x + y); y']`
* Use spacing with keyword arguments

**Bad**: `foo(x::Integer=1)` **Good**: `foo(x::Integer = 1)`

* Don't parenthesize conditions

**Bad**: `if (a == b)` **Good**: `if a == b`
### Naming
* Modules and Type names use capitilization and camel case, e.g. `module LinearAlgebra`, `struct ConvexSets`.
* Functions are lowercase and use underscores to seperate words, e.g. `has_key(x)`, `is_valid(y)`.
* Normal variables are lowercase and use underscores like functions, e.g. `convex_set`
* Constants are uppercase, e.g. `const MY_CONSTANT`
* **Always** append `!` to names of functions that modify their arguments.
* Function arguments that are mutated come first. Otherwise follow the rules layed out in Julia Base [Argument ordering](https://docs.julialang.org/en/v1.0/manual/style-guide/#Write-functions-with-argument-ordering-similar-to-Julia-Base-1)
* Files are named like functions, e.g. `my_new_file.jl`
### Syntax
* Use `1.0` instead of `1.`

### Git(hub)-specific conventions
* Branch names should be prepended with the initials of the creator and a forward slash, e.g. `mg/newIdea` instead of `newIdea`
* Commit messages should have the following format:
```
<#IssueId> Short (72 chars or less) summary

More detailed explanatory text. Wrap it to 72 characters. The blank
line separating the summary from the body is critical.

Imperative style for the commit message: "Fix bug" and not "Fixed
bug" or "Fixes bug."

The issue id can be ommitted if the commit does not related to a specific open issue
```