# Contributing

Contributions are always welcome:

* If you want to contribute features, bug fixes, etc, please take a look at our __Code Style Guide__ below
* Please report any issues and bugs that you encounter in [Issues](https://github.com/oxfordcontrol/ossdp/issues)
* As an open source project we are also interested in any projects and applications that use COSMO. Please let us know via email to: michael.garstka[at]eng.ox.ac.uk

## Code Style Guide

The code in this repository follows the naming and style conventions of [Julia Base](https://docs.julialang.org/en/v1.0/manual/style-guide/#Style-Guide-1) with a few modifications. This style guide is heavily "inspired" by the guides of [John Myles White](https://github.com/johnmyleswhite/Style.jl) and [JuMP](http://www.juliaopt.org/JuMP.jl/latest/style).

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