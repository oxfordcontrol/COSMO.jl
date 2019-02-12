var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "COSMO.jl is a Julia implementation of the Conic Operator Splitting Method. The underlying ADMM-algortihm is well-suited for large convex conic problems. COSMO solves the following problem:beginarrayll mboxminimize  textstylefrac12x^top Px + q^top x mboxsubject to  Ax + s = b   s in mathcalC endarraywith decision variables x in mathrmR^n, s in mathrmR^m and data matrices P=P^top succeq 0, q in mathrmR^n, A in mathrmR^m times n, and b in mathrmR^m. The convex set mathcalC  is a composition of convex sets and cones."
},

{
    "location": "#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "COSMO solves linear programs, quadratic programs, second-order cone programs and semidefinite programs\nSemi-definite programs with quadratic objective functions are natively supported\nInfeasible problems are detected without a homogeneous self-dual embedding of the problem\nSupport for MathOptInterface and the upcoming JuMP v0.19 release, which allows you to describe your problem in JuMP.\nCOSMO tries to decompose large structured PSD constraints using chordal decomposition techniques. This often results in a significant speedup compared to the original problem. (this feature is not yet available in Julia v1.0)\nCOSMO supports warm-starting of the decision variables."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "COSMO can be installed using the Julia package manager for Julia v1.0 and higher. Inside the Julia REPL, type ] to enter the Pkg REPL mode then runpkg> add COSMOIf you want to install the latest version from master runpkg> add COSMO#master"
},

{
    "location": "#Quick-Example-1",
    "page": "Home",
    "title": "Quick Example",
    "category": "section",
    "text": "Consider the following 2x2 semidefinite program with decision variable X:beginarrayll mboxminimize   texttr(CX)\nmboxsubject to   texttr(A X) = b \n                    X succeq 0\nendarraywith problem data A, b and C:A = beginbmatrix 1  0  5  2endbmatrix\nC = beginbmatrix 1  2  0  2endbmatrix\nb = 4where tr denotes the trace of a matrix. We can solve this problem either using COSMO\'s interface:using COSMO, LinearAlgebra\n\nC =  [1. 2; 0 2]\nA = [1. 0; 5 2]\nb = 4;\n\nmodel = COSMO.Model();\n\n# define the cost function\nP = zeros(4, 4)\nq = vec(C)\n\n# define the constraints\n# A x = b\ncs1 = COSMO.Constraint(vec(A)\', -b, COSMO.ZeroSet)\n# X in PSD cone\ncs2 = COSMO.Constraint(Matrix(1.0I, 4, 4), zeros(4), COSMO.PsdCone)\nconstraints = [cs1; cs2]\n\n# assemble and solve the model\nassemble!(model, P, q, constraints)\nresult = COSMO.optimize!(model);\n\nX_sol = reshape(result.x, 2, 2)\nobj_value = result.obj_valor we can describe the problem using JuMP and use COSMO as the backend solver:using COSMO, JuMP\n\nC =  [1 2; 0 2]\nA = [1 0; 5 2]\nb = 4;\n\nm = Model(with_optimizer(COSMO.Optimizer));\n@variable(m, X[1:2, 1:2], PSD)\n@objective(m, Min, tr(C * X));\n@constraint(m, tr(A * X) == b);\nstatus = JuMP.optimize!(m);\n\nX_sol = JuMP.value.(X)\nobj_value = JuMP.objective_value(m)"
},

{
    "location": "#Credits-1",
    "page": "Home",
    "title": "Credits",
    "category": "section",
    "text": "The following people are involved in the development of COSMO:Michael Garstka (main development)\nNikitas Rontsis (algorithm performance)\nPaul Goulart (code architecture, maths and algorithms)\nMark Cannon (maths and algorithms)*all contributors are affiliated with the University of Oxford.If this project is useful for your work please considerCiting the relevant paper\nLeaving a star on the GitHub repository"
},

{
    "location": "#Licence-1",
    "page": "Home",
    "title": "Licence",
    "category": "section",
    "text": "COSMO.jl is licensed under the Apache License 2.0. For more details click here."
},

{
    "location": "guide/#",
    "page": "User Guide",
    "title": "User Guide",
    "category": "page",
    "text": ""
},

{
    "location": "guide/#User-Guide-1",
    "page": "User Guide",
    "title": "User Guide",
    "category": "section",
    "text": "This user guide describes the basic structures and functions to define an optimisation problem, to solve the problem and to analyse the result. If you want to use JuMP to describe the problem, see the JuMP Interface section."
},

{
    "location": "guide/#Model-1",
    "page": "User Guide",
    "title": "Model",
    "category": "section",
    "text": "The problem data, user settings and workspace variables are all stored in a Model. To get started define an empty model:m = COSMO.Model()To initialize the model with an optimisation problem we need to define three more things:the objective function, i.e. the matrix P and the vector q in frac12x^top P x + q^top x\nan array of constraints\na Settings object that specifies how COSMO solves the problem (optional)"
},

{
    "location": "guide/#Objective-Function-1",
    "page": "User Guide",
    "title": "Objective Function",
    "category": "section",
    "text": "To set the objective function of your optimisation problem simply define the square positive semidefinite matrix P in mathrmR^ntimes n and the vector q in mathrmR^n. You might have to transform your optimisation problem for this step."
},

{
    "location": "guide/#COSMO.Constraint",
    "page": "User Guide",
    "title": "COSMO.Constraint",
    "category": "type",
    "text": "Constraint(A, b, convex_set, dim = 0, indices = 0:0)\n\nCreates a COSMO constraint: Ax + b ∈ convex_set.\n\nBy default the following convex sets are supported: ZeroSet, Nonnegatives, SecondOrderCone, PsdCone, PsdConeTriangle.\n\nExamples\n\njulia> Constraint([1 0;0 1], zeros(2), COSMO.Nonnegatives)\nConstraint\nSize of A: (2, 2)\nConvexSet: COSMO.Nonnegatives\n\n\n\nThe optional arguments dim and indices can be used to specify A and b for subparts of variable x. If x has dimension dim = 4, then x[2] and x[3] can be constrained to the zero cone in the following way:\n\nExamples\n\njulia> c = Constraint([1 0;0 1], zeros(2), COSMO.ZeroSet, 4, 2:3)\nConstraint\nSize of A: (2, 4)\nConvexSet: COSMO.ZeroSet\n\nNotice that extra columns of A have been added automatically.\n\njulia>Matrix(c.A)\n2×4 Array{Float64,2}:\n0.0  1.0  0.0  0.0\n0.0  0.0  1.0  0.0\n\n\n\n\n\n"
},

{
    "location": "guide/#Constraints-1",
    "page": "User Guide",
    "title": "Constraints",
    "category": "section",
    "text": "COSMO.Constraint"
},

{
    "location": "guide/#COSMO.ZeroSet",
    "page": "User Guide",
    "title": "COSMO.ZeroSet",
    "category": "type",
    "text": "ZeroSet(dim)\n\nCreates the zero set  0 ^dim of dimension dim. If x ∈ ZeroSet then all entries of x are zero.\n\n\n\n\n\n"
},

{
    "location": "guide/#COSMO.Nonnegatives",
    "page": "User Guide",
    "title": "COSMO.Nonnegatives",
    "category": "type",
    "text": "Nonnegatives(dim)\n\nCreates the nonnegative orthant  x in mathbbR^dim  x ge 0   of dimension dim.\n\n\n\n\n\n"
},

{
    "location": "guide/#COSMO.SecondOrderCone",
    "page": "User Guide",
    "title": "COSMO.SecondOrderCone",
    "category": "type",
    "text": "SecondOrderCone(dim)\n\nCreates the second-order cone (or Lorenz cone)  (tx) in mathrmR^dim   x _2  leq t .\n\n\n\n\n\n"
},

{
    "location": "guide/#COSMO.PsdCone",
    "page": "User Guide",
    "title": "COSMO.PsdCone",
    "category": "type",
    "text": "PsdCone(dim)\n\nCreates the cone of symmetric positive semidefinite matrices mathcalS_+^dim The entries of the matrix X are stored column-by-column in the vector x of dimension dim AccordinglyX \\in \\mathbb{S}+ \\Rightarrow x \\in \\mathcal{S}+^{dim} whereX = \\text{mat}(x)``.\n\n\n\n\n\n"
},

{
    "location": "guide/#COSMO.PsdConeTriangle",
    "page": "User Guide",
    "title": "COSMO.PsdConeTriangle",
    "category": "type",
    "text": "PsdConeTriangle(dim)\n\nCreates the cone of symmetric positive semidefinite matrices. The entries of the upper-triangular part of matrix X are stored in the vector x of dimension dim. A r times r matrix has r(r+1)2 upper triangular elements and results in a vector of mathrmdim = r(r+1)2.\n\nExamples\n\nThe matrix\n\nbeginbmatrix x_1  x_2  x_4 x_2  x_3  x_5 x_4  x_5  x_6 endbmatrix\n\nis transformed to the vector x_1 x_2 x_3 x_4 x_5 x_6^top with corresponding constraint  PsdConeTriangle(6).\n\n\n\n\n\n"
},

{
    "location": "guide/#Convex-Sets-1",
    "page": "User Guide",
    "title": "Convex Sets",
    "category": "section",
    "text": "COSMO.ZeroSet\nCOSMO.Nonnegatives\nCOSMO.SecondOrderCone\nCOSMO.PsdCone\nCOSMO.PsdConeTriangle"
},

{
    "location": "guide/#Settings-1",
    "page": "User Guide",
    "title": "Settings",
    "category": "section",
    "text": "Settings can be specified using the COSMO.Settings struct. The following settings are available:Argument Description Values (default)\nrho ADMM rho step 0.1\nsigma ADMM sigma step 1e-6\nalpha Relaxation parameter 1.6\neps_abs Absolute residual tolerance 1e-4\neps_rel Relative residual tolerance 1e-4\neps_prim_inf Primal infeasibility tolerance 1e-4\neps_dual_inf Dual infeasibility tolerance 1e-4\nmax_iter Maximum number of iterations 2500\nverbose Verbose printing false\nverbose_timing Verbose timing false\ncheck_termination Check termination interval 40\ncheck_infeasibility Check infeasibility interval 40\nscaling Number of scaling iterations 10\nadaptive_rho Automatic adaptation of step size parameter true\ntime_limit set solver time limit in s 0.0For more low-level settings, see the Settings type definition in /src/settings.jl."
},

{
    "location": "guide/#Results-1",
    "page": "User Guide",
    "title": "Results",
    "category": "section",
    "text": "After attempting to solve the problem, COSMO will return a result object with the following fields:Fieldname Type Description\nx Vector{Float64} Primal variable\ny Vector{Float64} Dual variable\ns Vector{Float64} (Primal) set variable\nobj_val Float64 Objective value\niter Int64 Number of iterations\nstatus Symbol Solution status\ninfo COSMO.ResultInfo Struct with more information\ntimes COSMO.ResultTimes Struct with several measured times"
},

{
    "location": "guide/#Status-Codes-1",
    "page": "User Guide",
    "title": "Status Codes",
    "category": "section",
    "text": "COSMO will return one of the following statuses:Status Code Description\n:Solved A optimal solution was found\n:Unsolved Default value\n:Max_iter_reached Solver reached iteration limit (set with Settings.max_iter)\n:Time_limit_reached Solver reached time limit (set with Settings.time_limit)\n:Primal_infeasible Problem is primal infeasible\n:Dual_infeasible Problem is dual infeasible"
},

{
    "location": "guide/#Timings-1",
    "page": "User Guide",
    "title": "Timings",
    "category": "section",
    "text": "If settings.verbose_timing is set to true, COSMO will report the following times in result.times:Time Name Description\nsolver_time Total time used to solve the problem\nsetup_time Setup time = graph_time + factor_time\ngraph_time Time used to perform chordal decomposition\nfactor_time Time used to factor the system of linear equations\niter_time Time spent in iteration loop\nproj_time Time spent in projection functions\npost_time Time used for post processingIt holds: solver_time = setup_time+ iter_time + post_time, setup_time = graph_time+ factor_time, proj_time is a subset of iter_time."
},

{
    "location": "jump/#",
    "page": "JuMP Interface",
    "title": "JuMP Interface",
    "category": "page",
    "text": ""
},

{
    "location": "jump/#JuMP-Interface-1",
    "page": "JuMP Interface",
    "title": "JuMP Interface",
    "category": "section",
    "text": "Our JuMP interface allows you to describe and modify your optimisation problem with JuMP and use COSMO as the backend solver. The interface is defined in /src/MOIWrapper.jl.note: Note\nCOSMO requires the upcoming JuMP v0.19 release that is based on the MathOptInterface package. Until this version is released we recommend using the latest beta version which can be downloaded via the Julia package manager withpkg> add JuMP#v0.19-beta2"
},

{
    "location": "jump/#Use-COSMO-1",
    "page": "JuMP Interface",
    "title": "Use COSMO",
    "category": "section",
    "text": "To specify COSMO as the solver for your JuMP model, load the solver module with using COSMO and then use the with_optimizer() function when initialising the JuMP model:m = Model(with_optimizer(COSMO.Optimizer);"
},

{
    "location": "jump/#Pass-Solver-Settings-1",
    "page": "JuMP Interface",
    "title": "Pass Solver Settings",
    "category": "section",
    "text": "Solver-specific settings can be passed after the COSMO.Optimizer object. For example, if you want to adjust the maximum number of iterations and turn on verbose printing usem = Model(with_optimizer(COSMO.Optimizer, max_iter = 5000, verbose = true);The full list of available settings can be found in the Settings section."
},

{
    "location": "jump/#Results-1",
    "page": "JuMP Interface",
    "title": "Results",
    "category": "section",
    "text": "After solving the problem the result can be obtained using the standard JuMP commands. To see if the optimisation was successful useJuMP.termination_status(m)\nJuMP.primal_status(m)If a solution is available, the optimal objective value can be queried usingJuMP.objective_value(m)and the value of a decision variable x can be obtained withJuMP.value.(x)For more information on how to use JuMP check the JuMP documentation."
},

{
    "location": "method/#",
    "page": "Method",
    "title": "Method",
    "category": "page",
    "text": ""
},

{
    "location": "method/#Method-1",
    "page": "Method",
    "title": "Method",
    "category": "section",
    "text": "Some information about the underlying ADMM algorithm will be added soon."
},

{
    "location": "method/#Algorithm-1",
    "page": "Method",
    "title": "Algorithm",
    "category": "section",
    "text": ""
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "Some example problems that are solved with COSMO can be found in the /examples folder.In the first example the native solver interface is used to define and solve a Linear Program. In the second example JuMP is used to describe the problem and COSMO is set as the solver backend."
},

{
    "location": "examples/#Linear-Program-1",
    "page": "Examples",
    "title": "Linear Program",
    "category": "section",
    "text": "We want to solve the following linear program with decision variable x:beginarrayll mboxminimize   c^top x\nmboxsubject to   A x leq b \n                    x geq 1 \n                    x_2 geq 5 \n                    x_1 + x_3 geq 4\nendarrayThe problem can be solved with COSMO in the following way:using COSMO, LinearAlgebra, SparseArrays, Test\n\nc = [1; 2; 3; 4.]\nA = Matrix(1.0I, 4, 4)\nb = [10; 10; 10; 10]\nn = 4\n# -------------------\n# create constraints A * x + b in set\n# -------------------\n# Ax <= b\nc1 = COSMO.Constraint(-A, b, COSMO.Nonnegatives)\n# x >= 1\nc2 = COSMO.Constraint(Matrix(1.0I, n, n), -ones(n), COSMO.Nonnegatives)\n# x2 >= 5\nc3 = COSMO.Constraint(1, -5, COSMO.Nonnegatives, n, 2:2)\n# x1 + x3 >= 4\nc4 = COSMO.Constraint([1 0 1 0], -4, COSMO.Nonnegatives)\n\n# -------------------\n# define cost function\n# -------------------\nP = spzeros(4, 4)\nq = c\n\n# -------------------\n# assemble solver model\n# -------------------\nsettings = COSMO.Settings(max_iter=2500, verbose=true, eps_abs = 1e-4, eps_rel = 1e-5)\nmodel = COSMO.Model()\nassemble!(model, P, q, [c1; c2; c3; c4], settings)\nres = COSMO.optimize!(model);\n\n@testset \"Linear Problem\" begin\n  @test isapprox(res.x[1:4], [3; 5; 1; 1], atol=1e-2, norm = (x -> norm(x, Inf)))\n  @test isapprox(res.obj_val, 20.0, atol=1e-2)\nend"
},

{
    "location": "examples/#Closest-Correlation-Matrix-1",
    "page": "Examples",
    "title": "Closest Correlation Matrix",
    "category": "section",
    "text": "We consider the problem of finding the closest correlation matrix X to a given random matrix C. With closest correlation matrix we mean a positive semidefinite matrix with ones on the diagonal. The problem is given by:beginarrayll mboxminimize   frac12X - C_F^2\nmboxsubject to   X_ii = 1 quad i=1dotsn \n                    X succeq 0\nendarrayNotice how JuMP is used to describe the problem. COSMO is chosen as the backend solver using JuMP\'s with_optimizer() function.using COSMO, JuMP, LinearAlgebra, SparseArrays, Test, Random\nrng = Random.MersenneTwister(12345);\n\n# create a random test matrix C\nn = 8\nC = -1 .+ rand(rng, n, n) .* 2;\nc = vec(C);\n\n# define problem in JuMP\nq = -vec(C);\nr = 0.5 * vec(C)\' * vec(C);\nm = Model(with_optimizer(COSMO.Optimizer, verbose=true, eps_abs = 1e-4));\n@variable(m, X[1:n, 1:n], PSD);\nx = vec(X);\n@objective(m, Min, 0.5 * x\' * x  + q\' * x + r)\nfor i = 1:n\n  @constraint(m, X[i, i] == 1.)\nend\n\n# solve and get results\nstatus = JuMP.optimize!(m)\nobj_val = JuMP.objective_value(m)\nX_sol = JuMP.value.(X)"
},

{
    "location": "citing/#",
    "page": "Citing COSMO",
    "title": "Citing COSMO",
    "category": "page",
    "text": ""
},

{
    "location": "citing/#Citing-1",
    "page": "Citing COSMO",
    "title": "Citing",
    "category": "section",
    "text": "If you find COSMO useful in your project, we kindly request that you cite the following paper:@article{garstka_2019,\n  author        = {Michael Garstka and Mark Cannon and Paul Goulart},\n  title         = {{COSMO}: A conic operator splitting method for convex conic problems},\n  journal       = {arXiv e-prints},\n  year          = {2019},\n  month         = jan,\n  archiveprefix = {arXiv},\n  eprint        = {1901.10887},\n  keywords      = {Mathematics - Optimization and Control},\n  primaryclass  = {math.OC},\n  url           = {https://arxiv.org/abs/1901.10887},\n}A preprint can be downloaded here."
},

{
    "location": "contributing/#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing/#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "Contributions are always welcome:If you want to contribute features, bug fixes, etc, please take a look at our Code Style Guide below\nPlease report any issues and bugs that you encounter in Issues\nAs an open source project we are also interested in any projects and applications that use COSMO. Please let us know via email to: michael.garstka[at]eng.ox.ac.uk"
},

{
    "location": "contributing/#Code-Style-Guide-1",
    "page": "Contributing",
    "title": "Code Style Guide",
    "category": "section",
    "text": "The code in this repository follows the naming and style conventions of Julia Base with a few modifications. This style guide is heavily \"inspired\" by the guides of John Myles White and JuMP."
},

{
    "location": "contributing/#Formatting-1",
    "page": "Contributing",
    "title": "Formatting",
    "category": "section",
    "text": "Use one tab when indenting a new block (except module)\nUse spaces between operators, except for ^, \', and :\nUse single space after commas and semicolons\nDon\'t use spaces around parentheses, or bracesBad: f(x,y) = [5*sin(x+y);y\'] Good: f(x, y) = [5 * sin(x + y); y\']Use spacing with keyword argumentsBad: foo(x::Integer=1) Good: foo(x::Integer = 1)Don\'t parenthesize conditionsBad: if (a == b) Good: if a == b"
},

{
    "location": "contributing/#Naming-1",
    "page": "Contributing",
    "title": "Naming",
    "category": "section",
    "text": "Modules and Type names use capitilization and camel case, e.g. module LinearAlgebra, struct ConvexSets.\nFunctions are lowercase and use underscores to seperate words, e.g. has_key(x), is_valid(y).\nNormal variables are lowercase and use underscores like functions, e.g. convex_set\nConstants are uppercase, e.g. const MY_CONSTANT\nAlways append ! to names of functions that modify their arguments.\nFunction arguments that are mutated come first. Otherwise follow the rules layed out in Julia Base Argument ordering\nFiles are named like functions, e.g. my_new_file.jl"
},

{
    "location": "contributing/#Syntax-1",
    "page": "Contributing",
    "title": "Syntax",
    "category": "section",
    "text": "Use 1.0 instead of 1."
},

{
    "location": "contributing/#Git(hub)-specific-conventions-1",
    "page": "Contributing",
    "title": "Git(hub)-specific conventions",
    "category": "section",
    "text": "Branch names should be prepended with the initials of the creator and a forward slash, e.g. mg/newIdea instead of newIdea\nCommit messages should have the following format:<#IssueId> Short (72 chars or less) summary\n\nMore detailed explanatory text. Wrap it to 72 characters. The blank\nline separating the summary from the body is critical.\n\nImperative style for the commit message: \"Fix bug\" and not \"Fixed\nbug\" or \"Fixes bug.\"\n\nThe issue id can be ommitted if the commit does not related to a specific open issue"
},

{
    "location": "api/#",
    "page": "API Reference",
    "title": "API Reference",
    "category": "page",
    "text": ""
},

{
    "location": "api/#API-Reference-1",
    "page": "API Reference",
    "title": "API Reference",
    "category": "section",
    "text": ""
},

{
    "location": "api/#COSMO.assemble!",
    "page": "API Reference",
    "title": "COSMO.assemble!",
    "category": "function",
    "text": "assemble!(model, P, q, constraint(s), [settings, x0, y0, s0])\n\nAssembles a COSMO.Model with a cost function defind by P and q, and a number of constraints.\n\nThe positive semidefinite matrix P and vector q are used to specify the cost function of the optimization problem:\n\nmin   1/2 x\'Px + q\'x\ns.t.  Ax + b ∈ C\n\nconstraints is a COSMO.Constraint or an array of COSMO.Constraint objects that are used to describe the constraints on x.\n\n\n\nThe optinal arguments x0, s0, and y0 can be used to provide the solver with warm starting values for the primal variable x, the primal slack variable s and the dual variable y. The optinal argument settings can be used to pass custom solver settings.\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.set!",
    "page": "API Reference",
    "title": "COSMO.set!",
    "category": "function",
    "text": "set!(model, P, q, A, b, convex_sets, [settings])\n\nSets model data directly based on provided fields.\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.warm_start_primal!",
    "page": "API Reference",
    "title": "COSMO.warm_start_primal!",
    "category": "function",
    "text": "warm_start_primal!(model, x0, [ind])\n\nProvides the COSMO.Model with warm starting values for the primal variable x. ind can be used to warm start certain components of x.\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.warm_start_slack!",
    "page": "API Reference",
    "title": "COSMO.warm_start_slack!",
    "category": "function",
    "text": "warm_start_slack!(model, s0, [ind])\n\nProvides the COSMO.Model with warm starting values for the primal slack variable s. ind can be used to warm start certain components of s.\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.warm_start_dual!",
    "page": "API Reference",
    "title": "COSMO.warm_start_dual!",
    "category": "function",
    "text": "warm_start_dual!(model, y0, [ind])\n\nProvides the COSMO.Model with warm starting values for the dual variable y. ind can be used to warm start certain components of y.\n\n\n\n\n\n"
},

{
    "location": "api/#Model-1",
    "page": "API Reference",
    "title": "Model",
    "category": "section",
    "text": "COSMO.assemble!\nCOSMO.set!\nCOSMO.warm_start_primal!\nCOSMO.warm_start_slack!\nCOSMO.warm_start_dual!"
},

{
    "location": "api/#COSMO.Settings",
    "page": "API Reference",
    "title": "COSMO.Settings",
    "category": "type",
    "text": "Settings(;arg=val)\n\nCreates a COSMO settings object that is used to pass user settings to the solver.\n\nArgument Description Values (default)\nrho ADMM rho step 0.1\nsigma ADMM sigma step 1e-6\nalpha Relaxation parameter 1.6\neps_abs Absolute residual tolerance 1e-4\neps_rel Relative residual tolerance 1e-4\neps_prim_inf Primal infeasibility tolerance 1e-4\neps_dual_inf Dual infeasibility tolerance 1e-4\nmax_iter Maximum number of iterations 2500\nverbose Verbose printing false\nverbose_timing Verbose timing false\ncheck_termination Check termination interval 40\ncheck_infeasibility Check infeasibility interval 40\nscaling Number of scaling iterations 10\nadaptive_rho Automatic adaptation of step size parameter true\ntime_limit set solver time limit in s 0\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.Result",
    "page": "API Reference",
    "title": "COSMO.Result",
    "category": "type",
    "text": "Result{T <: AbstractFloat}\n\nObject returned by the COSMO solver after calling optimize!(model). It has the following fields:\n\nFieldname Type Description\nx Vector{T} Primal variable\ny Vector{T} Dual variable\ns Vector{T} (Primal) set variable\nobj_val T Objective value\niter Int64 Number of iterations\nstatus Symbol Solution status\ninfo COSMO.ResultInfo Struct with more information\ntimes COSMO.ResultTimes Struct with several measured times\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.Model",
    "page": "API Reference",
    "title": "COSMO.Model",
    "category": "type",
    "text": "Model()\n\nInitializes an empty COSMO model that can be filled with problem data using assemble!(model, P, q,constraints, [settings]).\n\n\n\n\n\n"
},

{
    "location": "api/#Types-1",
    "page": "API Reference",
    "title": "Types",
    "category": "section",
    "text": "COSMO.Settings\nCOSMO.Result\nCOSMO.Model"
},

]}
