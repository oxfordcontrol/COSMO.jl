var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "COSMO.jl is a Julia implementation of the Conic Operator Splitting Method. The underlying ADMM-algorithm is well-suited for large convex conic problems. COSMO solves the following problem:beginarrayll mboxminimize  textstylefrac12x^top Px + q^top x mboxsubject to  Ax + s = b   s in mathcalK endarraywith decision variables x in mathbbR^n, s in mathbbR^m and data matrices P=P^top succeq 0, q in mathbbR^n, A in mathbbR^m times n, and b in mathbbR^m. The convex set mathcalK  is a composition of convex sets and cones."
},

{
    "location": "#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "Versatile: COSMO solves linear programs, quadratic programs, second-order cone programs and semidefinite programs\nQuad SDPs: Positive semidefinite programs with quadratic objective functions are natively supported\nInfeasibility detection: Infeasible problems are detected without a homogeneous self-dual embedding of the problem\nJuMP support: COSMO supports MathOptInterface and JuMP v0.19, which allows you to describe your problem in JuMP\nChordal decomposition: COSMO tries to decompose large structured PSD constraints using chordal decomposition techniques. This often results in a significant speedup compared to the original problem. (this feature is not yet available in Julia v1.0)\nWarm starting: COSMO supports warm starting of the decision variables\nOpen Source: Our code is available on GitHub and distributed under the Apache 2.0 Licence"
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
    "text": "This user guide describes the basic structures and functions to define an optimisation problem, to solve the problem and to analyse the result. If you want to use JuMP to describe the problem, see the JuMP Interface section.COSMO solves optimisation problems in the following format:beginarrayll mboxminimize  textstylefrac12x^top Px + q^top x mboxsubject to  Ax + s = b   s in mathcalK endarraywith decision variables x in mathbbR^n, s in mathbbR^m and data matrices P=P^top succeq 0, q in mathbbR^n, A in mathbbR^m times n, and b in mathbbR^m. The convex set mathcalK  is a composition of convex sets and cones."
},

{
    "location": "guide/#Model-1",
    "page": "User Guide",
    "title": "Model",
    "category": "section",
    "text": "The problem data, user settings and workspace variables are all stored in a Model. To get started define an empty model:model = COSMO.Model()To initialize the model with an optimisation problem we need to define three more things:the objective function, i.e. the matrix P and the vector q in frac12x^top P x + q^top x\nan array of constraints\na Settings object that specifies how COSMO solves the problem (optional)"
},

{
    "location": "guide/#Objective-Function-1",
    "page": "User Guide",
    "title": "Objective Function",
    "category": "section",
    "text": "To set the objective function of your optimisation problem simply define the square positive semidefinite matrix P in mathrmR^ntimes n and the vector q in mathrmR^n. You might have to transform your optimisation problem into a solver compatible format for this step."
},

{
    "location": "guide/#Constraints-1",
    "page": "User Guide",
    "title": "Constraints",
    "category": "section",
    "text": "The COSMO interface expects constraints to have the form A_i x + b_i in mathcalK_i, where mathcalK_i is one of the convex sets defined below:Convex Set Description\nZeroSet(dim) The set  0 ^dim that contains the origin\nNonnegatives(dim) The nonnegative orthant  x in mathbbR^dim  x_i ge 0 forall i=1dotsmathrmdim \nSecondOrderCone(dim) The second-order (Lorenz) cone  (tx) in mathbbR^dim    x_2   leq t \nPsdCone(dim) The vectorized positive semidefinite cone mathcalS_+^dim. x is the vector obtained by stacking the columns of the positive semidefinite matrix X, i.e. X in mathbbS^sqrtdim_+ Rightarrow textvec(X) = x in mathcalS_+^dim\nPsdConeTriangle(dim) The vectorized positive semidefinite cone mathcalS_+^dim. x is the vector obtained by stacking the columns of the upper triangular part of the positive semidefinite matrix X, i.e. X in mathbbS^d_+ Rightarrow textsvec(X) = x in mathcalS_+^dim where d=sqrt14 + 2 textdim - 12The constructor for a constraint expects a matrix A, a vector b and a convex_set.Lets consider a problem with a decision variable x in mathbbR^5. Suppose we want to create the two constraint x_2 + 5 geq 0 and x_3 - 3 geq 0. We can do this either by creating two constraints and adding them to an array:  constraint1 = COSMO.Constraint([0.0 1.0 0.0 0.0 0.0], 5.0, COSMO.Nonnegatives)\n  constraint2 = COSMO.Constraint([0.0 0.0 1.0 0.0 0.0], -3.0, COSMO.Nonnegatives)\n  constraints = [constraint1; constraint2]The second option is to include both in one constraint:constraint1 = COSMO.Constraint([0.0 1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0], [5.0; -3.0], COSMO.Nonnegatives)Another way to construct the constraint is to used the optional arguments dim, the dimension of x, and indices, the elements of x that appear in the constraint. When specifying these arguments, A and b only refer to the elements of x in indices:constraint1 = COSMO.Constraint([1.0 0.0; 0.0 1.0], [5.0; -3.0], COSMO.Nonnegatives, 5, 2:3)Consider as a second example the positive semidefinite constraint on a matrix X in  mathbbS_+^3. Our decision variable is the vector x obtained by stacking the columns of X. We can specify the constraint on x in the following way:I_9 x + 0_9 in mathcalS_+^9or in Julia:constraint1 = COSMO.Constraint(Matrix(1.0I, 9, 9), zeros(9), COSMO.PsdCone)Several constraints can be combined in an array:constraints = [constraint_1, constraint_2, ..., constraint_N]"
},

{
    "location": "guide/#Settings-1",
    "page": "User Guide",
    "title": "Settings",
    "category": "section",
    "text": "The solver settings are stored in a Settings object and can be adjusted by the user. To create a Settings object just call the constructor:settings = COSMO.Settings()The settings object holds the following options and default values:Argument Description Values (default)\nrho ADMM rho step 0.1\nsigma ADMM sigma step 1e-6\nalpha Relaxation parameter 1.6\neps_abs Absolute residual tolerance 1e-4\neps_rel Relative residual tolerance 1e-4\neps_prim_inf Primal infeasibility tolerance 1e-4\neps_dual_inf Dual infeasibility tolerance 1e-4\nmax_iter Maximum number of iterations 2500\nverbose Verbose printing false\nverbose_timing Verbose timing false\ncheck_termination Check termination interval 40\ncheck_infeasibility Check infeasibility interval 40\nscaling Number of scaling iterations 10\nadaptive_rho Automatic adaptation of step size parameter true\ntime_limit set solver time limit in s 0.0To adjust those values, either pass your preferred option and parameter as a key-value pair to the constructor or edit the corresponding field afterwards. For example if you want to enable verbose printing and increase the solver accuracy, you can typesettings = COSMO.Settings(verbose = true, eps_abs = 1e-5, eps_rel = 1e-5)\n# the following is equivalent\nsettings = COSMO.Settings()\nsettings.verbose = true\nsettings.eps_abs = 1e-5\nsettings.eps_rel = 1e-5"
},

{
    "location": "guide/#Assembling-the-model-1",
    "page": "User Guide",
    "title": "Assembling the model",
    "category": "section",
    "text": "Once the objective function and an array of constraints have been defined, we can assemble the model withCOSMO.assemble!(model, P, q, constraints)This simply sets the corresponding variables in the model and transforms the array of constraints into the problem format defined at the top of the page.If you want to change the default settings, you can pass your settings object custom_settings to the assemble! function:COSMO.assemble!(model, P, q, constraints, settings = custom_settings)"
},

{
    "location": "guide/#Warm-starting-1",
    "page": "User Guide",
    "title": "Warm starting",
    "category": "section",
    "text": "One of the advantages of ADMM-based solvers is that they can be easily warm started. By providing starting values for the primal variable x and/or the dual variable y in the vicinity of their optimal values, the number of iterations to convergence can often be dramatically decreased.Consider the case where you have a decision variable x in mathbbR^3 and a dual variable y in mathbbR^2. Assume you expect their optimal values to be close to x_0 = (1 5 3) and y_0 = (1 2). You can pass these values when assembling the model.x_0 = [1.0; 5.0; 3.0]\ny_0 = [1.0; 2.0]\nCOSMO.assemble!(model, P, q, constraints, x0 = x_0, y0 = y_0)Another option is to useCOSMO.assemble!(model, P, q, constraints)\nwarm_start_primal!(model, x_0)\nwarm_start_dual!(model, y_0)"
},

{
    "location": "guide/#Solving-1",
    "page": "User Guide",
    "title": "Solving",
    "category": "section",
    "text": "After the model has been assembled, we can solve the problem by typingresults = COSMO.optimize!(model)Once the solver algorithm terminates, it will return a Results object that gives information about the status of the solver. If successful, it contains the optimal objective value and optimal primal and dual variables. For more information see the following section."
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
    "text": "COSMO will return one of the following statuses:Status Code Description\n:Solved An optimal solution was found\n:Unsolved Default value\n:Max_iter_reached Solver reached iteration limit (set with Settings.max_iter)\n:Time_limit_reached Solver reached time limit (set with Settings.time_limit)\n:Primal_infeasible Problem is primal infeasible\n:Dual_infeasible Problem is dual infeasible"
},

{
    "location": "guide/#Timings-1",
    "page": "User Guide",
    "title": "Timings",
    "category": "section",
    "text": "If settings.verbose_timing is set to true, COSMO will report the following times in result.times:Time Name Description\nsolver_time Total time used to solve the problem\nsetup_time Setup time = graph_time + factor_time\ngraph_time Time used to perform chordal decomposition\nfactor_time Time used to factor the system of linear equations\niter_time Time spent in iteration loop\nproj_time Time spent in projection functions\npost_time Time used for post processingIt holds: solver_time = setup_time+ iter_time + post_time,setup_time = graph_time+ factor_time,proj_time is a subset of iter_time."
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
    "text": "Our JuMP interface allows you to describe and modify your optimisation problem with JuMP and use COSMO as the backend solver. The interface is defined in /src/MOIWrapper.jl.note: Note\nCOSMO requires the newest JuMP v0.19 release that is based on the MathOptInterface package."
},

{
    "location": "jump/#Use-COSMO-1",
    "page": "JuMP Interface",
    "title": "Use COSMO",
    "category": "section",
    "text": "To specify COSMO as the solver for your JuMP model, load the solver module with using COSMO and then use the with_optimizer() function when initialising the JuMP model:m = JuMP.Model(with_optimizer(COSMO.Optimizer);"
},

{
    "location": "jump/#Specify-Solver-Settings-1",
    "page": "JuMP Interface",
    "title": "Specify Solver Settings",
    "category": "section",
    "text": "Solver-specific settings can be passed after the COSMO.Optimizer object. For example, if you want to adjust the maximum number of iterations and turn on verbose printing usem = JuMP.Model(with_optimizer(COSMO.Optimizer, max_iter = 5000, verbose = true);The full list of available settings can be found in the Settings section."
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
    "text": "This section describes COSMO\'s underlying ADMM algorithm and how the user can use the settings to adjust this algorithm. For a more detailed explanation take a look at the associated publication in Citing COSMO.COSMO solves problems with quadratic objective function and a number of conic constraints in the following form:beginarrayll mboxminimize  textstylefrac12x^top Px + q^top x mboxsubject to  Ax + s = b   s in mathcalK endarraywith primal decision variable x in mathbbR^n, primal slack variable s in mathbbR^m. The objective function is defined by positive semidefinite matrix P=P^top succeq 0 and vector q in mathbbR^n. The constraints are defined by matrix A in mathbbR^m times n, vector b in mathbbR^m and a non-empty, closed convex set mathcalK. The convex set itself can be a Cartesian product of convex sets in the form:  mathcalK = mathcalK_1^m_1 times mathcalK_2^m_2 times cdots times mathcalK_N^m_NAccordingly, by an appropriate choice of convex sets one can represent any LP, QP, SOCP or SDP."
},

{
    "location": "method/#Dual-problem-1",
    "page": "Method",
    "title": "Dual problem",
    "category": "section",
    "text": "The dual problem of the optimisation problem above is given by:beginalignat*2\ntextmaximize    quad    -textstylefrac12x^top Px - b^top  y - textsup_s in mathcalK(-y^top s ) \ntextsubject to   Px + A^top  y = -q\n              y in (mathcalK^infty)^* \nendalignat*with dual variable y in mathbbR^m."
},

{
    "location": "method/#Algorithm-1",
    "page": "Method",
    "title": "Algorithm",
    "category": "section",
    "text": "The algorithm considers a slightly transformed problem. By introducing two dummy variables tildex = x and tildes = s we can rewrite the original problem:beginalignat*2\ntextminimize         textstylefrac12tildex^top P tildex + q^top tildex + mathcalI_Ax+s=b(tildextildes) + mathcalI_mathcalK(s)\ntextsubject to   (tildextildes) = (xs)\nendalignat*where indicator functions mathcalI were used to move the constraints into the objective function. The resulting problem is now in the right format to apply the alternating direction method of multipliers (ADMM). To apply ADMM we first find the augmented Lagrangian L:  L(xstildextildeslambday) = textstylefrac12tildex^top Ptildex + q^top tildex + mathcalI_Ax+s=b(tildextildes) + mathcalI_mathcalK(s) + fracsigma2 tildex - x + textstylefrac1sigma λ _2^2 + fracrho2  tildes - s + textstylefrac1rho y _2^2Minimizing the Lagrangian in an alternating fashion with respect to the two variable pairs (tildextildes) and (xs) yields the following algorithm steps:beginalign*\n    ( tildex^k+1tildes^k+1)  rightarrow undersettildextildestextargmin   Lleft( tildextildesx^ks^ky^k right)\n    x^k+1 leftarrow tildex^k+1  labeleqADMM1\n    s^k+1 leftarrow undersetstextargmin fracrho2   tildes^k+1  - s+textstylefrac1rhoy^k _2^2 + I_mathcalK(s) \n    y^k+1 leftarrow y^k + rho left( tildes^k+1 -s^k+1 right)\nendalign*By the construction of the ADMM method those iterates are converging to the global solution. These steps are executed in a loop until convergence. Two important parameters are the ADMM steps sizes rho (Settings.rho) and sigma (Settings.sigma) which can be adjusted via the solver settings.The two most important steps of the algorithm happen in the first and third line. The evaluation of the first line turns out to be an equality constrained quadratic program. We get a solution for ( tildex^k+1tildes^k+1) at every iteration by solving the following linear system:beginalign*\nbeginbmatrix\nP + sigma I  A^top A - frac1rhoI\n    endbmatrixbeginbmatrixtildex^k+1  nu^k+1endbmatrix= beginbmatrix-q+sigma x^k b-s^k+frac1rhoy^kendbmatrix\ntildes^k+1 = s^k - frac1rholeft(nu^k+1 + y^kright)\nendalign*Fortunately, the left hand matrix doesn\'t change, which is why COSMO only has to factor the matrix once at the start of the algorithm.The second important step in the algorithm is the update equation for s^k+1 which can be interpreted as a projection onto the constraint set mathcalK:s^k+1 = Pi_mathcalKleft( tildes^k+1 + frac1rhoy^kright)The computational cost of this projection is highly dependent on the constraints of the problem. While projections onto the zero set or the nonnegative orthant are inexpensive, projections onto the positive semidefinite cone of order N involve an eigen-decomposition. Since methods for eigen-decompositions have a complexity of mathcalO(N^3) the projection can become the computationally most expensive operation of the algorithm."
},

{
    "location": "method/#Scaling-1",
    "page": "Method",
    "title": "Scaling",
    "category": "section",
    "text": "The convergence of ADMM-based algorithms depends on the relative scaling of the problem data. Especially to improve the convergence of badly scaled problems, COSMO tries to rescale the data in a preprocessing step.We rescale the equality constraints with diagonal positive semidefinite matrices D and E. The scaled problem is given by:beginalignat*2\nlabeleqscaled\ntextminimize    quad    textstylefrac12 hatx^top hatP hatx + hatq^top hatx\ntextsubject to   hatA hatx + hats  = hatb  \n              hats in EmathcalK\nendalignat*with scaled problem data  hatP=DPD quad hatq=Dq  quadhatA=EAD quad hatb=Eband the scaled convex cone EmathcalK = Ev in mathbbR^m mid v in mathcalK . After solving the scaled problem the original solution is obtained by reversing the scaling:   x = Dhatx quad s = E^-1hats quad y = EhatyTo obtain the scaling matrices D and E we use a modified Ruiz equilibration algorithm which involves a certain number of scaling iterations to equilibrate the column norms of the data matrices P and A. The number of these iterations can be adjusted by the user with scaling. To disable the scaling step set scaling = 0."
},

{
    "location": "method/#Termination-criteria-1",
    "page": "Method",
    "title": "Termination criteria",
    "category": "section",
    "text": "The COSMO algorithm can terminate for four reasons:The maximum number of allowed iterations has been reached. The user can specify this value in the solver settings with max_iter.\nThe solver runtime reaches the time limit specified by the user (time_limit).\nCOSMO detects an infeasible problem.\nThe iterates fulfil the termination criteria for convergence.COSMO uses the primal and dual residuals of the problem to determine if the algorithm has converged. The primal and dual residuals are given by:beginalign*\nr_p = Ax + s -b\nr_d = Px + q + A^top y\nendalign*The solver terminates when the infty-norms of the residuals lie below a specified tolerance. COSMO uses the sum of an absolute and relative tolerance term:beginalign*\n  r_p^k _infty leq epsilon_mathrmabs + epsilon_mathrmrel  textmax left Ax^k _inftys^k_infty b_infty right\n   r_d^k_infty leq epsilon_mathrmabs + epsilon_mathrmrel  textmax leftPx^k_inftyq_infty A^top y^k_infty right\nendalign*The absolute and relative tolerances epsilon_mathrmabsand epsilon_mathrmrel can be set by the user by specifying eps_abs and eps_rel. Furthermore, the user can adjust the number of iterations after which the convergence criteria are checked (check_termination)."
},

{
    "location": "method/#Infeasibility-detection-1",
    "page": "Method",
    "title": "Infeasibility detection",
    "category": "section",
    "text": "COSMO uses conditions based on separating hyperplanes to detect infeasible problems. The conditions for COSMO\'s problem format have been developed in [1]. Define the convex set mathcalC = mathcal-K + b then we can use the following infeasibility conditions:beginalign*\nmathcalP = leftx in mathbbR^n mid  Px = 0  Ax in mathcalC^infty  langle  qx rangle  0  right \nmathcalD = lefty in mathbbR^m mid  A^top  y  = 0   sigma_mathcalC(y)  0 right\nendalign*The existence of some y in mathcalD is a certificate that the problem is primal infeasible, while the existence of some x in mathcalP is a certificate for dual infeasibility. COSMO regularly checks above conditions to detect infeasible problems. If the detection is successful, the solver terminates and returns the status codes :Primal_infeasible or :Dual_infeasible. COSMO checks the conditions every check_infeasibility iterations, which can be adjusted by the user."
},

{
    "location": "method/#References-1",
    "page": "Method",
    "title": "References",
    "category": "section",
    "text": "[1] Banjac, G. et al. Infeasibility detection in the alternating direction method of multipliers for convex optimization. Preprint, 2017."
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
    "text": "We want to solve the following linear program with decision variable x:beginarrayll mboxminimize   c^top x\nmboxsubject to   A x leq b \n                    x geq 1 \n                    x_2 geq 5 \n                    x_1 + x_3 geq 4\nendarrayThe problem can be solved with COSMO in the following way:using COSMO, LinearAlgebra, SparseArrays, Test\n\nc = [1; 2; 3; 4.]\nA = Matrix(1.0I, 4, 4)\nb = [10; 10; 10; 10]\nn = 4\n# -------------------\n# create constraints A * x + b in set\n# -------------------\n# Ax <= b\nc1 = COSMO.Constraint(-A, b, COSMO.Nonnegatives)\n# x >= 1\nc2 = COSMO.Constraint(Matrix(1.0I, n, n), -ones(n), COSMO.Nonnegatives)\n# x2 >= 5\nc3 = COSMO.Constraint(1, -5, COSMO.Nonnegatives, n, 2:2)\n# x1 + x3 >= 4\nc4 = COSMO.Constraint([1 0 1 0], -4, COSMO.Nonnegatives)\n\n# -------------------\n# define cost function\n# -------------------\nP = spzeros(4, 4)\nq = c\n\n# -------------------\n# assemble solver model\n# -------------------\nsettings = COSMO.Settings(max_iter=2500, verbose=true, eps_abs = 1e-4, eps_rel = 1e-5)\nmodel = COSMO.Model()\nassemble!(model, P, q, [c1; c2; c3; c4], settings = settings)\nres = COSMO.optimize!(model);\n\n@testset \"Linear Problem\" begin\n  @test isapprox(res.x[1:4], [3; 5; 1; 1], atol=1e-2, norm = (x -> norm(x, Inf)))\n  @test isapprox(res.obj_val, 20.0, atol=1e-2)\nend"
},

{
    "location": "examples/#Closest-Correlation-Matrix-1",
    "page": "Examples",
    "title": "Closest Correlation Matrix",
    "category": "section",
    "text": "We consider the problem of finding the closest correlation matrix X to a given random matrix C. With closest correlation matrix we mean a positive semidefinite matrix with ones on the diagonal. The problem is given by:beginarrayll mboxminimize   frac12X - C_F^2\nmboxsubject to   X_ii = 1 quad i=1dotsn \n                    X succeq 0\nendarrayNotice how JuMP is used to describe the problem. COSMO is chosen as the backend solver using JuMP\'s with_optimizer() function.using COSMO, JuMP, LinearAlgebra, SparseArrays, Test, Random\nrng = Random.MersenneTwister(12345);\n\n# create a random test matrix C\nn = 8\nC = -1 .+ rand(rng, n, n) .* 2;\nc = vec(C);\n\n# define problem in JuMP\nq = -vec(C);\nr = 0.5 * vec(C)\' * vec(C);\nm = JuMP.Model(with_optimizer(COSMO.Optimizer, verbose=true, eps_abs = 1e-4));\n@variable(m, X[1:n, 1:n], PSD);\nx = vec(X);\n@objective(m, Min, 0.5 * x\' * x  + q\' * x + r)\nfor i = 1:n\n  @constraint(m, X[i, i] == 1.)\nend\n\n# solve and get results\nstatus = JuMP.optimize!(m)\nobj_val = JuMP.objective_value(m)\nX_sol = JuMP.value.(X)"
},

{
    "location": "citing/#",
    "page": "Citing COSMO",
    "title": "Citing COSMO",
    "category": "page",
    "text": ""
},

{
    "location": "citing/#Citing-COSMO-1",
    "page": "Citing COSMO",
    "title": "Citing COSMO",
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
    "location": "api/#COSMO.Model",
    "page": "API Reference",
    "title": "COSMO.Model",
    "category": "type",
    "text": "Model()\n\nInitializes an empty COSMO model that can be filled with problem data using assemble!(model, P, q,constraints; [settings, x0, s0, y0]).\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.assemble!",
    "page": "API Reference",
    "title": "COSMO.assemble!",
    "category": "function",
    "text": "assemble!(model, P, q, constraint(s); [settings, x0, y0, s0])\n\nAssembles a COSMO.Model with a cost function defind by P and q, and a number of constraints.\n\nThe positive semidefinite matrix P and vector q are used to specify the cost function of the optimization problem:\n\nmin   1/2 x\'Px + q\'x\ns.t.  Ax + b ∈ C\n\nconstraints is a COSMO.Constraint or an array of COSMO.Constraint objects that are used to describe the constraints on x.\n\n\n\nThe optional keyword argument settings can be used to pass custom solver settings:\n\ncustom_settings = COSMO.Settings(verbose = true);\nassemble!(model, P, q, constraints, settings = custom_settings)\n\n\n\nThe optional keyword arguments x0, s0, and y0 can be used to provide the solver with warm starting values for the primal variable x, the primal slack variable s and the dual variable y.\n\nx_0 = [1.0; 5.0; 3.0]\nCOSMO.assemble!(model, P, q, constraints, x0 = x_0)\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.set!",
    "page": "API Reference",
    "title": "COSMO.set!",
    "category": "function",
    "text": "set!(model, P, q, A, b, convex_sets, [settings])\n\nSets model data directly based on provided fields.\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.empty_model!",
    "page": "API Reference",
    "title": "COSMO.empty_model!",
    "category": "function",
    "text": "empty_model!(model)\n\nResets all the fields of model to that of a model created with COSMO.Model() (apart from the settings).\n\n\n\n\n\n"
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
    "text": "COSMO.Model\nCOSMO.assemble!\nCOSMO.set!\nCOSMO.empty_model!\nCOSMO.warm_start_primal!\nCOSMO.warm_start_slack!\nCOSMO.warm_start_dual!"
},

{
    "location": "api/#COSMO.Constraint",
    "page": "API Reference",
    "title": "COSMO.Constraint",
    "category": "type",
    "text": "Constraint(A, b, convex_set, dim = 0, indices = 0:0)\n\nCreates a COSMO constraint: Ax + b ∈ convex_set.\n\nBy default the following convex sets are supported: ZeroSet, Nonnegatives, SecondOrderCone, PsdCone, PsdConeTriangle.\n\nExamples\n\njulia> Constraint([1 0;0 1], zeros(2), COSMO.Nonnegatives)\nConstraint\nSize of A: (2, 2)\nConvexSet: Nonnegatives{Float64}\n\n\n\nThe optional arguments dim and indices can be used to specify A and b for subparts of variable x. If x has dimension dim = 4, then x[2] and x[3] can be constrained to the zero cone in the following way:\n\nExamples\n\njulia> c = Constraint([1 0;0 1], zeros(2), COSMO.ZeroSet, 4, 2:3)\nConstraint\nSize of A: (2, 4)\nConvexSet: ZeroSet{Float64}\n\nNotice that extra columns of A have been added automatically.\n\njulia>Matrix(c.A)\n2×4 Array{Float64,2}:\n0.0  1.0  0.0  0.0\n0.0  0.0  1.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.ZeroSet",
    "page": "API Reference",
    "title": "COSMO.ZeroSet",
    "category": "type",
    "text": "ZeroSet(dim)\n\nCreates the zero set  0 ^dim of dimension dim. If x ∈ ZeroSet then all entries of x are zero.\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.Nonnegatives",
    "page": "API Reference",
    "title": "COSMO.Nonnegatives",
    "category": "type",
    "text": "Nonnegatives(dim)\n\nCreates the nonnegative orthant  x in mathbbR^dim  x ge 0   of dimension dim.\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.SecondOrderCone",
    "page": "API Reference",
    "title": "COSMO.SecondOrderCone",
    "category": "type",
    "text": "SecondOrderCone(dim)\n\nCreates the second-order cone (or Lorenz cone)  (tx) in mathrmR^dim   x _2  leq t .\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.PsdCone",
    "page": "API Reference",
    "title": "COSMO.PsdCone",
    "category": "type",
    "text": "PsdCone(dim)\n\nCreates the cone of symmetric positive semidefinite matrices mathcalS_+^dim. The entries of the matrix X are stored column-by-column in the vector x of dimension dim. Accordingly  X in mathbbS_+ Rightarrow x in mathcalS_+^dim, where X = textmat(x).\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.PsdConeTriangle",
    "page": "API Reference",
    "title": "COSMO.PsdConeTriangle",
    "category": "type",
    "text": "PsdConeTriangle(dim)\n\nCreates the cone of symmetric positive semidefinite matrices. The entries of the upper-triangular part of matrix X are stored in the vector x of dimension dim. A r times r matrix has r(r+1)2 upper triangular elements and results in a vector of mathrmdim = r(r+1)2.\n\nExamples\n\nThe matrix\n\nbeginbmatrix x_1  x_2  x_4 x_2  x_3  x_5 x_4  x_5  x_6 endbmatrix\n\nis transformed to the vector x_1 x_2 x_3 x_4 x_5 x_6^top with corresponding constraint  PsdConeTriangle(6).\n\n\n\n\n\n"
},

{
    "location": "api/#Constraints-1",
    "page": "API Reference",
    "title": "Constraints",
    "category": "section",
    "text": "COSMO.Constraint\nCOSMO.ZeroSet\nCOSMO.Nonnegatives\nCOSMO.SecondOrderCone\nCOSMO.PsdCone\nCOSMO.PsdConeTriangle"
},

{
    "location": "api/#COSMO.Settings",
    "page": "API Reference",
    "title": "COSMO.Settings",
    "category": "type",
    "text": "Settings(;arg=val)\n\nCreates a COSMO settings object that is used to pass user settings to the solver.\n\nArgument Description Values (default)\nrho ADMM rho step 0.1\nsigma ADMM sigma step 1e-6\nalpha Relaxation parameter 1.6\neps_abs Absolute residual tolerance 1e-4\neps_rel Relative residual tolerance 1e-4\neps_prim_inf Primal infeasibility tolerance 1e-4\neps_dual_inf Dual infeasibility tolerance 1e-4\nmax_iter Maximum number of iterations 2500\nverbose Verbose printing false\nverbose_timing Verbose timing false\ncheck_termination Check termination interval 40\ncheck_infeasibility Check infeasibility interval 40\nscaling Number of scaling iterations 10\nadaptive_rho Automatic adaptation of step size parameter true\ntime_limit set solver time limit in s 0\n\n\n\n\n\n"
},

{
    "location": "api/#Settings-1",
    "page": "API Reference",
    "title": "Settings",
    "category": "section",
    "text": "COSMO.Settings"
},

{
    "location": "api/#COSMO.Result",
    "page": "API Reference",
    "title": "COSMO.Result",
    "category": "type",
    "text": "Result{T <: AbstractFloat}\n\nObject returned by the COSMO solver after calling optimize!(model). It has the following fields:\n\nFieldname Type Description\nx Vector{T} Primal variable\ny Vector{T} Dual variable\ns Vector{T} (Primal) set variable\nobj_val T Objective value\niter Int64 Number of iterations\nstatus Symbol Solution status\ninfo COSMO.ResultInfo Struct with more information\ntimes COSMO.ResultTimes Struct with several measured times\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.ResultTimes",
    "page": "API Reference",
    "title": "COSMO.ResultTimes",
    "category": "type",
    "text": "ResultTimes{T <: AbstractFloat}\n\nPart of the Result object returned by the solver. ResultTimes contains timing results for certain parts of the algorithm:\n\nTime Name Description\nsolver_time Total time used to solve the problem\nsetup_time Setup time = graphtime + factortime\ngraph_time Time used to perform chordal decomposition\nfactor_time Time used to factor the system of linear equations\niter_time Time spent in iteration loop\nproj_time Time spent in projection functions\npost_time Time used for post processing\n\nBy default COSMO only measures solver_time, setup_time and proj_time. To measure the other times set verbose_timing = true.\n\n\n\n\n\n"
},

{
    "location": "api/#COSMO.ResultInfo",
    "page": "API Reference",
    "title": "COSMO.ResultInfo",
    "category": "type",
    "text": "ResultInfo{T <: AbstractFloat}\n\nObject that contains further information about the primal and dual residuals.\n\n\n\n\n\n"
},

{
    "location": "api/#Results-1",
    "page": "API Reference",
    "title": "Results",
    "category": "section",
    "text": "COSMO.Result\nCOSMO.ResultTimes\nCOSMO.ResultInfo"
},

]}
