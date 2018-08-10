# User facing functions / structs:

# ------------------------------
# Structs: (defined in Types.jl)
# ------------------------------
# 1. QOCS.Model
# 2. QOCS.Settings
# 3. QOCS.Results
# 4. QOCS.ConvexSet (Abstract)
# 4a. Subtypes of ConvexSet: ZeroCone, NonNegativeOrthant, Interval (Box), SecondOrderCone, PositiveSemidefiniteCone <: QOCS.ConvexSet, User defined ones
# Each Subtype of ConvexSet needs an associated projection function project!(set::Type{<:ConvexSet}). We'll provide the standard ones

# ------------------------------
# Functions: (defined in Interface.jl)
# ------------------------------

# 1. assemble!() --> To copy problem data and warm start data into the model
# Idea: two assemble! functions:
# 1a. Low level: Pass the problem data exactly in the current problem format(A,q,A,b,K) with K sedumi style
# 1b. High Level: Pass P,q directly, pass constraints as an array of ConvexSet Subtypes

# 2. optimize!(model::QOCS.Model,settings::QOCS.settings) --> Calls the solver based on the model data and returns a QOCS.Results object
# 3. reset!() --> Empty functions to reset each struct QOCS.Model, QOCS.Settings, QOCS.Results objects

# ------------------------------
# Open questions:
# ------------------------------
# 1. How best handle warm starting?
# 2. How best handle update of problem data without re-assambling the whole model?
# 3. Ideally work directly on the QOCS.Model data without copying anything (that means that matrixes will be changed through scaling) --> Problem?


