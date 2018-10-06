
# -------------------------------------
# Convex Set -  Definitions
# -------------------------------------
abstract type AbstractConvexSet end


mutable struct Zeros <:AbstractConvexSet
  dim::Int
  indices::UnitRange{Int64}
  project!::Function
  scale!::Function
  inDual::Function
  inRecc::Function
  function Zeros()
    return new(0,0:0,QOCS.Projections.zeroCone!,QOCS.scale!,QOCS.inFreeCone,QOCS.inPolRecConeZeros)
  end

end

mutable struct Nonnegatives <:AbstractConvexSet
  dim::Int
  indices::UnitRange{Int64}
  project!::Function
  scale!::Function
  inDual::Function
  inRecc::Function
  function Nonnegatives()
    return new(0,0:0,QOCS.Projections.nonNegativeOrthant!,QOCS.scale!,QOCS.inNonnegatives,QOCS.inPolRecConeNonnegatives)
  end

end


mutable struct Box <: AbstractConvexSet
    dim::Int
    indices::UnitRange{Int64}
    project!::Function
    scale!::Function
    inDual::Function
    inRecc::Function
    l::AbstractVector{<:Real}
    u::AbstractVector{<:Real}

    function Box(l::AbstractVector{<:Real},u::AbstractVector{<:Real})
      # FIXME: For now intervals are not handled directly, but as one projection on nonnegative orthant
      return new(0,0:0,QOCS.Projections.box!,QOCS.scale!,QOCS.inDualBox,QOCS.inPolRecConeBox,l,u)
  end
end


mutable struct SecondOrderCone <:AbstractConvexSet
  dim::Int
  indices::UnitRange{Int64}
  project!::Function
  scale!::Function
  inDual::Function
  inRecc::Function

   function SecondOrderCone()
    return new(0,0:0,QOCS.Projections.secondOrderCone!,QOCS.scale!,QOCS.inSOC,QOCS.inPolRecConeSOC)
  end

end

mutable struct PositiveSemidefiniteCone <:AbstractConvexSet
  dim::Int
  indices::UnitRange{Int64}
  project!::Function
  scale!::Function
  inDual::Function
  inRecc::Function
  # Block-Lanczos data
  subspace::Array{Float64,2}
  subspace_dimension::Int
  positive_subspace::Bool

  function PositiveSemidefiniteCone()
    return new(0,0:0,QOCS.Projections.sdcone!,QOCS.scale!,QOCS.inPSD,QOCS.inPolRecPSD,zeros(0, 0),-1,true)
  end
end



# -------------------------------------
# IsInDualCone -  Functions
# -------------------------------------

function inFreeCone(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
  return true
end

# Note: R+ is self-dual
function inNonnegatives(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
  for i in eachindex(x)
    if x[i] <= -tol
      return false
    end
  end
  return true
end

function inDualBox(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
  l = convexSet.l
  u = convexSet.u
  for i in eachindex(x)
    if x[i] >= l[i]-tol || x[i] <= u[i]+tol
      return false
    end
  end
  return true
end

# Note: Second-order-cone is self-dual
function inSOC(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
  if norm(x[2:end],2) - x[1] <= tol
    return true
  else
    return false
  end
end

# Note: Positive semidefinite cone is self-dual
function inPSD(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
  dim = Int(sqrt(length(x)))
  if minimum(real(eigen(reshape(Array(x),dim,dim)).values)) >= -tol
    return true
  else
    return false
  end
end

# -------------------------------------
# IsInPolarRecessionCone -  Functions
# -------------------------------------
function inPolRecConeZeros(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
   for i in eachindex(x)
    if abs(x[i]) > tol
      return false
    end
  end
  return true
end

function inPolRecConeNonnegatives(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
   for i in eachindex(x)
    if x[i] > tol
      return false
    end
  end
  return true
end

function inPolRecConeBox(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
  #  for i in eachindex(x)
  #   if x[i] > tol
  #     return false
  #   end
  # end
  return true
end

function inPolRecConeSOC(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
  if norm(x[2:end],2) + x[1] <= tol
    return true
  else
    return false
  end
end

function inPolRecPSD(x::SubArray{Float64},convexSet::AbstractConvexSet,tol::Float64)
  dim = Int(sqrt(length(x)))
  if maximum(real(eigen(reshape(Array(x),dim,dim)).values)) <= tol
    return true
  else
    return false
  end
end


scale!(set::Union{Zeros,Nonnegatives}) = false
scale!(set::Box) = false,set.l,set.u
scale!(set::Union{SecondOrderCone,PositiveSemidefiniteCone}) = true

