
# -------------------------------------
# Convex Set -  Definitions
# -------------------------------------
abstract type AbstractConvexSet end


mutable struct Zeros <:AbstractConvexSet
  dim::Int
  indices::UnitRange{Int64}
  function Zeros()
    return new(0,0:0)
  end

end

mutable struct Nonnegatives <:AbstractConvexSet
  dim::Int
  indices::UnitRange{Int64}
  function Nonnegatives()
    return new(0,0:0)
  end

end


mutable struct Box <: AbstractConvexSet
    dim::Int
    indices::UnitRange{Int64}
    l::AbstractVector{<:Real}
    u::AbstractVector{<:Real}

    function Box(l::AbstractVector{<:Real},u::AbstractVector{<:Real})
      # FIXME: For now intervals are not handled directly, but as one projection on nonnegative orthant
      return new(0,0:0,l,u)
  end
end


mutable struct SecondOrderCone <:AbstractConvexSet
  dim::Int
  indices::UnitRange{Int64}
   function SecondOrderCone()
    return new(0,0:0)
  end

end

mutable struct PositiveSemidefiniteCone <:AbstractConvexSet
  dim::Int
  indices::UnitRange{Int64}

   function PositiveSemidefiniteCone()
    return new(0,0:0)
  end
end



# -------------------------------------
# IsInDualCone -  Functions
# -------------------------------------

function inDual(x::SubArray{Float64},convexSet::Zeros,tol::Float64)
  return true
end

# Note: R+ is self-dual
function inDual(x::SubArray{Float64},convexSet::Nonnegatives,tol::Float64)
  for i in eachindex(x)
    if x[i] <= -tol
      return false
    end
  end
  return true
end

function inDual(x::SubArray{Float64},convexSet::Box,tol::Float64)
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
function inDual(x::SubArray{Float64},convexSet::SecondOrderCone,tol::Float64)
  if norm(x[2:end],2) - x[1] <= tol
    return true
  else
    return false
  end
end

# Note: Positive semidefinite cone is self-dual
function inDual(x::SubArray{Float64},convexSet::PositiveSemidefiniteCone,tol::Float64)
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
function inRecc(x::SubArray{Float64},convexSet::Zeros,tol::Float64)
   for i in eachindex(x)
    if abs(x[i]) > tol
      return false
    end
  end
  return true
end

function inRecc(x::SubArray{Float64},convexSet::Nonnegatives,tol::Float64)
   for i in eachindex(x)
    if x[i] > tol
      return false
    end
  end
  return true
end

function inRecc(x::SubArray{Float64},convexSet::Box,tol::Float64)
  #  for i in eachindex(x)
  #   if x[i] > tol
  #     return false
  #   end
  # end
  return true
end

function inRecc(x::SubArray{Float64},convexSet::SecondOrderCone,tol::Float64)
  if norm(x[2:end],2) + x[1] <= tol
    return true
  else
    return false
  end
end

function inRecc(x::SubArray{Float64},convexSet::PositiveSemidefiniteCone,tol::Float64)
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

