#=
This file contains code for separating a timeseries into its seasonal and residual
fluctuations. It interfaces SignalSeparation.jl and also uses some analytic knowledge
of the form of the incoming solar radiation
=#
using SignalDecomposition

export decompose, Fourier, MeanAlbedo, ProductInversion, fielddecompose, Sinusoidal
export nrmse, rmse, TimeAnomaly

# my own decomposition method that dispatches on dimensional array for input method

function fielddecompose(A::AbDimArray, method::Decomposition)
    t = dims(A, Time)
    fielddecompose(t, A, method)
end
function fielddecompose(A::AbDimArray, method::Fourier)
    t = 1:size(A, Time)
    fielddecompose(t, A, method)
end

function fielddecompose(t::AbstractVector, A::AbDimArray, method::Decomposition)
    seasonal = DimensionalArray(copy(Array(A)), dims(A), A.name*"x")
    residual = DimensionalArray(copy(Array(A)), dims(A), A.name*"r")
    for i in spatialidxs(A)
        y = Array(A[i...])
        sea, res = SignalDecomposition.decompose(t, y, method)
        seasonal[i...] .= sea
        residual[i...] .= res
    end
    return seasonal, residual
end

SignalDecomposition.decompose(t::Vector{<:TimeType}, s, m::Fourier) =
decompose(1:length(t), s, m)

SignalDecomposition.decompose(s::AbDimArray, m::Sinusoidal) =
decompose(dims(s, Time).val, Array(s), m)

function SignalDecomposition.decompose(t::AbstractVector{<:TimeType}, s, m::Sinusoidal)
    truetime = time_in_days(t)
    decompose(truetime, s, m)
end

#########################################################################
# Mean Albedo
#########################################################################
"""
    MeanAlbedo(S::Vector) <: Decomposition
    MeanAlbedo(t::Vector, φ::Real)
Decomposition that splits a timeseries `R` into its mean albedo component
``\\bar{\\alpha}S`` with ``S`` the Solar
insolation (either given directly as `S` or calculated analytically with given time
`t` and latitude `φ`) and ``\\bar{\\alpha} = `` `timeaverage(R)/timeaverage(S)`,
and the residuals of that.

It is supposed to approximately purely internal dynamical fluctuations.
"""
struct MeanAlbedo{T<:Real} <: Decomposition
    solar::Vector{T}
end
MeanAlbedo(t::Vector, φ::Real) = MeanAlbedo(insolation.(t, φ))
MeanAlbedo(x) = MeanAlbedo(Array(x))

# Extention for convenient usage in decompose below
ClimateBase.timeagg(f, t::Vector, a::AbDimArray{<:Any, 1}) = timeagg(f, t, a.data)

function SignalDecomposition.decompose(t, x, method::MeanAlbedo)
    @assert length(x) == length(method.solar)
    if eltype(t) <: TimeType
        meanα = timeagg(mean, t, x)/timeagg(mean, t, method.solar) # average albedo
    else
        meanα = mean(x)/mean(method.solar)
    end
    seasonal = meanα .* method.solar
    return seasonal, x .- seasonal
end

function meanalbedo_decompose(A::AbDimArray, S)
    solar = DimensionalArray(copy(Array(A)), dims(A), "MeanAlbedo αS")
    residual = DimensionalArray(copy(Array(A)), dims(A), "MeanAlbedo r")
    for i in spatialidxs(A)
        y = Array(A[i...])
        s = Array(S[i...])
        method = MeanAlbedo(s)
        sea, res = SignalDecomposition.decompose(y, method)
        solar[i...] .= sea
        residual[i...] .= res
    end
    return solar, residual
end
