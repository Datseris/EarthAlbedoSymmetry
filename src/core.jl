#=
Core processing file for this project containing general context agnostic tools:
* data conversion
* aggregation
* convenience tools for DateTime, e.g. ensuring yearly averages
=#

# TODO: Add examples to all docstrings

shortname(attrib) = attrib["CF_name"]*" [$(attrib["units"])]"
longname(attrib) = attrib["long_name"]*" [$(attrib["units"])]"

export shortname, longname

export totalaverage
totalaverage(x) = spacemean(timemean(x))

using Distances, StaticArrays
export findnearest, SVector, split_at_longitude

"""
    split_at_longitude(A, λ) → A1, A2
Split `A` at given longitude, taking into account the periodic nature of the longitude.
It is expected that input `A` covers all longitudes from 0 to 360.
"""
function split_at_longitude(A, λ)
    λi = findnearest(λ, dims(A, Lon))
    # here 179 explicitly assumes longitudinal spacing of 1 degree.
    # TODO: change it to use internal spacing, and do 180 - δλ
    i1 = mod1.(λi:λi+179, size(A, Lon))
    i2 = setdiff(1:size(A, Lon), i1)
    return A[Lon(i1)], A[Lon(i2)]
end

"""
    findnearest(val, A)
Return the index of `A` which has value nearest to `val`.
"""
function findnearest(val, A)
    i = 1
    d = evaluate(Euclidean(), val, A[i])
    @inbounds for j in 1:length(A)
        dd = evaluate(Euclidean(), val, A[j])
        if dd < d
            i = j
            d = dd
        end
    end
    return i
end

export timeseries

_correct(::Coord, r) = r
_correct(::Coord, r::Real) = SVector(0, r)
_correct(::Lat, r) = r[2]
_correct(::Lat, r::Real) = r

"""
    timeseries(F, g)
Return the timeseries of field `F` at point `g` (either latitude or (lon, lat)
vector). Works with both global `F` as well as zonally-averaged.
"""
function timeseries(F, r)
    if hasdim(F, Lon)
        p = _correct(Coord(), r)
        i = findnearest(p[1], dims(F, 1))
        j = findnearest(p[2], dims(F, 2))
        given = SVector(dims(F, 1)[i], dims(F, 2)[j])
        a = Array(F[i, j, :])
    else
        p = _correct(dims(F, 1), r) #assumes spatial dim is first dim
        i = findnearest(p, dims(F, 1))
        given = dims(F, 1)[i]
        a = Array(F[i, :])
    end
    return a, given
end


using Dates: TimeType

"""
    monthly_insolation(t::TimeType, args...)
Average `insolation(τ, ...)` for `τ` in the month of given time.
"""
function monthly_insolation(t::TimeType, args...)
    d = monthspan(t)
    mean(insolation(τ, args...) for τ in d)
end

export linreg
"""
    linreg(x, y) -> s, o, fit, x
Return best linear fit so that y ≈ s*x + o.
If `x` is a date vector, `1:length(x)` is used instead.
"""
linreg(x::AbstractArray{<:Dates.AbstractDateTime}, y) = linreg(1:length(x), y)
function linreg(x, y)
    X = ones(length(x), 2)
    X[:, 1] .= x
    c = X\y
    s, o = c
    f = s .* x .+ o
    return s, o, f, x
end
