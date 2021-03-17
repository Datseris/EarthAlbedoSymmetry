#=
Function for basic statistically-relevant properties, e.g. mutual information
fit errors, etc.
=#
using Statistics
export meanstd, standard
meanstd(x) = mean(x), std(x)
standard(x) = (x .- mean(x)) ./ std(x)

#########################################################################
# Basic statistical properties
#########################################################################
export yearly_variability
using Dates: AbstractDateTime

"""
    yearly_variability(t, x) → dates, res
Decompose some vector `x` to its yearly variability with respect to time vector `t`,
i.e. to all the values of `x` that correspond to a unique day and month.
Return `dates`, which are the unique dates (tuples (day, month)) and `res`,
where `res[i]` are all the values of `x` that correspond to date `dates[i]`.
"""
function yearly_variability(t, xs)
    time = Array(t)
    dates = unique(daymonth.(time))
    res = [Float64[] for u in dates]
    for (t, x) in zip(time, xs)
        i = findfirst(isequal(daymonth(t)), dates)
        push!(res[i], x)
    end
    return dates, res
end

using Statistics
export periodic_correlation
"""
	periodic_correlation(x, y, lags, demean=true)
Calculate the correlation function of periodic signals both of which have period
their length (e.g. functions of longitude).
"""
function periodic_correlation(x,y,τ, demean=true)
	if length(unique(x)) == length(unique(y)) == 1
		# Dividing by `std` if `x` is constant leads to `Inf` results.
		return ones(eltype(x), length(τ))
	end
    if demean
        x = x .- mean(x)
        y = y .- mean(y)
    end
    N = length(x)
    r = [sum(x[n]*y[mod1(n+m, N)] for n in 1:N)/N for m in τ]
    r ./ (std(x)*std(y)) # follow normalization of `StatsBase.crosscor`
end
