#=
Functions for decomposing albedo into contributions, or estimating effective
albedo fields.
=#
export surface_atmosphere_contributions, cloud_effective_albedo, fill_missings
export total_albedo

"""
    cloud_effective_albedo(f, τ, g = 0.9; frequencies = [1, 2]) → αC
Estimate an effective cloud albedo field using the formula
```math
\\alpha_C \\approx f \\frac{\\sqrt{3}(1-g)\\tau}{2 + \\sqrt{3}(1-g)\\tau}
```
with ``f`` cloud fraction, ``\\tau`` cloud optical depth and ``g`` the asymmetry factor
for the cloud particle phase function. The formula is eq. (19) of [^Lacis1974] multiplied
with the cloud fraction.

`τ` has missing values (typically), and thus a regularization procedure is performed on
it first by fitting cosines.

[^Lacis1974]: [Lacis, A. A., & Hansen, J. (1974), A Parameterization for the Absorption of Solar Radiation in the Earth’s Atmosphere](https://doi.org/10.1175/1520-0469(1974)031<0118:APFTAO>2.0.CO;2)
"""
function cloud_effective_albedo(f, τ, g = 0.9; frequencies = [1.0, 2.0])
    if any(ismissing, τ) # missing values, perform regularization
        @assert hasdim(τ, Time)
        T = sinusoidal_continuation(τ, frequencies)
    else
        T = τ
    end
    p = @. sqrt(3)*(1 - g)
    R = p .* T ./ (2 .+ p .* T)
    cloud_frac = any(x -> x > 1, f) ? f ./ 100 : f
    eca = cloud_frac .* R
    gstr = g isa Real ? string(g) : "array"
    return ClimArray(eca.data, eca.dims, "effective cloud albedo, g = $(gstr)")
end
