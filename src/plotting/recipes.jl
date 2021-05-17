include("style.jl")

# Uncomment the following to install the necessary packages for plotting.
# Pkg.add(["Conda", "PyCall"])
# using Conda
# Conda.add("cartopy")

using PyCall, PyPlot, ClimateBase
ccrs = pyimport("cartopy.crs")
LONLAT = ccrs.PlateCarree()
DEFPROJ = ccrs.Mollweide()
Time = ClimateBase.Ti

"""
    earthscatter(A::ClimArray, projection = ccrs.Mollweide(); kwargs...) → fig, ax, cb
Plot a `ClimArray` with space type `GaussianEqualArea` as a scatter plot with the color
of the points being the value of `A` at these points. This requires that `A` has only
one dimension, the coordinates.

Keyword values are propagated to `scatter` and keys of interest are `cmap, vmin, vmax`.
"""
function earthscatter(A::ClimArray, projection = DEFPROJ; kwargs...)
    coords = dims(A, Coord)
    @assert length(dims(A)) == 1 "Input ClimArray must only have one `Coord` dimension."
    eqlon = [l[1] for l in coords]
    eqlat = [l[2] for l in coords]
    fig = figure()
    fig.tight_layout()
    ax = subplot(111, projection=projection)
    ax, cb = earthscatter!(ax, eqlon, eqlat, A.data, projection; kwargs...)
    if A.name ≠ Symbol("")
        ax.set_title(string(A.name))
    end
    return fig, ax, cb
end

function earthscatter!(ax, A, projection = DEFPROJ; kwargs...)
    coords = dims(A, Coord)
    lon = [l[1] for l in coords]
    lat = [l[2] for l in coords]
    earthscatter!(ax, lon, lat, A, projection; kwargs...)
end

function earthscatter!(ax, lon, lat, A, projection = DEFPROJ;
    coast = true, vmin = minimum(A), vmax = maximum(A), levels = 11,
    ticks = nothing, kwargs...)
    lon = Array(lon)
    lat = Array(lat)
    lvls = range(vmin, vmax, length = levels)
    cmap = matplotlib.cm.get_cmap(get(kwargs, :cmap, :viridis), length(lvls)-1)
    norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
    cax, kw = matplotlib.colorbar.make_axes(ax,location="right",pad=0.02,shrink=0.8)
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)

    ax.scatter(lon, lat; transform = LONLAT, c=A, kwargs..., cmap = cmap, vmin = vmin, vmax = vmax)
    coast && ax.coastlines()
    if projection == LONLAT
        ax.set_xticks(-180:60:180, crs=LONLAT)
        ax.set_yticks(-90:30:90, crs=LONLAT)
    end
    ax.gridlines(alpha = 0.25)
    cb.set_ticks(isnothing(ticks) ? lvls[1:max(1, levels÷5):end] : ticks)
    return ax, cb
end

"""
    earthsurface(A::ClimArray, projection = ccrs.Mollweide(); kwargs...) → fig, ax, cb
Plot a `ClimArray` with space type `LonLatGrid` as a surface plot with the color
of the surface being the value of `A`. This requires that `A` has exactly two dimensions,
`Lon, Lat`.

Keyword values are `coast=true` to enable coastlines and also `cmap, vmin, vmax, levels`.
"""
function earthsurface(A, projection = DEFPROJ;
    coast = true, kwargs...)
    fig = figure()
    fig.tight_layout()
    ax = subplot(111, projection=projection)
    ax, cb = earthsurface!(ax, A; kwargs...)
    if projection == LONLAT
        ax.set_xticks(-180:60:180, crs=LONLAT)
        ax.set_yticks(-90:30:90, crs=LONLAT)
    end
    if A.name ≠ Symbol("")
        ax.set_title(string(A.name))
    end
    return fig, ax, cb
end

function earthsurface!(ax, A::ClimArray; kwargs...)
    ax.gridlines(alpha = 0.25)
    lon, lat, s = _getwrappeddata(A)
    earthsurface!(ax, lon, lat, s; kwargs...)
end

function earthsurface!(ax, lon, lat, A;
    vmin = minimum(A), vmax = maximum(A), levels = 11,
    coast = true, specific_contour = nothing, kwargs...)

    lon, lat, A = _getwrappeddata(copy(Array(lon)), copy(Array(lat)), A)
    lvls = range(vmin, vmax, length = levels)

    # TODO: Add colorbar as a separate function
    lvls = range(vmin, vmax, length = levels)
    cmap = matplotlib.cm.get_cmap(get(kwargs, :cmap, :viridis), length(lvls)-1)
    # cmap.set_under("k")
    norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
    cax, kw = matplotlib.colorbar.make_axes(ax,location="right",pad=0.02,shrink=0.8)
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    coast && ax.coastlines(; color = "gray")

    ax.contourf(Array(lon), Array(lat), clamp.(A', vmin, vmax), lvls, cmap = cmap,
                transform=LONLAT, vmin=vmin, vmax = vmax)

    if specific_contour isa Array
        ax.contour(Array(lon), Array(lat), clamp.(A', vmin, vmax);
        levels = specific_contour, colors = "r", transform=LONLAT)
    end
    return ax, cb
end

function _getwrappeddata(A::ClimArray)
    lon = copy(Array(dims(A, Lon).val))
    lat = copy(Array(dims(A, Lat).val))
    s = copy(Array(A.data)) # assumes data have Lon x Lat dimensions
    _getwrappeddata(lon, lat, s)
end

function _getwrappeddata(lon, lat, s)
    if lon[end] == 359.5 && lon[1] == 0.5
        # add one more entry to have perfect overlap in plotting
        x = s[1, :]
        s = vcat(Array(s), x')
        push!(lon, 360.5)
    end
    return lon, lat, s
end


_locator = matplotlib.dates.AutoDateLocator(interval_multiples=true)
_locator.intervald
dateformater = matplotlib.dates.ConciseDateFormatter(
    matplotlib.dates.AutoDateLocator(interval_multiples=true)
)

PyPlot.plot(A::ClimArray{T, 1}; kwargs...) where {T} = plot!(gca(), A; kwargs...)
function plot!(ax, A::ClimArray{T, 1}; kwargs...) where {T}
    r = ax.plot(dims(A, 1).val, A.data; kwargs...)
    ax.set_xlabel(_xlabel(A))
    return r
end
function plot!(ax, A::ClimArray{T, 1, <: Tuple{<:Time}}; addlabels = true, kwargs...) where {T}
    t = dims(A, 1)
    x = 1:length(t.val)
    la = ["$(month(a))/$(year(a))" for a in t]
    r = ax.plot(x, A.data; kwargs...)
    ax.set_xlabel(_xlabel(A))
    i = [1, length(la)÷3, 2length(la)÷3, length(la)]
    ax.set_xticks(x[i])
    ax.set_xticklabels(la[i])
    return r
end

function plot!(ax, x, y; kwargs...) where {T}
    r = ax.plot(x, y; kwargs...)
    ax.set_xlabel(_xlabel(x))
end

function heatmap!(ax, A::ClimArray{T, 2}; kwargs...) where {T}
    x = dims(A, 2).val
    y = dims(A, 1).val
    a = A.data#[1:end-1, 1:end-1]
    m = ax.pcolormesh(x, y, a; kwargs...)
    ax.set_xlabel(_xlabel(dims(A)[2]))
    ax.set_ylabel(_xlabel(dims(A)[1]))
    return m
end

_xlabel(x) = ""
_xlabel(x::ClimArray) = _xlabel(dims(x)[1])
_xlabel(::Lon) = "longitude λ [ᵒ]"
_xlabel(::Lat) = "latitude φ [ᵒ]"
_xlabel(::Time) = "time [months]"


LANDMARKS_LON = Dict(
    350 => "West Africa",
    15 => "Europe",
    87 => "Himalaya",
    43 => "Saudi Ar.",
    76 => "India",
    135 => "Australia",
    282 => "Andes",
    300 => "Amazon",
    318 => "Greenland",
    334 => "Atlantic",
    209 => "Alaska",
    165 => "Pacific (W)",
    246 => "Pacific (E)",
    260 => "USA",
)

LANDMARKS_LAT = Dict(
    24 => "North\nAfrica",
    -25 => "Australia\nAndes",
    7 => "ITCZ",
    -5 => "Amazon\nCentral Africa",
    48 => "Europe\nMongolia",
    35 => "Himalaya",
    -2 => "Indian oc.",
    62 => "Canada\nScandinavia",
    -75 => "Antartica",
)

LANDMARKS_LAT_ABS = Dict(abs.(keys(LANDMARKS_LAT)) .=> values(LANDMARKS_LAT))

function landmarks!(ax, lmars)
    i = 0
    for (l, s) in lmars
        ax.axvline(l, ls = "dashed", lw = 1.0, alpha = 0.25, color = "k")
        ax.text(l, 0.98 - (i%3)*0.02, s, transform = ax.get_xaxis_transform(),
        size = 10, va="top", ha = "center")
        i+=1
    end
end

function plot_meanstd!(ax, y; color = "C0", kwargs...)
    m = mean(y); s = std(y)
    ax.axhline(m, lw = 1.5, color = color, ls = "dotted")
    ax.axhspan(m-s, m+s, alpha = 0.1, color = color,
        label = "mean±std: $(rdspl(m)) ± $(rdspl(s))"
    )
    for v in (m-s, m+s)
        ax.axhline(v, lw = 1.0, alpha = 0.5, ls = "solid",
        color = color)
    end
    ax.legend(loc = "lower left")
    return m, s
end


function plottrend!(ax, x, y; label="s", addlabel=true, kwargs...)
    if eltype(x) <: TimeType
        s, o, f, xf = linreg(1:length(x), y)
    else
        s, o, f, xf = linreg(x, y)
    end
    if addlabel
        ax.plot(x, f; label = label*"=$(round(s;sigdigits=3))", kwargs...)
        ax.legend()
    else
        ax.plot(x, f; kwargs...)
    end
    return s, o, f
end

function latitudinal_axis!(ax, llats = -70:20:70)
    ax.set_xticks(sind.(llats))
    ax.set_xticklabels(llats)
    ax.set_xlim(-1,1)
end

#########################################################################
# Color coding
#########################################################################
FIELDCOLORS = Dict(
    "toa_sw_all_mon" => ("C0", "viridis"),
    "toa_sw_clr_t_mon" => ("C3", "YlOrBr_r"),
    "toa_sw_clr_c_mon" => ("C3", "YlOrBr_r"),
    "solar_mon" => ("C5", "plasma"),
    "cldarea_total_daynight_mon" => ("C2", "Blues_r"),
    "sfc_sw_up_all_mon" => ("C4", "RdPu_r"),
    "sfc_sw_down_clr_c_mon" => ("C4", "PuRd_r"),
    "toa_cre_sw_mon" => ("C2", "bone"),
    "toa_sw_all_mon_per" => ("C2", "Blues_r"),
    "seasonal" => ("C2", "Blues_r"),
    "toa_sw_all_mon_res" => ("C4", "plasma"),
    "residual" => ("C4", "plasma"),
    "original" => ("C0", "viridis"),
)

fieldcolors(f) = get(FIELDCOLORS, string(f), ("C0", "binary_r"))


#########################################################################
# Color coding
#########################################################################
function crosscorrelation_plot(x, y, lags; ax = (figure(); gca()))
    cc = crosscor(x, y, lags)
    ax.plot(lags, cc, color = "C1", label = "CrossCor", ls = "dashed")
    ax.set_ylabel("crosscorrelation")
    mi, center, σ = mutualinformation(Array(x), Array(y), lags)
    ax.axhspan(center-3σ, center+3σ, label = "95% null", color = "C1", alpha = 0.25, zorder = -9)
    ax.axhline(center, lw = 1.0, ls = "dashed", color = "C1")
    ax.plot(lags, mi, label = "MI", color = c)
    ax.set_xlabel("delay [months]")
    ax.grid("on")
    ax.legend()
    return ax
end

function hemispheric_component_plot(C)
    Cn, Cs = hemispheric_means(zonalaverage(C))
    hemispheric_component_plot(Cn, Cs)
end

function hemispheric_component_plot(Cn, Cs)
    fig, ax = subplots()
    hemispheric_component_plot!(ax, Cn, Cs)
    return fig, ax
end

function hemispheric_component_plot!(ax, Cn, Cs)
    Cg = 0.5(Cn .+ Cs)
    plot!(ax, Cn, label = "NH, μ=$(rdspl(timeaverage(Cn)))")
    plot!(ax, Cs, label = "SH, μ=$(rdspl(timeaverage(Cs)))", ls = "--")
    plot!(ax, Cg, label = "GL, μ=$(rdspl(timeaverage(Cg)))", ls = ":")
    ax.legend()
    return ax
end
