#=
Style file for plotting
=#
using PyPlot, PyCall, DrWatson
using Printf
DrWatson._wsave(s, fig::Figure) = fig.savefig(s, dpi = 600, transparent = false)

"""
    rdspl(x, n = 4)
Round `x` with `n` sigdigits for display purposes.
"""
rdspl(x::Real, n = 4) = round(x, sigdigits=n)
rdspl(x::AbstractVector, n = 3) = Tuple((round.(Float64.(x); sigdigits=n)))

# Font related
PyPlot.rc("font", size = 26) # set default fontsize
PyPlot.rc("axes", labelsize = 32)
PyPlot.rc("legend", fontsize = 22)
# PyPlot.rc("font", family = "Times New Roman") # Serif main font
PyPlot.rc("font", family = "DejaVu Sans") # sans main font
# PyPlot.rc("mathtext", rm = "sanserif", fontset="dejavusans") # sans math font
PyPlot.rc("mathtext", rm = "serif", fontset="cm") # serif math font

for z in ("x", "y")
    PyPlot.rc("$(z)tick.major", size = 7, width = 1.2)
    PyPlot.rc("$(z)tick.minor", size = 3, visible = false)
end

figx = 16 # width correspoding to full text width
figy = 8 # height corresponding to 2 lines of plots
# default figure size full page width, 1 plot height
PyPlot.rc("figure", figsize = (figx, figy))

if false # test font in default figure size
    figure(); plot(1:10, label = "\$a=\\int_0^\\infty xdx\$")
    xlabel("x label"); ylabel("\$H_2\$"); legend(); tight_layout()
end

# Plot related
PyPlot.rc("errorbar", capsize = 6)
PyPlot.rc("lines", lw = 2.0, markersize = 7.5)
PyPlot.rc("axes", grid = true)
PyPlot.rc("grid", color = "0.75", alpha = 0.5)

struct CyclicContainer <: AbstractVector{String}
    c::Vector{String}
    n::Int
end
CyclicContainer(c) = CyclicContainer(c, 0)

Base.length(c::CyclicContainer) = typemax(Int)
Base.size(c::CyclicContainer) = size(c.c)
Base.iterate(c::CyclicContainer, state=1) = Base.iterate(c.c, state)
Base.getindex(c::CyclicContainer, i) = c.c[(i-1)%length(c.c) + 1]
function Base.getindex(c::CyclicContainer)
    c.n += 1
    c[c.n]
end
Base.iterate(c::CyclicContainer, i = 1) = iterate(c.c, i)

COLORS = CyclicContainer(["#008080", "#101636", "#1e8fff",
"#d07b17", "#6f0d4d", "#add54d"])
MARKERS = CyclicContainer(["o", "s", "^", "p", "P", "D", "X"])
LINES = CyclicContainer(["-", "--", ":", "-."])
# Also set default color cycle
PyPlot.rc("axes", prop_cycle = matplotlib.cycler(color=COLORS.c))

# Default colormap
CMAP = matplotlib.cm.get_cmap("magma")

# cool sequential maps are:
# :Tokyo, :Devon, :Bamako, :Batlow, :Bilbao, :Davos, :Turku
# sequential(s) =
# palettable = pyimport("palettable")
# TOKYO = palettable.scientific.sequential.Tokyo_20.mpl_colormap
# DEVON = palettable.scientific.sequential.Devon_20.mpl_colormap
# DEVON = palettable.scientific.sequential.Bamako_20.mpl_colormap
#
# # Cool diverging maps are:
# VIK = palettable.scientific.diverging.Vik_20.mpl_colormap


export COLORS, MARKERS, LINES, CMAP

if false
    figure(figsize = (20, 15)) # show colors
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(212)
    lw = 60
    for (i, c) in enumerate(COLORS)
        chsv = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(c))
        ax1.plot([0, 1], [0, 0] .+ i, color = c, lw = lw)
        ax1.set_title("color")
        ax2.plot([0, 1], [0, 0] .+ i, color = string(chsv[3]), lw = lw)
        ax2.set_title("brightness")
        ax3.plot([0, 1], [0, 0] .+ i, color = string(chsv[2]), lw = lw)
        ax3.set_title("saturation")
        local x = 0:0.05:10π
        ax4.plot(x, cos.(x .+ i/2) .+ rand(length(x))/2, lw = 2)
    end
end

bbox = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 1.0)
bbox2 = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 0.5)
function add_heatmap_legend!(ax, s; loc = (0.5, 0.9))
    ax.text(loc..., s, ha="center",
    transform = ax.transAxes, bbox = bbox2, zorder = 99)
end

function add_identifiers!(fig = gcf())
    bbox = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 1.0)
    for (i, ax) in enumerate(fig.get_axes())
        l = collect('a':'z')[i]
        loc = (0.985, 0.975)
        try
            ax.text(loc..., "$(l)",
            transform = ax.transAxes, bbox = bbox, zorder = 99)
        catch err
            ax.text(loc..., 1, "$(l)",
            transform = ax.transAxes, bbox = bbox, zorder = 99)
        end
    end
end

function coolfill!(ax, x, y, dy, c, label = "")
    α = 0.25
    ax.plot(x, y, label = label, color = c, lw = 2.0)
    ax.fill_between(x, y .- dy, y .+ dy, alpha = α, color = c)
    lw = 0.5
    α2 = 0.5
    ax.plot(x, y .+ dy,  color = c, lw = lw, alpha = α2)
    ax.plot(x, y .- dy,  color = c, lw = lw, alpha = α2)
end

function nice_arrow!(ax, xc, yc, xspan, yspan;
    style = "<->", tex = "", xo = 0.2, yo = -0.2)
    ax.annotate("",  xy=(xc-xspan/2, yc - yspan/2), xycoords="data",
                xytext=(xc+xspan/2, yc + yspan/2), textcoords="data",
                arrowprops = (Dict(:arrowstyle=>style,
                                :connectionstyle=>"arc3",
                                :lw=>1.5, :facecolor => "black")), zorder = 99)
    if tex != ""
        ax.text(xc + xo, yc + yo, tex, size = 24)
    end
end

function coolhist!(ax, data, bins, color, label = "", alpha = 0.25)
    h, b, = ax.hist(data, bins, density = true, color = color,
    alpha = alpha)

    b = 0.5(b[1:end-1] .+ b[2:end])
    ax.plot(b, h, color = color, lw = 1.0, label = label)
end

function add_grid!(ax, nx::Int, ny = nx; kwargs...)
    @assert n ≥ 3
    x = ax.get_xlim()
    y = ax.get_ylim()
    dx = (x[2]-x[1])/nx; dy = (y[2]-y[1])/ny
    for i in 1:(n-1)
        ax.axhline(y[1]+n*dy, color = "gray", alpha = 0.5, kwargs...)
        ax.axvline(x[1]+n*dx, color = "gray", alpha = 0.5, kwargs...)
    end
end

"""
    axis_zoomin!(zoomin, origin, zbox, rbox, co)
Create a zoomin box connecting two axes, the `zoomin` and `origin`.
The `zoomin` box is in the origin axis, while the `zoomin` axis is the
box. `rbox` is the enclosing box of the `zooming` axes, while `zbox`
is the small "zoom-in" box of the `origin` axis. They must be in the form
((x1, y1), (x2, y2)). `co` is color.
"""
function axis_zoomin!(zoomin, origin, zbox, rbox, co)
    # plot box in zoomin axis
    line, = zoomin.plot(
    [rbox[1][1], rbox[2][1], rbox[2][1], rbox[1][1], rbox[1][1]],
    [rbox[1][2], rbox[1][2], rbox[2][2], rbox[2][2], rbox[1][2]],
    color=co,  lw = 2.0)
    line.set_clip_on(false)

    line, = origin.plot(
    [zbox[1][1], zbox[2][1], zbox[2][1], zbox[1][1], zbox[1][1]],
    [zbox[1][2], zbox[1][2], zbox[2][2], zbox[2][2], zbox[1][2]],
    color=co,  lw = 2.0)
    line.set_clip_on(false)

    for e in 1:2
        con = matplotlib.patches.ConnectionPatch(
        xyA = (rbox[1][1], rbox[e][2]), xyB=(zbox[2][1], zbox[e][2]),
        coordsA="data", coordsB="data",
        axesA = zoomin, axesB=origin, color=co, lw = 2.0)

        zoomin.add_artist(con)
    end
end

function add_colorbar!(ax, lvls; c = :viridis, location = "right", name = nothing)
    cmap = matplotlib.cm.get_cmap(c, length(unique(lvls))-1)
    norm = matplotlib.colors.Normalize(vmin=minimum(lvls), vmax=maximum(lvls))
    or = location == "right" ? "vertical" : "horizontal"
    cax, kw = matplotlib.colorbar.make_axes(ax,
        pad=0.01, shrink=0.8, location = location
    )
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation = or)
    location == "top" && cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_tick_params(size = 3, labelsize = 10)
    if !isnothing(name)
        if location == "top"
            cax.text(-0.01, 0.5, name, transform = cax.transAxes,
            va = "center", ha = "right")
        elseif location == "right"
            cb.set_label(name)
        end
    end
    return cax, cb, cmap
end


colorsys = pyimport("colorsys")
"""
    lighten(color, amount = 0.5)
Lighten the color by the given amount, which will darken if less than 1.0.
(i.e. amount multiplies current luminosity)
"""
function lighten(color, amount=0.5)
    try
        c = matplotlib.colors.cnames[color]
    catch err
        c = color
    end
    c = colorsys.rgb_to_hls(matplotlib.colors.to_rgb(c)...)
    return colorsys.hls_to_rgb(c[1], max(0, min(1, amount * c[2])), c[3])
end
