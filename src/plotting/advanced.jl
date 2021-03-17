include("recipes.jl")
# reduce sizes a bit
PyPlot.rc("font", size = 20) # set default fontsize
PyPlot.rc("axes", labelsize = 22)
PyPlot.rc("legend", fontsize = 18)

using Statistics, StatsBase
# using KernelDensity

function field_comparison_timeplot(f1, f2, l1, l2;
    c = "C0", cmap = "viridis", savename = nothing,
    titlestr = "$l1 vs $l2")

    fig, axs = subplots(2, 3, figsize = (23.39, 16.53)) # A2 page size
    fig.suptitle("temporal analysis | "*titlestr)

    for (i, (v, l)) in enumerate(zip((f1, f2), (l1, l2)))
        isnothing(v) && continue
        co = i == 1 ? c : "C1"
        ls = i == 1 ? "solid" : "dashed"
        gv = globalaverage(v)
        # global timeseries
        ax = axs[1, 1]
        μ, σ = timeagg(mean, gv), timeagg(std, gv)
        plot!(ax, gv, label = l*" μ±σ = $(rdspl(μ,4))±$(rdspl(σ,4))", c = co, ls = ls)
        ax.legend()
        ax.set_ylabel("spatial mean")
        ax.set_xlabel("")
        ax.set_xticklabels([])

        ax = axs[2, 1]
        σv = globalagg(std, v)
        μ, σ = timeagg(mean, σv), timeagg(std, σv)
        plot!(ax, σv, label = l*" μ±σ = $(rdspl(μ,4))±$(rdspl(σ,4))", c = co, ls = ls)
        ax.legend()
        ax.set_ylabel("spatial std")
        i == 1 && ax.set_xlabel("")

        # Yearly variability
        ax = axs[2, 3]
        yearly_variability_plot!(ax, globalaverage(v), i; color = co, ls=ls,
        label = l)
        i == 1 && ax.set_xlabel("")

        # Autocorrelation
        ax = axs[1, 2]
        lags = 0:15
        a = autocor(gv, lags)
        ax.plot(lags, a, label = l, c = co, ls = ls)
        ax.set_ylabel("autocorrelation")
        ax.legend()
        ax.grid("on")
    end

    # cross correlation / mutual information
    ax = axs[2, 2]
    gf1, gf2 = globalaverage.((f1, f2))
    lags = -15:15
    cc = crosscor(gf1, gf2, lags)
    ax.plot(lags, cc, color = "C1", label = "Cor($l1, $l2+lag)", ls = "dashed")
    ax.set_ylabel("crosscorrelation")
    mi, center, σ = mutualinformation(Array(gf1), Array(gf2), lags)
    ax.axhspan(center-3σ, center+3σ, label = "95% null", color = "C1", alpha = 0.25, zorder = -9)
    ax.axhline(center, lw = 1.0, ls = "dashed", color = "C1")
    ax.plot(lags, mi, label = "MI($l1, $l2+lag)", color = c)
    ax.set_xlabel("delay [months]")
    ax.grid("on")
    ax.legend()

    # Time trajectory
    ax = axs[1, 3]
    ax.scatter(gf2, gf1, c = 1:length(t), cmap = cmap, edgecolors = "k")
    ax.plot(gf2, gf1, color = c, alpha = 0.75, lw = 0.5, zorder = -9)
    ax.set_ylabel(l1)
    ax.set_xlabel(l2)
    # correlation textboxes:
    rMI = relativeMI(gf2, gf1)
    ρ = corspearman(gf2, gf1)
    ax.text(
        0.95, 0.95, "time trajectory\nρ=$(rdspl(ρ))\nrMI=$(rdspl(rMI))", va = "top",
        bbox = bbox2, transform = ax.transAxes, ha = "right"
    )

    fig.tight_layout()
    fig.subplots_adjust(hspace = 0.13, top = 0.93)

    cax, cb, = add_colorbar!(axs[1,3], 1:length(t);
    name = "months", location = "top", c = cmap)
    !isnothing(savename) && wsave(savename*"_time.png", fig)
    return fig
end

function yearly_variability_plot!(ax, v, i = 1, t = Array(dims(v, Time));
    doscatter=true, shift=true, ploterror = true, kwargs...)
    ut, fluct = yearly_variability(t, v)
    means = mean.(fluct)
    stds = std.(fluct)
    meanσ = mean(stds)
    max = (1:12)
    shift && (max .*= (i-1)*0.1)
    M = maximum(length(l) for l in fluct)
    if ploterror == true
        ax.errorbar(max, means, yerr = stds; marker = "s", ms = 5,
        label = get(kwargs, :label, "")*" \$\\langle\\sigma_\\mathrm{month}\\rangle\$ = $(rdspl(meanσ))",
        color = get(kwargs, :color, "C$(i-1)"),
        kwargs...)
    elseif ploterror == :fill
        ax.plot(max, means;
        label = get(kwargs, :label, "")*" \$\\langle\\sigma_\\mathrm{month}\\rangle\$ = $(rdspl(meanσ))",
        color = get(kwargs, :color, "C$(i-1)"),
        kwargs...)
        ax.fill_between(max, means .- stds, means .+ stds;
        color = get(kwargs, :color, "C$(i-1)"), alpha = 0.5,
        kwargs...)
    else
        ax.plot(max, means;
        label = get(kwargs, :label, "")*" \$\\langle\\sigma_\\mathrm{month}\\rangle\$ = $(rdspl(meanσ))",
        color = get(kwargs, :color, "C$(i-1)"),
        kwargs...)
    end
    if doscatter
        for (j, m) in enumerate(max)
            f = fluct[j]
            ax.plot(fill(m, length(f)), f, marker = "o", ls = "none",
            ms = 6, color = get(kwargs, :color, "C$(i-1)"),
            alpha = 5/M, mew = 0, #fillstyle = "none",
            )
        end
    end
    ax.legend()
    ax.set_ylabel("yearly variability")
    ax.set_xticks(2:3:12)
    ax.set_xticklabels([monthname(τ)[1:3] for τ in t[2:3:12]])
end

function field_comparison_spaceplot(f1, f2, l1, l2;
    c = "C0", cmap = "viridis", savename = nothing, reverselat = false,
    titlestr = "$l1 vs $l2")

    fig, axs = subplots(2, 3, figsize = (23.39, 16.53)) # A2 page size
    fig.suptitle("spatial analysis | "*titlestr)

    for (i, (v, l)) in enumerate(zip((f1, f2), (l1, l2)))
        co = i == 1 ? c : "C1"
        ls = i == 1 ? "solid" : "dashed"
        tv = timeaverage(v)
        if hasdim(v, Lat)
            lv = lataverage(tv)
            # global timeseries
            ax = axs[1, 1]
            μ, σ = mean(lv), std(lv)
            plot!(ax, lv, label = l*", μ±σ = $(rdspl(μ, 4))±$(rdspl(σ))", c = co, ls = ls)
            ax.legend()
            ax.set_ylabel("time & lat mean")
            i == 1 && ax.set_xlabel("")

            # Pdfs
            ax = axs[2, 1]
            val = lv .- mean(lv)
            kde = KernelDensity.kde(val; npoints = 2048*4)
            N = 50
            fs = range(minimum(val), maximum(val); length = 100)
            ax.hist(val.data, N, density = true, color = co, alpha = 0.5)
            ax.plot(fs, KernelDensity.pdf(kde, fs), color = co, label = l, ls = ls)
            ax.legend()
            ax.set_yticks([])
            ax.set_ylabel("pdf of time & lat mean")
        end

        # latitudes
        ax = axs[1, 2]
        zv = zonalaverage(tv)
        ax.plot(sind.(Array(dims(zv, Lat))), zv, ls = ls, color=co, label = l)
        ax.set_xlabel("")
        ax.set_ylabel("time & lon mean")
        ax.legend()

        ax = axs[i, 3]
        yv = timeaverage(v)
        if spacestructure(yv) == Grid()
            lats = sind.(Array(dims(yv, Lat)))
            i == 2 && reverselat && (lats .*= -1)
            p = ax.pcolormesh(Array(dims(yv, Lon)), lats, Array(yv)', cmap = cmap)
            ax.set_aspect("auto")
            fig.colorbar(p, ax = ax)
        elseif spacestructure(yv) == EqArea()
            lons = [x[1] for x in dims(yv, Coord)]
            lats = [x[2] for x in dims(yv, Coord)]
            i == 2 && reverselat && (lats .*= -1)
            ax.scatter(lons, lats, c = Array(yv), cmap = cmap)
            cax, cb = add_colorbar!(ax, range(minimum(yv), maximum(yv); length = 10))
        end
        μ = globalagg(mean, yv)
        σ = globalagg(std, yv)
        ax.text(-0.1, 1.0, "time mean "*l*", μ±σ = $(rdspl(μ))±$(rdspl(σ))",
        transform = ax.transAxes, va = "bottom")
        ax.set_ylabel("\$\\sin(\\phi)\$")
        i == 2 && ax.set_xlabel("longitude λ [ᵒ]")
    end

    # latitude difference
    f1lat = zonalaverage(timeaverage(f1)) |> Array
    f2lat = zonalaverage(timeaverage(f2)) |> Array
    lat = dims(zonalaverage(f1), Lat) |> Array
    ax = axs[2, 2]
    ax.plot(sind.(lat), f1lat .- f2lat, color = c, label = "\$\\psi\$ = NH - SH")
    ax.plot(sind.(lat), (f1lat .- f2lat) .* cosd.(lat), color = "C1",
            alpha = 0.75, label = "\$\\psi \\cos\\phi\$")
    ax.set_xlabel("\$\\sin(\\phi)\$")
    latavg = lataverage(DimensionalArray(f1lat .- f2lat, (Lat(lat),)))
    hemavg = sum(globalaverage.(timeaverage.((f1, f2))))/2
    ax.axhline(latavg, ls = "dashed", color = "C1",
    label = "\$\\langle\\psi\\rangle_\\phi\$ = $(rdspl(latavg)) ($(rdspl(100*abs(latavg)/hemavg))%)")
    ax.legend()
    ax.set_ylabel("time & lon mean")

    fig.tight_layout()
    fig.subplots_adjust(hspace = 0.22, top = 0.91, bottom = 0.08)
    !isnothing(savename) && wsave(savename*"_space.png", fig)
    return fig
end

using FFTW
function surrogate_plot(args...; normalize = false)
    fig = figure()
    ax1 = subplot(2, 1, 1)
    ax1.set_ylabel("\$x\$")
    ax2 = subplot(2, 3, 4)
    ax2.set_ylabel("AC(\$x\$)")
    ax2.set_xlabel("lag")
    ax3 = subplot(2, 3, 5)
    ax3.set_ylabel("\$|\\mathcal{F}(x) |^2\$")
    ax3.set_xlabel("period (in points)")
    ax4 = subplot(2, 3, 6)
    ax4.set_ylabel("pdf(\$x\$)")
    ax4.set_xlabel("\$x\$")
    for (i, x) in enumerate(args)
        normalize && (x = (x .- mean(x))./std(x))
        c = "C$(i-1)"
        ax1.plot(x, color = c, alpha = 0.75, ls = LINES[i])
        ax1.axhline(mean(x), color = c, alpha = 0.75, ls = LINES[i])

        a = autocor(x)
        ax2.plot(0:length(a)-1, a, c=c, alpha = 0.75, ls = LINES[i])

        P = abs.(rfft(x .- mean(x)))
        f = rfftfreq(length(x))
        ps = 1 ./ Float64.(f)
        r = 10:(length(ps))
        ax3.plot(ps[r], P[r] ./ maximum(P[r]), marker = "o",
        color = c, zorder = -99, ms = 4, ls = LINES[i], alpha = 0.75)
        ax3.set_yticks([])

        # kde = KernelDensity.kde(x; npoints = 2048*2)
        fs = range(minimum(x), maximum(x); length = 100)
        ax4.hist(x, 50, density = true, color = c, alpha = 0.5)
        # ax4.plot(fs, KernelDensity.pdf(kde, fs), color = c, alpha = 0.75, ls = LINES[i])
        ax4.set_yticks([])
    end
    fig.tight_layout()
    return fig, ax1
end

using Statistics, StatsBase
function diff_corr_plot(resnhA, resshA, t = 1:length(resnhA))
    lags, lagsa = -15:15, 0:15
    fig = figure()
    ax1 = subplot(2,1,1)
    ax1.plot(t, resnhA .- resshA, label = "difference")
    ax1.axhline(0, label="0", color="C1")
    ax1.set_ylabel("difference"); ax1.legend()

    ax2 = subplot(2,2,3)
    ax2.plot(lagsa, autocor(resnhA, lagsa), label = "AC, NH")
    ax2.plot(lagsa, autocor(resshA, lagsa), label = "AC, SH", ls = "dashed")
    ax2.set_ylabel("correlation"); ax2.legend()
    ax3 = subplot(2,2,4)
    ax3.plot(lags, crosscor(resnhA, resshA, lags), label = "CC, SH")
    mi, center, σ = mutualinformation(Array(resnhA), Array(resshA), lags)
    ax3.axhspan(center-2σ, center+2σ, label = "95% null", color = "C1", alpha = 0.25, zorder = -9)
    ax3.axhline(center, lw = 1.0, ls = "dashed", color = "C1")
    ax3.plot(lags, mi, label = "MI")
    ax3.grid("on")
    ax3.legend()
end

function corr_plot(resnhA, resshA, t = 1:length(resnhA); lx = "x", ly = "y")
    lags, lagsa = -15:15, 0:15
    fig = figure()
    ax1 = subplot(2,1,1)
    x = Array((resnhA .- mean(resnhA)) ./ std(resnhA))
    y = Array((resshA .- mean(resshA)) ./ std(resshA))
    ax1.plot(t, x, label = lx)
    ax1.plot(t, y, label = ly, ls = "--")
    ax1.set_ylabel("standarized signal"); ax1.legend()

    ax2 = subplot(2,2,3)
    ax2.plot(lagsa, autocor(x, lagsa), label = "AC, $(lx)")
    ax2.plot(lagsa, autocor(y, lagsa), label = "AC, $(ly)", ls = "dashed")
    ax2.set_ylabel("correlation"); ax2.legend()
    ax3 = subplot(2,2,4)
    ax3.plot(lags, crosscor(x, y, lags), label = "CC")
    mi, center, σ = mutualinformation(Array(x), Array(y), lags)
    ax3.axhspan(center-2σ, center+2σ, label = "95% null", color = "C1", alpha = 0.25, zorder = -9)
    ax3.axhline(center, lw = 1.0, ls = "dashed", color = "C1")
    ax3.plot(lags, mi, label = "MI")
    ax3.grid("on")
    ax3.legend()
end
