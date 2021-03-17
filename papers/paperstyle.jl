using DrWatson
@quickactivate :AlbedoProperties

include(srcdir("plotting", "recipes.jl"))

# width correspoding to half-column width
figx = 12
figy = 8
PyPlot.rc("figure", figsize = (figx, figy))
PyPlot.rc("legend", fontsize = 24, framealpha = 0.8, handlelength = 1)
PyPlot.rc("grid", color = "0.75", alpha = 0.5, lw = 1.5)
