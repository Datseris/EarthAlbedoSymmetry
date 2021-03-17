module AlbedoProperties

using Reexport
@reexport using Dates, Statistics, DrWatson, ClimateBase

using NCDatasets: NCDataset, dimnames, NCDatasets
export NCDataset, dimnames

Time = ClimateBase.Time

using DimensionalData
export DimensionalData
AbDimArray = DimensionalData.AbstractDimArray
export AbDimArray
include("core.jl")
include("statistics.jl")
include("residuals.jl")
include("albedo.jl")

end
