module ProteinDiffusion

import OffsetArrays:
	OffsetArray,
	OffsetVector,
	Origin
import LinearAlgebra: SymTridiagonal
import QuadGK: quadgk
import NumericalIntegration: integrate
import Interpolations:
	interpolate,
	Gridded,
	LinearInterpolation,
	Linear
# import Intervals:
# 	Interval,
# 	Open

import Plots:
	palette,
	RGB
import RecipesBase:
	RecipesBase,
	@recipe,
	@userplot

export Membrane
export FusionFC
export FusionKR
export DiffusionFC
export DiffusionKR

include("preamble.jl")
include("auxiliary.jl")
include("fusion.jl")
include("fem.jl")
include("diffusion.jl")
include("plot.jl")

end # module