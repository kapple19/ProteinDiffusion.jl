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
import IntervalArithmetic: (..)
import Plots: palette
import RecipesBase:
	RecipesBase,
	@recipe

export Membrane
export FullFusion
export KNRFusion

include("preamble.jl")
include("auxiliary.jl")
include("fem.jl")
include("diffusion.jl")
include("fusion.jl")
include("plot.jl")

end # module