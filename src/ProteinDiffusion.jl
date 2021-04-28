module ProteinDiffusion

import OffsetArrays: OffsetArray, OffsetVector, Origin
import LinearAlgebra: SymTridiagonal
import QuadGK: quadgk
import NumericalIntegration: integrate
import Interpolations: interpolate, Gridded, LinearInterpolation, Linear
import IntervalArithmetic: (..)
import Plots: palette
import RecipesBase: RecipesBase, @recipe

export full_fusion
export knr_fusion
export fusion
export Comparison

include("preamble.jl")
include("auxiliary.jl")
include("fem.jl")
include("fusion.jl")
include("diffusion.jl")
include("plot.jl")

end # module