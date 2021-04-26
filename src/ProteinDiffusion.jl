module ProteinDiffusion

import OffsetArrays: OffsetArray, OffsetVector, Origin
import LinearAlgebra: SymTridiagonal
import QuadGK: quadgk
import Interpolations: Gridded, Linear, interpolate
import IntervalArithmetic: (..)
import Plots: palette
import RecipesBase: RecipesBase, @recipe

export full_fusion
export knr_fusion
export fusion

include("preamble.jl")
include("auxiliary.jl")
include("fem.jl")
include("fusion.jl")
include("diffusion.jl")
include("plot.jl")

end # module