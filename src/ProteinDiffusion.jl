module ProteinDiffusion
using OffsetArrays: OffsetArray, OffsetVector, Origin
using LinearAlgebra: SymTridiagonal
using QuadGK: quadgk
using Interpolations: Gridded, Linear, interpolate
using Plots: cgrad
using RecipesBase

export full_fusion
export knr_fusion
export fusion

abstract type PD <: Any end

Base.show(io::IO, pd::PD) = println(io, pd |> typeof |> string)

include("auxiliary.jl")
include("fem.jl")
include("diffusion.jl")
include("fusion.jl")
include("plot.jl")

end