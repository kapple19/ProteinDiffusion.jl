"""
`PD`

Generic supertype for all types involved in this protein diffusion model.
"""
abstract type PD end

function Base.show(io::IO, pd::PD)
	print(io, pd |> typeof |> string)
	if hasproperty(pd, :fusion)
		println(": ", pd.fusion)
	else
		println()
	end
	println(io, pd |> propertynames |> string)
end

"""
`NOVector64`

Nested offset vector of double precision floating point values, i.e. an offset vector of offset vectors each of float64s.

For instance `V`, elements are called as `V[n][p]`.
"""
NOVector64 = OffsetVector{OffsetVector{Float64, Vector{Float64}}, Vector{OffsetVector{Float64, Vector{Float64}}}}

"""
`OVector64`

Offset vector of float64s.
"""
OVector64 = OffsetVector{Float64}
