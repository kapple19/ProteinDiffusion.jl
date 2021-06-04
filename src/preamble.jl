Vector64 = Vector{Float64}

NVector64 = Vector{Vector64}

OVector64 = OffsetVector{Float64}

NOVector64 = OffsetVector{
	OffsetVector{Float64, Vector64},
	Vector{OffsetVector{Float64, Vector64}}
}

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