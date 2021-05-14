##
using RecipesBase
using Plots

##
struct Vec
	x::Vector{Float64}
	y::Vector{Float64}

	function Vec()
		x = LinRange(0, 10, 13)
		y = randn(13)
		new(x, y)
	end
end

struct Fcn
	ω::Float64
	xlims::NTuple{2, Float64}
	f::Function

	function Fcn()
		ω = 1 + randn()/2
		f(x) = sin(ω*x)
		xlims = (0.0, 10.0)
		new(ω, xlims, f)
	end
end

(fcn::Fcn)(args...) = fcn.f(args...)

@recipe function plot_recipe(vec::Vec)
	linecolor --> :blue
	xlims --> extrema(vec.x)
	vec.x, vec.y
end

@recipe function plot_recipe(::Type{Fcn}, fcn::Fcn)
	linecolor --> :red
	xlims --> fcn.xlims
	fcn.f
end

##
a = Vec()
f = Fcn()

@show a.x
@show a.y
@show f(0)

##
plot(f)
scatter!(a)
