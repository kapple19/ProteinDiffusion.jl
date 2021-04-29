##
using Plots
using RecipesBase

##
struct Vesicle
	Rv::Float64
	Rc::Float64
end

##
@recipe function plot(v::Vesicle)
	xlabel("x")
	ylabel("z")
end

##
plot

##