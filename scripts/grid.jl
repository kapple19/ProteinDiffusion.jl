pj = 6
P = 11
sj = π/3
sP = π/2

function spatial_grid(p)
	p == pj && return sj
	p < pj && return sj - sj / pj^2 * (p - pj)^2
	p > pj && return sj + (sP - sj) / (P - pj)^2 * (p - pj)^2
	return NaN
end

s = [spatial_grid(p) for p ∈ 0:P]

##
using Plots

plot(0:P, s)
scatter!([0, pj, P], [0, sj, sP])

##
scatter(0:P, s)
