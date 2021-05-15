Nt′ = 100

@recipe function plot(fc::FusionFC)
	ϕv = LinRange(-fc.ϕj, fc.ϕj, 101)
	ϕc = LinRange(fc.ϕj, 2π - fc.ϕj, 101)

	xv = @. fc.R * cos(π/2 - ϕv)
	zv = @. fc.R * sin(π/2 - ϕv)
	xc = @. fc.R * cos(π/2 - ϕc)
	zc = @. fc.R * sin(π/2 - ϕc)

	color_palette --> palette(:default)[1:2]
	fill := (0, 0, :white)

	label := ["Vesicle" "Cell"]
	aspect_ratio := 1
	
	[xv, xc], [zv, zc]
end

@recipe function plot(kr::FusionKR)
	φv = LinRange(π - kr.φv, π + kr.φv, 101)
	ψc = LinRange(π - kr.ψc, π + kr.ψc, 101)

	junction_offset = kr.Rc * cos(π - kr.ψc) - kr.Rv * cos(π - kr.φv)

	xv = @. kr.Rv * cos(π/2 - φv)
	zv = @. kr.Rv * sin(π/2 - φv) + junction_offset
	xc = @. kr.Rc * cos(π/2 - ψc)
	zc = @. kr.Rc * sin(π/2 - ψc)

	label := ["Vesicle" "Cell"]
	aspect_ratio := 1

	[xv, xc], [zv, zc]
end

@recipe function plot(
	raw::RawOutput;
	n = 0:min(Nt′, length(raw.U))-1)

	pal = [c for c ∈ palette(:blues, Nt′)[Nt′ .- n]] |> reverse
	color_palette := pal
	# seriescolor --> pal
	label --> ""
	xmax = 2raw.s[raw.pj]
	xmaxperc = round(100xmax/raw.s[end], digits = 2)
	xlims --> (0.0, xmax)
	yguide --> "Relative Concentration"
	xguide --> "Arc Length ($xmaxperc% of span)"
	title --> raw.mode * ": Raw Data\n" * "(timesteps $(n[1]) to $(n[end]))"
	raw.s, [raw.U[n]]
end

@recipe function plot(arc::ArcLength)
	t = LinRange(0, arc.tmax, Nt′)
	label --> ""
	title --> arc.mode * ": Arc Length\n" * "($Nt′ equidistant timesteps)"
	color_palette := palette(:blues, length(t))
	xmax = 2arc.sj
	xmaxperc = round(2arc.sj/arc.smax, digits = 2)
	xlims --> (0.0, xmax)
	xguide --> "Arc Length ($xmaxperc% of span)"
	yguide --> "Relative Concentration"
	[s -> arc.u(s, t′) for t′ ∈ t]
end

# TODO: Plots WRT Angle

# @recipe function plot(ang::AngleFC)

# 	ang, π
# end

# @recipe function plot(ang::AngleKR)

# end

# @recipe function plot(::DiffusionAngleSolution, u::Function, ϕmax)
	
# end

@recipe function plot(int::Intensity)
	membranes = [:u, :v, :c]
	label --> ["Total" "Vesicle" "Cell"]
	title --> int.mode * ": Integrated Concentration"
	xguide --> "Time"
	yguide --> "Integrated Concentration"
	xlims --> (0, int.tmax)
	[getproperty(int, m) for m ∈ membranes]
end

density_thickness_ratio = 1/10

# TODO: change to heatmap.
# TODO: sphere or thickness to represent density?

@userplot DensitySliceFC

@recipe function plot(dsfc::DensitySliceFC)
	if length(dsfc.args) ≠ 2
		error("Diffusion plot must have two inputs.")
	end
	
	fc, t = dsfc.args
	if !(fc isa DiffusionFC)
		error("First argument must be a DiffusionFC instance.")
	end
	if !(t isa Float64)
		error("Second argument must be a Float64.")
	end
	
	v(ϕ, t) = fc.ang.v(abs(ϕ), t)
	c(ϕ, t) = fc.ang.c(abs(ϕ), t)
	u(ϕ, t) = abs(ϕ) ≤ fc.ang.ϕj ? v(ϕ, t) : c(ϕ, t)
	
	x(ϕ) = fc.fus.R * (1 + u(ϕ, t)*density_thickness_ratio) * cos(π/2 - ϕ)
	z(ϕ) = fc.fus.R * (1 + u(ϕ, t)*density_thickness_ratio) * sin(π/2 - ϕ)

	# x′(ϕ) = fc.fus.R * cos(π/2 - ϕ)
	# z′(ϕ) = fc.fus.R * sin(π/2 - ϕ)

	density_colour = :purple
	seriescolor := density_colour
	fill := (0, 0.1, density_colour)

	aspect_ratio := 1
	rlims = fc.fus.R * (1 + 2*density_thickness_ratio) .* (-1, 1)
	xlims := 1.5 .* rlims
	ylims := rlims

	axis := nothing
	framestyle := :none

	label := "Density"
	
	x, z, -π, π
end

@userplot DensitySliceKRv

@recipe function plot(dskr::DensitySliceKRv)
	if length(dskr.args) ≠ 2
		error("Diffusion plot must have two inputs.")
	end

	kr, t = dskr.args
	if !(kr isa DiffusionKR)
		error("First argument must be a DiffusionKR instance.")
	end
	if !(t isa Float64)
		error("Second argument must be a Float64.")
	end

	v(φ, t) = kr.ang.v(abs(π - φ), t)

	junction_offset = kr.fus.Rc * cos(π - kr.fus.ψc) - kr.fus.Rv * cos(π - kr.fus.φv)

	xv(φ) = kr.fus.Rv * (1 - v(φ, t)*density_thickness_ratio) * cos(π/2 - φ)
	zv(φ) = kr.fus.Rv * (1 - v(φ, t)*density_thickness_ratio) * sin(π/2 - φ) + junction_offset

	label := "Vesicle"
	aspect_ratio := 1

	xv, zv, π - kr.fus.φv, π + kr.fus.φv
end

@userplot DensitySliceKRc

@recipe function plot(dskr::DensitySliceKRc)
	if length(dskr.args) ≠ 2
		error("Diffusion plot must have two inputs.")
	end

	kr, t = dskr.args
	if !(kr isa DiffusionKR)
		error("First argument must be a DiffusionKR instance.")
	end
	if !(t isa Float64)
		error("Second argument must be a Float64.")
	end

	c(ψ, t) = kr.ang.c(π - abs(π - ψ), t)

	xc(ψ) = kr.fus.Rc * (1 + c(ψ, t)*density_thickness_ratio) * cos(π/2 - ψ)
	zc(ψ) = kr.fus.Rc * (1 + c(ψ, t)*density_thickness_ratio) * sin(π/2 - ψ)

	aspect_ratio := 1
	rlims = kr.fus.Rc * (1 + 2*density_thickness_ratio) .* (-1, 1)
	xlims := 1.5 .* rlims
	ylims := rlims

	label := "Cell"

	xc, zc, π - kr.fus.ψc, π + kr.fus.ψc
end