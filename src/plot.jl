Nt = 101

## Full Fusion
@recipe function f(raw::FullDiffusionRaw)
	title := "Full Fusion: Raw Data"
	yguide := "Relative Density"
	xguide := "Polar Angle [rad]"
	legend := false
	
	N = length(raw.t) - 1
	T = raw.t[end]
	ns = [
		findnearest(raw.t, t)[1]
		for t ∈ LinRange(0, T, Nt)
	] |> unique |> sort
	
	color_palette := cgrad(:blues, length(ns))

	raw.ϕ, [raw.U[n] for n ∈ ns]
end

function select_concentration_full_fusion(membrane, ang)
	sys = (:system, :sys, :s)
	ves = (:vesicle, :ves, :v)
	cel = (:cell, :cel, :c)

	if membrane ∉ (sys..., ves..., cel...)
		error("Unrecognised name for membrane.")
	end
	
	function concentration(s, t)
		membrane ∈ sys && return ang.u(s, t)
		membrane ∈ ves && return ang.v(s, t)
		membrane ∈ cel && return ang.c(s, t)
	end
	
	return [
		s -> concentration(s, t)
		for t = LinRange(0.0, ang.tmax, Nt)
	]
end

@recipe function f(
	ang::FullDiffusionAngle,
	membrane = :system;
	xlim = (0, 2ang.ϕj),
	ylim = (0, 1))

	title := "Full Fusion: WRT Angle"
	yguide := "Relative Concentration"
	xguide := "Polar Angle [rad]  (first $(200ang.ϕj/π) % of domain)"
	legend := false

	xlims := xlim
	ylims := ylim

	color_palette := cgrad(:blues, Nt)

	select_concentration_full_fusion(membrane, ang)
end

## KNR Fusion
@recipe function f(raw::KNRDiffusionRaw)
	title := "KNR Fusion: Raw Data"
	yguide := "Relative Density"
	xguide := "Arc Length [rad]"
	legend := false
	
	N = length(raw.t) - 1
	T = raw.t[end]
	ns = [
		findnearest(raw.t, t)[1]
		for t ∈ LinRange(0, T, Nt)
	] |> unique |> sort
	
	color_palette := cgrad(:blues, length(ns))

	raw.s, [raw.U[n] for n ∈ ns]
end

function select_concentration_knr_fusion(membrane, arc)
	sys = (:system, :sys, :s)
	ves = (:vesicle, :ves, :v)
	cel = (:cell, :cel, :c)

	if membrane ∉ (sys..., ves..., cel...)
		error("Unrecognised name for membrane.")
	end

	function concentration(s, t)
		membrane ∈ sys && return arc.u(s, t)
		membrane ∈ ves && return arc.v(s, t)
		membrane ∈ cel && return arc.c(s, t)
	end
	
	return [
		s -> concentration(s, t)
		for t = LinRange(0.0, arc.tmax, Nt)
	]
end

@recipe function f(
	arc::KNRDiffusionArc,
	membrane = :system;
	xlim = (0, 2arc.sj),
	ylim = (0, 1))

	title := "KNR Fusion: WRT Arc Length"
	yguide := "Relative Concentration"
	xguide := "Arc Length (first $(200arc.sj/arc.S) % of domain)"
	legend := false

	xlims := xlim
	ylims := ylim

	color_palette := palette(:blues, Nt)

	select_concentration_knr_fusion(membrane, arc)
end

function intensity_title(int::PDIntensity)
	int isa FullDiffusionIntensity && return "Full Fusion: Intensity"
	int isa KNRDiffusionIntensity && return "KNR Fusion: Intensity"
	error("Unrecognised struct.")
end

## Both
@recipe function f(int::PDIntensity)
	title := intensity_title(int)
	yguide := "Integrated Concentration"
	xguide := "Time"
	legend := :outerbottom
	label := ["Total" "Vesicle" "Cell"]

	xlims := (0, int.tmax)

	[int.I, int.v, int.c]
end
