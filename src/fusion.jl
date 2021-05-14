struct Membrane <: PD
	R::Float64
	D::Float64
end

abstract type FusionMode <: PD end

struct FusionFC <: FusionMode
	ves::Membrane
	cel::Membrane
	R::Float64
	ϕj::Float64

	function FusionFC(ves::Membrane, cel::Membrane)
		Rv = ves.R
		Rc = cel.R

		Rv² = Rv^2
		Rc² = Rc^2
		R′² = Rv² + Rc²
		R′ = √R′²
		# Rs = 2Rv*Rc/R′
		ϕj = acos((Rc² - Rv²)/R′²)

		return new(ves, cel, R′, ϕj)
	end
end

struct FusionKR <: FusionMode
	ves::Membrane
	cel::Membrane
	Rj::Float64
	Rv::Float64
	Rc::Float64
	φv::Float64
	ψc::Float64

	function FusionKR(ves::Membrane, cel::Membrane, Rj::Float64)
		Rv = ves.R
		Rc = cel.R

		Rv² = Rv^2
		Rc² = Rc^2
		Rj² = Rj^2
		Rv′² = 4Rv²^2/(4Rv² - Rj²)
		Rc′² = 4Rc²^2/(4Rc² - Rj²)
		Rv′ = √Rv′²
		Rc′ = √Rc′²
		φv = π - asin(Rj/Rv′)
		ψc = π - asin(Rj/Rc′)

		return new(ves, cel, Rj, Rv′, Rc′, φv, ψc)
	end
end