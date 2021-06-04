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

function transforms(fus::FusionFC)
	# polar angle to arc-length
	ϕ2s(ϕ) = ϕ * fus.R
	s2ϕ(s) = s / fus.R

	# polar angle to vertical
	ϕ2z(ϕ) = fus.R * cos(ϕ)
	z2ϕ(z) = acos(z / fus.R)

	# vertical to shifted vertical
	z2ζ(z) = fus.R - z
	ζ2z(ζ) = fus.R - ζ

	# vertical to horizontal
	z2x(z) = √(fus.R^2 - z^2)
	x2z⁻(x) = √(fus.R^2 - x^2)
	x2z⁺(x) = -x2z⁻(x)

	# two-step
	s2z(s) = s |> s2ϕ |> ϕ2z
	z2s(z) = z |> z2ϕ |> ϕ2s
	ϕ2ζ(ϕ) = ϕ |> ϕ2z |> z2ζ
	ζ2ϕ(ζ) = ζ |> ζ2z |> z2ϕ
	ϕ2x(ϕ) = ϕ |> ϕ2z |> z2x
	x2ϕ⁻(x) = x |> x2z⁻ |> z2ϕ
	x2ϕ⁺(x) = x |> x2z⁺ |> z2ϕ
	ζ2x(ζ) = ζ |> ζ2z |> z2x
	x2ζ⁻(x) = x |> x2z⁻ |> z2ζ
	x2ζ⁺(x) = x |> x2z⁺ |> z2ζ

	# three-step
	s2ζ(s) = s |> s2z |> z2ζ
	ζ2s(ζ) = ζ |> ζ2ϕ |> ϕ2s
	s2x(s) = s |> s2z |> z2x
	x2s⁻(x) = x |> x2ϕ⁻ |> ϕ2s
	x2s⁺(x) = x |> x2ϕ⁺ |> ϕ2s
	
	FC = (
		ϕ2s = ϕ2s,
		s2ϕ = s2ϕ,
		ϕ2z = ϕ2z,
		z2ϕ = z2ϕ,
		z2ζ = z2ζ,
		ζ2z = ζ2z,
		z2x = z2x,
		x2z⁻ = x2z⁻,
		x2z⁺ = x2z⁺,
		s2z = s2z,
		z2s = z2s,
		ϕ2ζ = ϕ2ζ,
		ζ2ϕ = ζ2ϕ,
		ϕ2x = ϕ2x,
		x2ϕ⁻ = x2ϕ⁻,
		x2ϕ⁺ = x2ϕ⁺,
		ζ2x = ζ2x,
		x2ζ⁻ = x2ζ⁻,
		x2ζ⁺ = x2ζ⁺,
		s2ζ = s2ζ,
		ζ2s = ζ2s,
		s2x = s2x,
		x2s⁻ = x2s⁻,
		x2s⁺ = x2s⁺
	)
end

function transformsKRv(fus::FusionKR)
	# polar angle to arc-length
	φ2s(φ) = φ * fus.Rv
	s2φ(s) = s / fus.Rv

	# polar angle to vertical
	φ2z(φ) = fus.Rv * cos(φ)
	z2φ(z) = acos(z / fus.Rv)

	# vertical to shifted vertical
	z2ζ(z) = z - fus.Rv * cos(π - fus.φv)
	ζ2z(ζ) = ζ + fus.Rv * cos(π - fus.φv)

	# vertical to horizontal
	z2x(z) = √(fus.Rv^2 - z^2)
	x2z⁺(x) = √(fus.Rv^2 - x^2)
	x2z⁻(x) = -x2z⁺_v(x)

	# two-step
	s2z(s) = s |> s2φ |> φ2z
	z2s(z) = z |> z2φ |> φ2s
	φ2ζ(φ) = φ |> φ2z |> z2ζ
	ζ2φ(ζ) = ζ |> ζ2z |> z2φ
	φ2x(φ) = φ |> φ2z |> z2x
	x2φ⁻(x) = x |> x2z⁻ |> z2φ
	x2φ⁺(x) = x |> x2z⁺ |> z2φ
	ζ2x(ζ) = ζ |> ζ2z |> z2x
	x2ζ⁻(x) = x |> x2z⁻ |> z2ζ
	x2ζ⁺(x) = x |> x2z⁺ |> z2ζ

	# three-step
	s2ζ(s) = s |> s2z |> z2ζ
	ζ2s(ζ) = ζ |> s2φ |> φ2s
	s2x(s) = s |> s2z |> z2x
	x2s⁻(x) = x |> x2φ⁻ |> φ2s
	x2s⁺(x) = x |> x2φ⁺ |> φ2s

	KRv = (
		φ2s = φ2s,
		s2φ = s2φ,
		φ2z = φ2z,
		z2φ = z2φ,
		z2ζ = z2ζ,
		ζ2z = ζ2z,
		z2x = z2x,
		x2z⁺ = x2z⁺,
		x2z⁻ = x2z⁻,
		s2z = s2z,
		z2s = z2s,
		φ2ζ = φ2ζ,
		ζ2φ = ζ2φ,
		φ2x = φ2x,
		x2φ⁻ = x2φ⁻,
		x2φ⁺ = x2φ⁺,
		ζ2x = ζ2x,
		x2ζ⁻ = x2ζ⁻,
		x2ζ⁺ = x2ζ⁺,
		s2ζ = s2ζ,
		ζ2s = ζ2s,
		s2x = s2x,
		x2s⁻ = x2s⁻,
		x2s⁺ = x2s⁺
	)
end

function transformsKRc(fus::FusionKR)
	# polar angle to arc-length
	ψ2s(ψ) = fus.φv * fus.Rv + fus.Rc * (ψ + fus.ψc - π)
	s2ψ(s) = (s - fus.φv * fus.Rv) / fus.Rc - fus.ψc + π

	# polar angle to vertical
	ψ2z(ψ) = fus.Rc * cos(ψ)
	z2ψ(s) = acos(s / fus.Rc)

	# vertical to shifted vertical
	z2ζ(z) = fus.Rc * cos(π - fus.ψc) - z
	ζ2z(ζ) = fus.Rc * cos(π - fus.ψc) - ζ

	# vertical to horizontal
	z2x(z) = √(fus.Rc^2 - z^2)
	x2z⁻(x) = √(fus.Rc^2 - x^2)
	x2z⁺(x) = -x2z⁻(x)

	# two-step
	s2z(s) = s |> s2ψ |> ψ2z
	z2s(z) = z |> z2ψ |> ψ2s
	ψ2ζ(ψ) = ψ |> ψ2z |> z2ζ
	ζ2ψ(ζ) = ζ |> ζ2z |> z2ψ
	ψ2x(ψ) = ψ |> ψ2z |> z2x
	x2ψ⁻(x) = x |> x2z⁻ |> z2ψ
	x2ψ⁺(x) = x |> x2z⁺ |> z2ψ
	ζ2x(ζ) = ζ |> ζ2z |> z2x
	x2ζ⁻(x) = x |> x2z⁻ |> z2ζ
	x2ζ⁺(x) = x |> x2z⁺ |> z2ζ

	# three-step
	s2ζ(s) = s |> s2z |> z2ζ
	ζ2s(ζ) = ζ |> s2ψ |> ψ2s
	s2x(s) = s |> s2z |> z2x
	x2s⁻(x) = x |> x2ψ⁻ |> ψ2s
	x2s⁺(x) = x |> x2ψ⁺ |> ψ2s

	KRc = (
		ψ2s = ψ2s,
		s2ψ = s2ψ,
		ψ2z = ψ2z,
		z2ψ = z2ψ,
		z2ζ = z2ζ,
		ζ2z = ζ2z,
		z2x = z2x,
		x2z⁻ = x2z⁻,
		x2z⁺ = x2z⁺,
		s2z = s2z,
		z2s = z2s,
		ψ2ζ = ψ2ζ,
		ζ2ψ = ζ2ψ,
		ψ2x = ψ2x,
		x2ψ⁻ = x2ψ⁻,
		x2ψ⁺ = x2ψ⁺,
		ζ2x = ζ2x,
		x2ζ⁻ = x2ζ⁻,
		x2ζ⁺ = x2ζ⁺,
		s2ζ = s2ζ,
		ζ2s = ζ2s,
		s2x = s2x,
		x2s⁻ = x2s⁻,
		x2s⁺ = x2s⁺
	)
end

transforms(fus::FusionKR) = (transformsKRv(fus), transformsKRc(fus))