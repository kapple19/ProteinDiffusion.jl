struct Raw <: PD
	fusion::String
	s::OVector64
	t::OVector64
	U::NOVector64
	pj::Int64
end

struct Diffusion <: PD
	fusion::String
	raw::Raw
end