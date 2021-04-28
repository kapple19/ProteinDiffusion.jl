### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ fc567100-a675-11eb-1f44-f1a8c9822102
begin
	using ProteinDiffusion
	using PlutoUI
	using Plots
end

# ╔═╡ 1c479d20-2611-4a9e-975b-aaf7dbf867f0
md"Vesicle Radius (log scale): $(@bind pRv Slider(0:0.1:3, default = 1))"

# ╔═╡ 3d9cd9a3-a595-4b40-af33-76cedc7bfde1
Rv = 10.0^pRv

# ╔═╡ 90454198-763c-48cb-9bc7-9c2f57b0d9b4
md"Cell Radius (log scale): $(@bind pRc Slider(1:0.1:3, default = 1))"

# ╔═╡ bf481f6d-44e9-4a9b-8cbd-0950e0e5fad4
Rc = 10.0^pRc

# ╔═╡ d295c4f2-cd80-4d86-9c8c-17f6ec9dbc9d
md"Junction Radius (log scale): $(@bind pRj Slider(0:0.1:3, default = 1))"

# ╔═╡ f5236744-0e1e-4b99-92c8-b2a7c0279186
Rj = 10.0^pRj

# ╔═╡ 12ccf1c2-de78-4528-bf7e-d078360a2dae
md"Vesicle Diffusivity: $(@bind Dv NumberField(0.001:10))"

# ╔═╡ 52fdd5b0-bb92-464d-a79b-5f956441528a
md"Cell Diffusivity: $(@bind Dc NumberField(0.001:10))"

# ╔═╡ 32a76355-0e85-48e8-8a98-dc209398e34c
f = full_fusion(Rv, Rc, Dv, Dc)

# ╔═╡ 8d61f6cc-784e-4cc8-9e21-070da8c84a1b
k = knr_fusion(Rv, Rc, Rj, Dv, Dc)

# ╔═╡ 60719c7d-8e19-4cdf-aeda-3ff83f446922
plot(
	plot(f.arc),
	plot(k.arc),
	layout = (1, 2)
)

# ╔═╡ 15460682-a64a-41d6-934c-c66f2fb40f75
plot(f.int)

# ╔═╡ de169125-f359-408c-83af-5aebb5fcf2cc
plot(k.int)

# ╔═╡ Cell order:
# ╠═fc567100-a675-11eb-1f44-f1a8c9822102
# ╠═1c479d20-2611-4a9e-975b-aaf7dbf867f0
# ╟─3d9cd9a3-a595-4b40-af33-76cedc7bfde1
# ╠═90454198-763c-48cb-9bc7-9c2f57b0d9b4
# ╟─bf481f6d-44e9-4a9b-8cbd-0950e0e5fad4
# ╠═d295c4f2-cd80-4d86-9c8c-17f6ec9dbc9d
# ╟─f5236744-0e1e-4b99-92c8-b2a7c0279186
# ╟─12ccf1c2-de78-4528-bf7e-d078360a2dae
# ╟─52fdd5b0-bb92-464d-a79b-5f956441528a
# ╠═32a76355-0e85-48e8-8a98-dc209398e34c
# ╠═8d61f6cc-784e-4cc8-9e21-070da8c84a1b
# ╠═60719c7d-8e19-4cdf-aeda-3ff83f446922
# ╠═15460682-a64a-41d6-934c-c66f2fb40f75
# ╠═de169125-f359-408c-83af-5aebb5fcf2cc
