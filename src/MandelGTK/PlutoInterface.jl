### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 3a8b34f3-e4b3-4fd3-b94c-900b62e7dcdc
begin
    using Pkg
    Pkg.activate(Base.current_project())
	Pkg.instantiate()

    include("./src/MandelGTK.jl")
    using .MandelGTK

    viewer = MandelGTK.viewer
end

# ╔═╡ 1cb7de2e-39d0-4e30-a668-f01281b35682
f(z,c) = z^2 + c

# ╔═╡ 174b726c-5c63-423b-b07f-21c8fcebf01a
viewer(f, 1.0im)

# ╔═╡ Cell order:
# ╠═3a8b34f3-e4b3-4fd3-b94c-900b62e7dcdc
# ╠═1cb7de2e-39d0-4e30-a668-f01281b35682
# ╠═174b726c-5c63-423b-b07f-21c8fcebf01a
