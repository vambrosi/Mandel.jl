### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 4a085b76-da75-11ee-2b50-07309082fb46
begin
    # The lines below activate the environment defined in Project.toml
    # It will install the necessary packages the first time it runs.
    using Pkg
    Pkg.activate(@__DIR__)

    # Includes the file with definitions and imports the relevant modules
    include("./src/MandelMakie.jl")
    using .MandelMakie

    # The next lines are necessary because Pluto does not expose the exported variables of a local package. You can skip those lines if you are loading this package into the REPL or Jupyter notebook.
    Viewer = MandelMakie.Viewer
    Viewer3D = MandelMakie.Viewer3D
    Julia3D = MandelMakie.Julia3D
    Fatou3D = MandelMakie.Fatou3D
    set_parameter! = MandelMakie.set_parameter!
    get_attractors = MandelMakie.get_attractors
end

# ╔═╡ 1ba72fef-d5c4-4da4-9c17-4ba0096bf968
md"""
# Dylan's families
"""

# ╔═╡ fb6c0e42-b256-47aa-a496-aa5e3e4a57cc
md"""
## The quadratic family
"""

# ╔═╡ f49a6cef-86c9-4a47-a9ff-30588b3de27c
q(z, c) = z^2 + c

# ╔═╡ 6b1c8547-06de-4bb7-9062-0826ad04b92e
md"""
Then, call `Viewer` or `Viewer3D` to open a window where you see the Mandelbrot and Julia sets associated with that family of functions. (`Viewer3D` should only be used with rational functions.)
"""

# ╔═╡ 70f82c8f-5b1b-46c7-87b3-5fad91ba1094
# ╠═╡ disabled = true
#=╠═╡
Viewer(q, mandel_center = -0.5,convergence_criterion=:near_attractor)
  ╠═╡ =#

# ╔═╡ 3556d294-8f27-4c1c-8568-d2c548869e79
# ╠═╡ disabled = true
#=╠═╡
v = Viewer((z,λ) -> z^2 + λ*z, crit=λ -> -λ/2, convergence_criterion=:near_attractor)
  ╠═╡ =#

# ╔═╡ 781ae5e9-545e-4c6b-9285-86f01a91e8ff
# ╠═╡ disabled = true
#=╠═╡
MandelMakie.get_attractors(v)
  ╠═╡ =#

# ╔═╡ 5cdbb1b3-e783-4a1f-839d-b15c566bc1f3
# ╠═╡ disabled = true
#=╠═╡
Viewer3D(q)
  ╠═╡ =#

# ╔═╡ 63627492-a386-4137-ad07-05b687df03c5
md"""
## The family Per``_2``
"""

# ╔═╡ dc859310-0b21-4ede-a740-2b285f52a3df
p2(z, c) = (z^2 - c) / (z^2 - 1)

# ╔═╡ 8848c897-3fb7-4a7d-aa5f-d6b34438decb
# ╠═╡ disabled = true
#=╠═╡
v = Viewer(p2, mandel_center = -0.5, mandel_diameter = 6.0, julia_diameter = 6.0, convergence_criterion=:near_attractor)
  ╠═╡ =#

# ╔═╡ ac6a0d56-8990-45a8-90e0-4b876024c1c1
# ╠═╡ disabled = true
#=╠═╡
Viewer(
	p2,
	mandel_center = -0.5,
	mandel_diameter = 6.0,
	julia_diameter = 6.0,
	coloring_method=:mod_period
)
  ╠═╡ =#

# ╔═╡ 130ce106-d086-4c4e-af9d-9369208e2562
# ╠═╡ disabled = true
#=╠═╡
Viewer3D(p2)
  ╠═╡ =#

# ╔═╡ fa18fc4e-c7a7-4268-ab77-368a277d2ae6
md"""
# The ``Per_3(0)`` family
"""

# ╔═╡ 363f547e-715c-4b5e-b657-5c41f8e674ca
per3(z, c) = (z^2 - 1 - c + c^3) / (z^2 - c^2)

# ╔═╡ dd5a428d-c4af-4121-8905-902ca0ec1310
# ╠═╡ disabled = true
#=╠═╡
Viewer(per3, mandel_diameter = 6.0, coloring_method=:convergence_time)
  ╠═╡ =#

# ╔═╡ 744228c4-23f4-409a-805a-faa672e6fa9b
# ╠═╡ disabled = true
#=╠═╡
Viewer(per3, mandel_diameter = 6.0, coloring_method=:mod_period)
  ╠═╡ =#

# ╔═╡ ef15c50e-4158-4716-994b-c8885d9f81af
# ╠═╡ disabled = true
#=╠═╡
Viewer3D(per3)
  ╠═╡ =#

# ╔═╡ a9c743b4-7214-4b4d-b229-39aea38215bb
md"""
# The ``\mathrm{Per}_4(0)`` family
"""

# ╔═╡ b72a25ca-92a4-4950-88a9-e027335fbe9a
md"""
This is the parametrization from Milnor's 93 paper *Quadratic Rational Maps*. The marked cycle is
```math
0 \overset{2}{\longmapsto} \infty \longmapsto 1 \longmapsto \rho \longmapsto 0
```
and the free (active) critical point is at
```math
c = \frac{4\rho^2-2\rho}{-1+\rho+\rho^2}.
```
Note that $c$ is a quadratic function of $\rho$.
"""

# ╔═╡ 93eb1b62-e7ca-4969-8a1e-5902d509b8a4
per4(z, ρ) = (z - ρ) * (z - (2ρ - 1) / (ρ - 1)) / z^2

# ╔═╡ bea223da-13a5-4662-9e1f-ed6b27cf6d54
per4crit(ρ) = (4ρ^2 - 2ρ) / (-1 + ρ + ρ^2)

# ╔═╡ 2962904d-2c4d-4fe7-9c59-f735ab953bfa
# ╠═╡ disabled = true
#=╠═╡
Viewer(
	per4,
	crit = per4crit,
	coloring_method=:convergence_time,
	projective_metric=true
)
  ╠═╡ =#

# ╔═╡ 43af11b5-9593-4950-b853-8fa8869e6ace
# ╠═╡ disabled = true
#=╠═╡
Viewer(
	per4,
	crit = per4crit,
	coloring_method=:mod_period,
	projective_metric=true
)
  ╠═╡ =#

# ╔═╡ 9f9585b3-011c-4ae9-b5ce-069ecef3a23e
md"""
## Census maps
"""

# ╔═╡ f9c1ce04-81de-4458-8bf1-8225956b1011
f49(z) = (z / 3 - 1)^3 / (z - 1)

# ╔═╡ 1991d152-f45b-4b00-ab29-32f6d80c8997
# ╠═╡ disabled = true
#=╠═╡
Fatou3D(f49, show_critical_points = true)
  ╠═╡ =#

# ╔═╡ f966a1a4-8b1a-4f84-8452-50ac7b0024e9
f51(z) = 4 / 27 * (z - 1)^3 / z

# ╔═╡ 22b75feb-5c5e-42a4-aecc-b2f5ef4a938b
# ╠═╡ disabled = true
#=╠═╡
Julia3D(f51)
  ╠═╡ =#

# ╔═╡ 79cbbb8b-8add-4d95-87f2-a9d825494636
md"""
## The Devaney (2,2) family
"""

# ╔═╡ fb1f9c3a-ec39-4d6f-acb3-1ddfd56d614b
devaney22(z, λ) = z^2 + λ / z^2

# ╔═╡ 36714bb2-3622-47eb-a7b5-f4c99145bd59
devaney22crit(λ) = λ^(1 / 4)

# ╔═╡ 1c99e13e-e034-4d55-91be-ec22f1125a3d
# ╠═╡ disabled = true
#=╠═╡
Viewer(devaney22, crit = devaney22crit, mandel_diameter = 1.0)
  ╠═╡ =#

# ╔═╡ 29c2e25e-77f0-4d9f-94b0-163f8e58c7f8
md"""
## Devaney (2,1): A cubic symmetry locus
"""

# ╔═╡ 0e78c206-3a75-419c-80d6-5ef4a8af6dc4
md"""
Cubic maps with $f(\infty) = \infty$, $f'(\infty) = 0$ and the symmetry $f(\omega z) = \omega^2 f(z)$ where $\omega=\omega_3 = e^{2\pi i/3}$.
"""

# ╔═╡ 6cecbde5-42bc-4d74-8f68-2b458912585a
devaney21(z,λ) = z^2 + λ/z

# ╔═╡ 1388f3fa-c79e-4e4e-9455-b789fdb51378
devaney21crit(λ) = (λ/2)^(1/3)

# ╔═╡ 3c9b33d6-19a5-4cea-a787-cb143c9c82ee
# ╠═╡ disabled = true
#=╠═╡
Viewer(devaney21, crit = devaney21crit, mandel_diameter=1.5,c=-16/27)
  ╠═╡ =#

# ╔═╡ 60c35275-4216-40c1-91f0-ae35d0e8b31a
md"""
Special values: for $\lambda_1 = -16/27$, the Julia set is a Sierpiński gasket (with preperiodic free crit pt). For $\lambda_2=2/27$, the Julia set is the double of a Sierpiński gasket, this time hyperbolic.
"""

# ╔═╡ c70accb2-d5ba-405d-96cf-745dd8daba06
gasket1(z) = devaney21(z,-16/27)

# ╔═╡ bd6f8d54-1fde-4b73-a9f2-196e1df43eb2
gasket2(z) = devaney21(z, 2/27)

# ╔═╡ 35a3f7e6-18df-4733-96ec-e746b3adc38e
# ╠═╡ disabled = true
#=╠═╡
Julia3D(gasket2)
  ╠═╡ =#

# ╔═╡ dc85ce65-38a8-4990-86f5-ed6894b96695
md"""
## Instantiating Package
"""

# ╔═╡ d04798ec-c7c7-4868-b761-f9bb65e9a269
#=╠═╡
MandelMakie.get_attractors(v)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─1ba72fef-d5c4-4da4-9c17-4ba0096bf968
# ╟─fb6c0e42-b256-47aa-a496-aa5e3e4a57cc
# ╠═f49a6cef-86c9-4a47-a9ff-30588b3de27c
# ╟─6b1c8547-06de-4bb7-9062-0826ad04b92e
# ╠═70f82c8f-5b1b-46c7-87b3-5fad91ba1094
# ╠═3556d294-8f27-4c1c-8568-d2c548869e79
# ╠═781ae5e9-545e-4c6b-9285-86f01a91e8ff
# ╠═5cdbb1b3-e783-4a1f-839d-b15c566bc1f3
# ╟─63627492-a386-4137-ad07-05b687df03c5
# ╠═dc859310-0b21-4ede-a740-2b285f52a3df
# ╠═8848c897-3fb7-4a7d-aa5f-d6b34438decb
# ╠═ac6a0d56-8990-45a8-90e0-4b876024c1c1
# ╠═130ce106-d086-4c4e-af9d-9369208e2562
# ╟─fa18fc4e-c7a7-4268-ab77-368a277d2ae6
# ╠═363f547e-715c-4b5e-b657-5c41f8e674ca
# ╠═dd5a428d-c4af-4121-8905-902ca0ec1310
# ╠═744228c4-23f4-409a-805a-faa672e6fa9b
# ╠═ef15c50e-4158-4716-994b-c8885d9f81af
# ╟─a9c743b4-7214-4b4d-b229-39aea38215bb
# ╟─b72a25ca-92a4-4950-88a9-e027335fbe9a
# ╠═93eb1b62-e7ca-4969-8a1e-5902d509b8a4
# ╠═bea223da-13a5-4662-9e1f-ed6b27cf6d54
# ╠═2962904d-2c4d-4fe7-9c59-f735ab953bfa
# ╠═43af11b5-9593-4950-b853-8fa8869e6ace
# ╠═9f9585b3-011c-4ae9-b5ce-069ecef3a23e
# ╠═f9c1ce04-81de-4458-8bf1-8225956b1011
# ╠═1991d152-f45b-4b00-ab29-32f6d80c8997
# ╠═f966a1a4-8b1a-4f84-8452-50ac7b0024e9
# ╠═22b75feb-5c5e-42a4-aecc-b2f5ef4a938b
# ╟─79cbbb8b-8add-4d95-87f2-a9d825494636
# ╠═fb1f9c3a-ec39-4d6f-acb3-1ddfd56d614b
# ╠═36714bb2-3622-47eb-a7b5-f4c99145bd59
# ╠═1c99e13e-e034-4d55-91be-ec22f1125a3d
# ╠═29c2e25e-77f0-4d9f-94b0-163f8e58c7f8
# ╟─0e78c206-3a75-419c-80d6-5ef4a8af6dc4
# ╠═6cecbde5-42bc-4d74-8f68-2b458912585a
# ╠═1388f3fa-c79e-4e4e-9455-b789fdb51378
# ╠═3c9b33d6-19a5-4cea-a787-cb143c9c82ee
# ╟─60c35275-4216-40c1-91f0-ae35d0e8b31a
# ╠═c70accb2-d5ba-405d-96cf-745dd8daba06
# ╠═bd6f8d54-1fde-4b73-a9f2-196e1df43eb2
# ╠═35a3f7e6-18df-4733-96ec-e746b3adc38e
# ╟─dc85ce65-38a8-4990-86f5-ed6894b96695
# ╟─4a085b76-da75-11ee-2b50-07309082fb46
# ╠═d04798ec-c7c7-4868-b761-f9bb65e9a269
