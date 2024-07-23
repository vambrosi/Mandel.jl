### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 4a085b76-da75-11ee-2b50-07309082fb46
begin
	# The lines below activate the environment defined in Project.toml
	# It will install the necessary packages the first time it runs.
	using Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()

	# Includes the file with definitions and imports the relevant modules
	include("./src/MandelMakie.jl")
	using .MandelMakie

	# The next lines are necessary because Pluto does not expose the exported variables of a local package. You can skip those lines if you are loading this package into the REPL or Jupyter notebook.
	Viewer = MandelMakie.Viewer
	Viewer3D = MandelMakie.Viewer3D
	Julia3D = MandelMakie.Julia3D
	Fatou3D = MandelMakie.Fatou3D
	to_complex = MandelMakie.to_complex
	set_parameter! = MandelMakie.set_parameter!
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
Viewer(q, mandel_center=-0.5)
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
p2(z,c) = (z^2-c)/(z^2-1)

# ╔═╡ 8848c897-3fb7-4a7d-aa5f-d6b34438decb
# ╠═╡ disabled = true
#=╠═╡
Viewer(p2, mandel_center=-1)
  ╠═╡ =#

# ╔═╡ 130ce106-d086-4c4e-af9d-9369208e2562
Viewer3D(p2)

# ╔═╡ 79cbbb8b-8add-4d95-87f2-a9d825494636
md"""
## The Devaney (2,2) family
"""

# ╔═╡ fb1f9c3a-ec39-4d6f-acb3-1ddfd56d614b
devaney22(z,λ) = z^2 + λ/z^2

# ╔═╡ 36714bb2-3622-47eb-a7b5-f4c99145bd59
devaneycrit(λ) = λ^(1/4)

# ╔═╡ 1c99e13e-e034-4d55-91be-ec22f1125a3d
# ╠═╡ disabled = true
#=╠═╡
Viewer(devaney22,crit=(λ->λ^(1/4)))
  ╠═╡ =#

# ╔═╡ fa18fc4e-c7a7-4268-ab77-368a277d2ae6
md"""
# The ``Per_3(0)`` family
"""

# ╔═╡ 363f547e-715c-4b5e-b657-5c41f8e674ca
per3(z,c) = (z^2 - 1  - c + c^3)/(z^2-c^2)

# ╔═╡ dd5a428d-c4af-4121-8905-902ca0ec1310
# ╠═╡ disabled = true
#=╠═╡
Viewer(per3)
  ╠═╡ =#

# ╔═╡ ef15c50e-4158-4716-994b-c8885d9f81af
# ╠═╡ disabled = true
#=╠═╡
Viewer3D(per3)
  ╠═╡ =#

# ╔═╡ 9f9585b3-011c-4ae9-b5ce-069ecef3a23e
md"""
## Census manifolds
"""

# ╔═╡ f9c1ce04-81de-4458-8bf1-8225956b1011
f49(z) = (z/3-1)^3/(z-1)

# ╔═╡ 1991d152-f45b-4b00-ab29-32f6d80c8997
Fatou3D(f49, show_critical_points=true)

# ╔═╡ f966a1a4-8b1a-4f84-8452-50ac7b0024e9
f51(z) = 4/27 * (z-1)^3/z

# ╔═╡ 22b75feb-5c5e-42a4-aecc-b2f5ef4a938b
# ╠═╡ disabled = true
#=╠═╡
Julia3D(f51)
  ╠═╡ =#

# ╔═╡ dc85ce65-38a8-4990-86f5-ed6894b96695
md"""
## Instantiating Package
"""

# ╔═╡ Cell order:
# ╟─1ba72fef-d5c4-4da4-9c17-4ba0096bf968
# ╟─fb6c0e42-b256-47aa-a496-aa5e3e4a57cc
# ╠═f49a6cef-86c9-4a47-a9ff-30588b3de27c
# ╟─6b1c8547-06de-4bb7-9062-0826ad04b92e
# ╠═70f82c8f-5b1b-46c7-87b3-5fad91ba1094
# ╠═5cdbb1b3-e783-4a1f-839d-b15c566bc1f3
# ╠═63627492-a386-4137-ad07-05b687df03c5
# ╠═dc859310-0b21-4ede-a740-2b285f52a3df
# ╠═8848c897-3fb7-4a7d-aa5f-d6b34438decb
# ╠═130ce106-d086-4c4e-af9d-9369208e2562
# ╟─79cbbb8b-8add-4d95-87f2-a9d825494636
# ╠═fb1f9c3a-ec39-4d6f-acb3-1ddfd56d614b
# ╠═36714bb2-3622-47eb-a7b5-f4c99145bd59
# ╠═1c99e13e-e034-4d55-91be-ec22f1125a3d
# ╠═fa18fc4e-c7a7-4268-ab77-368a277d2ae6
# ╠═363f547e-715c-4b5e-b657-5c41f8e674ca
# ╠═dd5a428d-c4af-4121-8905-902ca0ec1310
# ╠═ef15c50e-4158-4716-994b-c8885d9f81af
# ╠═9f9585b3-011c-4ae9-b5ce-069ecef3a23e
# ╠═f9c1ce04-81de-4458-8bf1-8225956b1011
# ╠═1991d152-f45b-4b00-ab29-32f6d80c8997
# ╠═f966a1a4-8b1a-4f84-8452-50ac7b0024e9
# ╠═22b75feb-5c5e-42a4-aecc-b2f5ef4a938b
# ╟─dc85ce65-38a8-4990-86f5-ed6894b96695
# ╠═4a085b76-da75-11ee-2b50-07309082fb46
