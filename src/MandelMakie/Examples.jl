### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 4a085b76-da75-11ee-2b50-07309082fb46
begin
    # The lines below activate the environment defined in Project.toml
    # It will install the necessary packages the first time it runs.
    using Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    # Includes the file with definitions and imports the relevant modules
    include("./src/MandelMakie.jl")
    using .MandelMakie

    # The next lines are necessary because Pluto does not expose the exported variables of a local package. You can skip those lines if you are loading this package into the REPL or Jupyter notebook.
    Viewer = MandelMakie.Viewer
    Viewer3D = MandelMakie.Viewer3D
    Julia3D = MandelMakie.Julia3D
    Fatou3D = MandelMakie.Fatou3D
    set_parameter! = MandelMakie.set_parameter!
    get_parameter = MandelMakie.get_parameter
    get_attractors = MandelMakie.get_attractors
end

# ╔═╡ 1ba72fef-d5c4-4da4-9c17-4ba0096bf968
md"""
# Makie Examples

In this notebook, we have some examples of how to use the Makie user interface of `Mandel.jl`. The last cell in this notebook imports all the necessary files and functions. If you prefer, you can add that cell to your own Pluto notebook (or make a copy of the notebook), then follow the instructions here to use the GUI.
"""

# ╔═╡ fb6c0e42-b256-47aa-a496-aa5e3e4a57cc
md"""
First, you define a family of complex maps
"""

# ╔═╡ f49a6cef-86c9-4a47-a9ff-30588b3de27c
f(z, c) = z^2 + c

# ╔═╡ 89426560-0874-45a2-a2dd-6b2a044faabb
crit(c) = c

# ╔═╡ 6b1c8547-06de-4bb7-9062-0826ad04b92e
md"""
Then, call `Viewer` or `Viewer3D` to open a window where you see the Mandelbrot and Julia sets associated with that family of functions. (`Viewer3D` should only be used with rational functions.)
"""

# ╔═╡ 70f82c8f-5b1b-46c7-87b3-5fad91ba1094
Viewer(f; crit = crit, mandel_diameter = 1.0)

# ╔═╡ 048320b4-12f2-4e7f-b267-8b15ad2485e0
# ╠═╡ disabled = true
#=╠═╡
Viewer3D(f; crit=crit, c=0.1)
  ╠═╡ =#

# ╔═╡ 7fa7d5e2-f8d3-41c1-9d52-08dc10150fcf
md"""
Another example would be the family ``g(z, \lambda)`` given by applying the Newton's method to a polynomial of degree 3 (with roots at ``1``, ``\lambda``, and ``-\lambda``).
"""

# ╔═╡ ea959fbb-23d7-484f-932d-a64fe7130046
g(z, λ) =
    z -
    (z - 1) * (z - λ) * (z + λ) /
    ((z - 1) * (z - λ) + (z - 1) * (z + λ) + (z + λ) * (z - λ))

# ╔═╡ 259aebbe-a6ac-4b86-8a27-f12e104c3092
# ╠═╡ disabled = true
#=╠═╡
Viewer3D(g; crit=1/3)
  ╠═╡ =#

# ╔═╡ 654fb826-65c8-4852-83ca-7e7f72351325
md"""
The subsections below give more details. Cells are disabled so they don't run at startup.
"""

# ╔═╡ fc92cb1c-0310-47ad-a74d-63363c5c6181
md"""
## Complex Plane Plots
"""

# ╔═╡ 879a14f8-207c-4510-b94d-0e7f6f392175
md"""
To get more details about the `Viewer` function just type `?Viewer` in an empty cell.
"""

# ╔═╡ 9ce5fbe4-2985-4665-8a30-682d93ac52ba
# ╠═╡ disabled = true
#=╠═╡
viewer = Viewer(g; crit=1/3, coloring_method=:convergence_time)
  ╠═╡ =#

# ╔═╡ 64508179-e449-4d9a-b669-a08cfd5b6df1
md"""
You can also see the attracting cycles for the map shown in a `Viewer`.

The cells below need to be manually updated to reflect the last changes in the viewer.
"""

# ╔═╡ b42b3753-33ba-4e25-a01b-b316dbd86787
# ╠═╡ disabled = true
#=╠═╡
get_parameter(viewer)
  ╠═╡ =#

# ╔═╡ f1a18d08-d81f-4667-b87a-8d19989b201b
# ╠═╡ disabled = true
#=╠═╡
get_attractors(viewer)
  ╠═╡ =#

# ╔═╡ 8a161348-be92-47cc-8e3c-62f74b8d96d2
md"""
## Complex Projective Line Plots
"""

# ╔═╡ 3ba77924-33c1-4d47-8367-379766a58a95
md"""
You can pick the Julia set parameter for `Viewer3D` by setting `c` as shown in the first example in the notebook (its default value is `0`). Or you can modify it by using the buttons on the viewer, or by calling `set_parameter!`, as shown below.
"""

# ╔═╡ 105f8b1b-de30-4f98-986f-9b667b87949b
# ╠═╡ disabled = true
#=╠═╡
viewer3D = Viewer3D((z,c) -> z^2 + c)
  ╠═╡ =#

# ╔═╡ 21dbad3b-86b2-43f7-912e-de5e794640d6
# ╠═╡ disabled = true
#=╠═╡
set_parameter!(viewer3D, -0.12256116687665 -0.74486176661974im)
  ╠═╡ =#

# ╔═╡ 6b30f798-21b3-44f0-a8a0-c8f06da21b33
md"""
### Julia Plots using Critical Orbits
"""

# ╔═╡ f0551a80-f208-470f-932a-7031e58b2d06
md"""
The function `Fatou3D(f, c)` plots the Julia set of the family `f` with parameter `c` and colors each Fatou component according to its limit set.
"""

# ╔═╡ ebbe4d7f-9ab0-47f2-85e4-c0de7cb141c0
# ╠═╡ disabled = true
#=╠═╡
Fatou3D(f, -0.160129312880546 + 0.13705188507903926im)
  ╠═╡ =#

# ╔═╡ 9e104891-fcb5-48fd-93a7-b2889381d2a6
# ╠═╡ disabled = true
#=╠═╡
Fatou3D(g, im)
  ╠═╡ =#

# ╔═╡ dc85ce65-38a8-4990-86f5-ed6894b96695
md"""
## Instantiating Package
"""

# ╔═╡ 0bbfff4b-7ee9-47f5-b19d-c97494438d4d
md"""
## Benchmarks
"""

# ╔═╡ 6e3d25a8-08a5-4e4d-a26c-6fce6623cef9
# ╠═╡ disabled = true
#=╠═╡
fatou = Fatou3D(g, im)
  ╠═╡ =#

# ╔═╡ cab73282-0e64-45eb-958b-92f03fb4c8ca
# ╠═╡ disabled = true
#=╠═╡
@time fatou.reset.clicks[] += 1
  ╠═╡ =#

# ╔═╡ acc7ac26-fe04-4da3-a329-aa8568bd72ef
# ╠═╡ disabled = true
#=╠═╡
julia = Julia3D(g, im)
  ╠═╡ =#

# ╔═╡ 949646ca-3cf9-46b9-9f0d-ced29cd2740f
# ╠═╡ disabled = true
#=╠═╡
@time julia.reset.clicks[] += 1
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─1ba72fef-d5c4-4da4-9c17-4ba0096bf968
# ╟─fb6c0e42-b256-47aa-a496-aa5e3e4a57cc
# ╠═f49a6cef-86c9-4a47-a9ff-30588b3de27c
# ╠═89426560-0874-45a2-a2dd-6b2a044faabb
# ╟─6b1c8547-06de-4bb7-9062-0826ad04b92e
# ╠═70f82c8f-5b1b-46c7-87b3-5fad91ba1094
# ╠═048320b4-12f2-4e7f-b267-8b15ad2485e0
# ╟─7fa7d5e2-f8d3-41c1-9d52-08dc10150fcf
# ╠═ea959fbb-23d7-484f-932d-a64fe7130046
# ╠═259aebbe-a6ac-4b86-8a27-f12e104c3092
# ╟─654fb826-65c8-4852-83ca-7e7f72351325
# ╟─fc92cb1c-0310-47ad-a74d-63363c5c6181
# ╟─879a14f8-207c-4510-b94d-0e7f6f392175
# ╠═9ce5fbe4-2985-4665-8a30-682d93ac52ba
# ╟─64508179-e449-4d9a-b669-a08cfd5b6df1
# ╠═b42b3753-33ba-4e25-a01b-b316dbd86787
# ╠═f1a18d08-d81f-4667-b87a-8d19989b201b
# ╟─8a161348-be92-47cc-8e3c-62f74b8d96d2
# ╟─3ba77924-33c1-4d47-8367-379766a58a95
# ╠═105f8b1b-de30-4f98-986f-9b667b87949b
# ╠═21dbad3b-86b2-43f7-912e-de5e794640d6
# ╟─6b30f798-21b3-44f0-a8a0-c8f06da21b33
# ╟─f0551a80-f208-470f-932a-7031e58b2d06
# ╠═ebbe4d7f-9ab0-47f2-85e4-c0de7cb141c0
# ╠═9e104891-fcb5-48fd-93a7-b2889381d2a6
# ╟─dc85ce65-38a8-4990-86f5-ed6894b96695
# ╠═4a085b76-da75-11ee-2b50-07309082fb46
# ╟─0bbfff4b-7ee9-47f5-b19d-c97494438d4d
# ╠═6e3d25a8-08a5-4e4d-a26c-6fce6623cef9
# ╠═cab73282-0e64-45eb-958b-92f03fb4c8ca
# ╠═acc7ac26-fe04-4da3-a329-aa8568bd72ef
# ╠═949646ca-3cf9-46b9-9f0d-ced29cd2740f
