### A Pluto.jl notebook ###
# v0.19.40

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
end

# ╔═╡ 1ba72fef-d5c4-4da4-9c17-4ba0096bf968
md"""
# Makie Examples

Below we have some examples of how to use the Makie user interface of `Mandel.jl`. This user interface can be run from the REPL, from Jupyter notebooks, or from Pluto. The last cell in this notebook imports all the necessary files and functions.

If you prefer, you can add that cell to your own Pluto notebook (or make a copy of the whole notebook), then follow the instructions here to test the GUI.
"""

# ╔═╡ fb6c0e42-b256-47aa-a496-aa5e3e4a57cc
md"""
First, you define a family of complex maps (don't put type annotations)
"""

# ╔═╡ f49a6cef-86c9-4a47-a9ff-30588b3de27c
f(z, λ) = z^2 + λ/z^2

# ╔═╡ 89426560-0874-45a2-a2dd-6b2a044faabb
crit(λ) = λ^(1/4)

# ╔═╡ 6b1c8547-06de-4bb7-9062-0826ad04b92e
md"""
Then, call `Viewer` or `Viewer3D` to open a window where you see the Mandelbrot and Julia sets associated with that family of functions. (`Viewer3D` should only be used with rational functions.)
"""

# ╔═╡ 70f82c8f-5b1b-46c7-87b3-5fad91ba1094
Viewer(f; crit=crit, mandel_diam=1.0)

# ╔═╡ 5cdbb1b3-e783-4a1f-839d-b15c566bc1f3
Viewer3D(f; crit=crit)

# ╔═╡ dc85ce65-38a8-4990-86f5-ed6894b96695
md"""
## Instantiating Package
"""

# ╔═╡ Cell order:
# ╟─1ba72fef-d5c4-4da4-9c17-4ba0096bf968
# ╟─fb6c0e42-b256-47aa-a496-aa5e3e4a57cc
# ╠═f49a6cef-86c9-4a47-a9ff-30588b3de27c
# ╠═89426560-0874-45a2-a2dd-6b2a044faabb
# ╟─6b1c8547-06de-4bb7-9062-0826ad04b92e
# ╠═70f82c8f-5b1b-46c7-87b3-5fad91ba1094
# ╠═5cdbb1b3-e783-4a1f-839d-b15c566bc1f3
# ╟─dc85ce65-38a8-4990-86f5-ed6894b96695
# ╠═4a085b76-da75-11ee-2b50-07309082fb46
