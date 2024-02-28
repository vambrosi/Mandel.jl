### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 7ca0bea0-d4e6-11ee-0acc-639605e39050
begin
	using Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	
	include("./src/Mandel.jl")
	using .Mandel

	mandel = Mandel.mandel
	julia = Mandel.julia
end

# ╔═╡ fe433fb7-db21-4f24-b786-246dbd6ad280
md"""
# Embedded.jl

Example file to show how to use the embedded user interface. Copy this file and change the functions to test it out. The command `git pull` will automatically update functionalities on the copy.
"""

# ╔═╡ d9d42585-e924-4175-9fe1-d7aa30baace9
md"""
First, define a 1-parameter family of functions to iterate as the following example:
"""

# ╔═╡ 6ea317f3-3332-4485-a62a-f74d3a5eae01
f(z, c) = z^2 + c

# ╔═╡ 23b6ed24-13e1-41ee-98a7-66bd0e2a8a4f
md"""
Then, call `julia(f)` or `mandel(f)` to plot the respective sets.
"""

# ╔═╡ 1d93dadf-964b-4f20-904b-94e26ec59918
mandel(f; center=-0.5)

# ╔═╡ c80bcd93-4721-4507-9c17-00b8972090ef
julia(f; c=-0.12711+0.75706im)

# ╔═╡ e7e56282-cc36-498a-8bd6-82bd1fa81cbd
julia(f; c=-0.12711+0.75706im, coloring_alg=:projective_convergence)

# ╔═╡ 13d8745b-fb7d-481f-8f5d-67b72c806da2
md"""
You can choose the initial center, initial diameter, and the parameter of the Julia set by choosing `center`, `diam`, and `c`. 

The variable `coloring_alg` chooses the coloring algorithm for the plot. The options are:
- `:escape_time` 
- `:projective_convergence`.
"""

# ╔═╡ 66d4ccf1-ac02-42a3-b02c-6992efb8a84b
md"""
## Instantiating Environment

This section imports the relevant functions.
"""

# ╔═╡ Cell order:
# ╟─fe433fb7-db21-4f24-b786-246dbd6ad280
# ╟─d9d42585-e924-4175-9fe1-d7aa30baace9
# ╠═6ea317f3-3332-4485-a62a-f74d3a5eae01
# ╟─23b6ed24-13e1-41ee-98a7-66bd0e2a8a4f
# ╠═1d93dadf-964b-4f20-904b-94e26ec59918
# ╠═c80bcd93-4721-4507-9c17-00b8972090ef
# ╠═e7e56282-cc36-498a-8bd6-82bd1fa81cbd
# ╟─13d8745b-fb7d-481f-8f5d-67b72c806da2
# ╟─66d4ccf1-ac02-42a3-b02c-6992efb8a84b
# ╠═7ca0bea0-d4e6-11ee-0acc-639605e39050
