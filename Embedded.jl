### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 7ca0bea0-d4e6-11ee-0acc-639605e39050
begin
	using Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	
	include("./src/Mandel.jl")
	using PlutoUI, .Mandel

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

# ╔═╡ 26ed5acd-6513-4406-ad6b-e4dda954fb0c
julia(f; c=-0.1145 + 0.7428im)

# ╔═╡ fe15ced3-6678-46d5-9b8f-7496e1707a3a
julia(f; c=-0.1145 + 0.7428im, coloring_alg=:projective_convergence)

# ╔═╡ 13d8745b-fb7d-481f-8f5d-67b72c806da2
md"""
You can choose the initial center, initial diameter, and the parameter of the Julia set by choosing `center`, `diam`, and `c`. 

The variable `coloring_alg` chooses the coloring algorithm for the plot. The options are:
- `:escape_time` 
- `:projective_convergence`.
"""

# ╔═╡ d25169ce-0698-4af3-b402-1c89a1f0655b
md"""
You can also connect the tracker of the `mandel` plot to the parameter of the `julia` plot by using the pattern below. (This interface will be improved soon.)
"""

# ╔═╡ 1d93dadf-964b-4f20-904b-94e26ec59918
@bind point mandel(f; center=-0.5)

# ╔═╡ c80bcd93-4721-4507-9c17-00b8972090ef
julia(f; c=point)

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
# ╠═26ed5acd-6513-4406-ad6b-e4dda954fb0c
# ╠═fe15ced3-6678-46d5-9b8f-7496e1707a3a
# ╟─13d8745b-fb7d-481f-8f5d-67b72c806da2
# ╟─d25169ce-0698-4af3-b402-1c89a1f0655b
# ╠═1d93dadf-964b-4f20-904b-94e26ec59918
# ╠═c80bcd93-4721-4507-9c17-00b8972090ef
# ╟─66d4ccf1-ac02-42a3-b02c-6992efb8a84b
# ╠═7ca0bea0-d4e6-11ee-0acc-639605e39050
