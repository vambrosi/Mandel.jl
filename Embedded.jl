### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 7ca0bea0-d4e6-11ee-0acc-639605e39050
begin
	using Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	
	using PlutoPlotly, HypertextLiteral, AbstractPlutoDingetjes

	md"Instatiates project" 
end

# ╔═╡ fe433fb7-db21-4f24-b786-246dbd6ad280
md"""
# Embedded.jl

Viewer embedded in Pluto to explore dynamical systems.
"""

# ╔═╡ 66d4ccf1-ac02-42a3-b02c-6992efb8a84b
md"""
## Definitions
"""

# ╔═╡ 4b35056f-00aa-459d-93f0-dd8422288187
function continuous_iter(n::Integer, z::Complex)
	return mod((n + 1.0 - log2(log(abs(z)))) / 64.0, 1.0)
end

# ╔═╡ ea71085e-4f4e-4ae7-a492-3852bf266d14
function escape_time(z::Complex, c::Complex)
	for n in 1:512
		z = z^2 + c
		if abs2(z) > 1000.0
			return continuous_iter(n, z)
		end
	end

	return 0.5
end

# ╔═╡ 498da509-562e-44cf-aaa0-77563d85e525
to_range(raw_range, n) = range(minimum(raw_range), maximum(raw_range); length=n)

# ╔═╡ 675bd95b-b86a-4479-81f2-b3d3e79f4dc6
function mandel_slice!(grid, ix, x, yrange)
	@inbounds for (iy, y) in enumerate(yrange)
		grid[iy, ix] = escape_time(0.0im, Complex(x,y))
	end
end

# ╔═╡ 4d34509f-b36b-42e4-9e68-8f2540eddeae
function mandel_grid(xrange, yrange)
	grid = Matrix{Float64}(undef, length(xrange), length(yrange))
	futures = Vector{Task}(undef, length(xrange))
	@inbounds for (ix, x) in enumerate(xrange)
		futures[ix] = Threads.@spawn mandel_slice!(grid, ix, x, yrange)
	end
	wait.(futures)
	return grid
end

# ╔═╡ 199f62a9-2b3c-43c5-a4f9-46d82eb38a9c
function mandel()
	n = 650
	x = to_range([-2.5, 1.5], n)
	y = to_range([-2, 2], n)
	z = mandel_grid(x,y)
	data = heatmap(;x,y,z,
		colorscale=PlutoPlotly.colors.cyclical[:twilight],
		zmin = 0.0, zmax = 1.0,
		showscale = false,
	)
	p = plot(
		data, 
		Layout(
			margin = attr(t=0, b=0, l=0, r=0),
			xaxis = attr(;
				range = [-2.5, 1.5],
				showgrid = false,
				showticklabels = false,
				zeroline = false,
			),
			yaxis = attr(;
				range = [-2,2],
				scaleanchor = "x",
				scaleratio = 1,
				showgrid = false,
				showticklabels = false,
				zeroline = false,
			),
			showlegend = false,
			height = n,
			width = n,
			dragmode = "pan",
		)
	)
	add_plotly_listener!(p, "plotly_relayout", htl_js("async (e) => {
		let xrange = PLOT._fullLayout.xaxis.range
		let yrange = PLOT._fullLayout.yaxis.range
		if (PLOT.updated_data == undefined) {return}
		let update = await PLOT.updated_data({x: xrange, y: yrange})
		Plotly.restyle(PLOT, update, 0)
	}"))
	f(d) = let
		x = range(d["x"]..., n) |> collect
		y = range(d["y"]..., n) |> collect
		z = mandel_grid(x, y) |> PlutoPlotly._preprocess
		out = (;
			x = [x],
			y = [y],
			z = [z],
		)
		out
	end
	return @htl("""
	$p
	<script>
		const plt = currentScript.closest('pluto-cell').querySelector('.js-plotly-plot')
		plt.updated_data = $(AbstractPlutoDingetjes.Display.with_js_link(f))
	</script>
	""")
end

# ╔═╡ 1d93dadf-964b-4f20-904b-94e26ec59918
mandel()

# ╔═╡ 1861ac13-6087-4938-9705-f23ce940d7b4
function julia_slice!(grid, ix, x, c, yrange)
	@inbounds for (iy, y) in enumerate(yrange)
		grid[iy, ix] = escape_time(Complex(x,y), c)
	end
end

# ╔═╡ ee9a255d-152a-4727-9034-372a61947023
function julia_grid(xrange, yrange, c)
	grid = Matrix{Float64}(undef, length(xrange), length(yrange))
	futures = Vector{Task}(undef, length(xrange))
	@inbounds for (ix, x) in enumerate(xrange)
		futures[ix] = Threads.@spawn julia_slice!(grid, ix, x, c, yrange)
	end
	wait.(futures)
	return grid
end

# ╔═╡ 8ece44ee-a8c9-4a9a-a7cc-cbcf0fbe741b
function julia(c)
	n = 650
	x = to_range([-2, 2], n)
	y = to_range([-2, 2], n)
	z = julia_grid(x,y, convert(ComplexF64, c))
	data = heatmap(;x,y,z,
		colorscale=PlutoPlotly.colors.cyclical[:twilight],
		zmin = 0.0, zmax = 1.0,
		showscale = false,
	)
	p = plot(
		data, 
		Layout(
			margin = attr(t=0, b=0, l=0, r=0),
			xaxis = attr(;
				range = [-2, 2],
				showgrid = false,
				showticklabels = false,
				zeroline = false,
			),
			yaxis = attr(;
				range = [-2,2],
				scaleanchor = "x",
				scaleratio = 1,
				showgrid = false,
				showticklabels = false,
				zeroline = false,
			),
			showlegend = false,
			height = n,
			width = n,
			dragmode = "pan",
		)
	)
	add_plotly_listener!(p, "plotly_relayout", htl_js("async (e) => {
		let xrange = PLOT._fullLayout.xaxis.range
		let yrange = PLOT._fullLayout.yaxis.range
		if (PLOT.updated_data == undefined) {return}
		let update = await PLOT.updated_data({x: xrange, y: yrange})
		Plotly.restyle(PLOT, update, 0)
	}"))
	f(d) = let
		x = range(d["x"]..., n) |> collect
		y = range(d["y"]..., n) |> collect
		z = julia_grid(x, y, convert(ComplexF64, c)) |> PlutoPlotly._preprocess
		out = (;
			x = [x],
			y = [y],
			z = [z],
		)
		out
	end
	return @htl("""
	$p
	<script>
		const plt = currentScript.closest('pluto-cell').querySelector('.js-plotly-plot')
		plt.updated_data = $(AbstractPlutoDingetjes.Display.with_js_link(f))
	</script>
	""")
end

# ╔═╡ e7e56282-cc36-498a-8bd6-82bd1fa81cbd
julia(1.0im)

# ╔═╡ Cell order:
# ╟─fe433fb7-db21-4f24-b786-246dbd6ad280
# ╠═1d93dadf-964b-4f20-904b-94e26ec59918
# ╠═e7e56282-cc36-498a-8bd6-82bd1fa81cbd
# ╟─66d4ccf1-ac02-42a3-b02c-6992efb8a84b
# ╟─7ca0bea0-d4e6-11ee-0acc-639605e39050
# ╟─4b35056f-00aa-459d-93f0-dd8422288187
# ╟─ea71085e-4f4e-4ae7-a492-3852bf266d14
# ╟─498da509-562e-44cf-aaa0-77563d85e525
# ╟─4d34509f-b36b-42e4-9e68-8f2540eddeae
# ╟─675bd95b-b86a-4479-81f2-b3d3e79f4dc6
# ╟─199f62a9-2b3c-43c5-a4f9-46d82eb38a9c
# ╟─ee9a255d-152a-4727-9034-372a61947023
# ╟─1861ac13-6087-4938-9705-f23ce940d7b4
# ╟─8ece44ee-a8c9-4a9a-a7cc-cbcf0fbe741b
