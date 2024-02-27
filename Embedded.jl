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

# ╔═╡ 6ea317f3-3332-4485-a62a-f74d3a5eae01
f(z, c) = z^2 + c

# ╔═╡ c3bc70c0-fc7d-408c-a947-0c9a55e04537
g(z, c) = z^2 + c/z^2

# ╔═╡ db2c2bd5-3e0f-4e26-84f9-d1de720cb10c
h(z, c) = (z^2 - c) / (z^2 - 1)

# ╔═╡ 66d4ccf1-ac02-42a3-b02c-6992efb8a84b
md"""
## Definitions
"""

# ╔═╡ 4b35056f-00aa-459d-93f0-dd8422288187
function continuous_iter(n::Integer, z::Complex)
	return mod((n + 1.0 - log2(log(abs(z)))) / 64.0, 1.0)
end

# ╔═╡ ea71085e-4f4e-4ae7-a492-3852bf266d14
function escape_time(f::Function, z::Complex, c::Complex)
	for n in 1:512
		z = f(z, c)
		if abs2(z) > 1000.0
			return continuous_iter(n, z)
		end
	end

	return 0.5
end

# ╔═╡ 498da509-562e-44cf-aaa0-77563d85e525
to_range(raw_range, n) = range(minimum(raw_range), maximum(raw_range); length=n)

# ╔═╡ 675bd95b-b86a-4479-81f2-b3d3e79f4dc6
function mandel_slice!(grid, f, ix, x, yrange)
	@inbounds for (iy, y) in enumerate(yrange)
		grid[iy, ix] = escape_time(f, 0.0im, Complex(x,y))
	end
end

# ╔═╡ 1861ac13-6087-4938-9705-f23ce940d7b4
function julia_slice!(grid, f, ix, x, c, yrange)
	@inbounds for (iy, y) in enumerate(yrange)
		grid[iy, ix] = escape_time(f, Complex(x,y), c)
	end
end

# ╔═╡ 4d34509f-b36b-42e4-9e68-8f2540eddeae
function grid(f, c, xrange, yrange, is_mandel)
	grid = Matrix{Float64}(undef, length(xrange), length(yrange))
	futures = Vector{Task}(undef, length(xrange))
	if is_mandel
		@inbounds for (ix, x) in enumerate(xrange)
			futures[ix] = Threads.@spawn mandel_slice!(grid, f, ix, x, yrange)
		end
	else
		@inbounds for (ix, x) in enumerate(xrange)
			futures[ix] = Threads.@spawn julia_slice!(grid, f, ix, x, c, yrange)
		end
	end
	wait.(futures)
	return grid
end

# ╔═╡ 199f62a9-2b3c-43c5-a4f9-46d82eb38a9c
function view(
	f::Function, center::ComplexF64, diam::Real, c::ComplexF64, is_mandel::Bool
)
	n = 650
	xlo, xhi = real(center) - diam / 2, real(center) + diam / 2
	ylo, yhi = imag(center) - diam / 2, imag(center) + diam / 2
	
	x = to_range([xlo, xhi], n)
	y = to_range([ylo, yhi], n)
	z = grid(f, c, x, y, is_mandel)
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
				range = [xlo, xhi],
				showgrid = false,
				showticklabels = false,
				zeroline = false,
			),
			yaxis = attr(;
				range = [ylo, yhi],
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
		),
		config = PlotConfig(responsive=false, modeBarButtonsToRemove=["autoscale"])
	)
	add_plotly_listener!(p, "plotly_relayout", htl_js("async (e) => {
		let xrange = PLOT._fullLayout.xaxis.range
		let yrange = PLOT._fullLayout.yaxis.range
		if (PLOT.updated_data == undefined) {return}
		let update = await PLOT.updated_data({x: xrange, y: yrange})
		Plotly.restyle(PLOT, update, 0)
	}"))
	new_grid(d) = let
		x = range(d["x"]..., n) |> collect
		y = range(d["y"]..., n) |> collect
		z = grid(f, c, x, y, is_mandel) |> PlutoPlotly._preprocess
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
		plt.updated_data = $(AbstractPlutoDingetjes.Display.with_js_link(new_grid))
	</script>
	""")
end

# ╔═╡ 8ece44ee-a8c9-4a9a-a7cc-cbcf0fbe741b
julia(f; center=0.0im, diam=4.0, c=0.0im) = view(
	f, convert(ComplexF64, center), diam, convert(ComplexF64, c), false
)

# ╔═╡ e7e56282-cc36-498a-8bd6-82bd1fa81cbd
julia(f; c=1.0im)

# ╔═╡ ac46ec26-c395-4aa0-bd6c-1e66fefea8d3
julia(g; c=-0.05)

# ╔═╡ 37edc9a2-e864-4963-959b-d7760184d0ba
mandel(f; center=0.0im, diam=4.0, c=0.0im) = view(
	f, convert(ComplexF64, center), diam, convert(ComplexF64, c), true
)

# ╔═╡ 1d93dadf-964b-4f20-904b-94e26ec59918
mandel(f; center=-0.5)

# ╔═╡ 3bdde31c-9fac-4475-824e-09371fba5df6
mandel(h; diam=6.0)

# ╔═╡ Cell order:
# ╟─fe433fb7-db21-4f24-b786-246dbd6ad280
# ╠═6ea317f3-3332-4485-a62a-f74d3a5eae01
# ╠═1d93dadf-964b-4f20-904b-94e26ec59918
# ╠═e7e56282-cc36-498a-8bd6-82bd1fa81cbd
# ╠═c3bc70c0-fc7d-408c-a947-0c9a55e04537
# ╠═ac46ec26-c395-4aa0-bd6c-1e66fefea8d3
# ╠═db2c2bd5-3e0f-4e26-84f9-d1de720cb10c
# ╠═3bdde31c-9fac-4475-824e-09371fba5df6
# ╟─66d4ccf1-ac02-42a3-b02c-6992efb8a84b
# ╟─7ca0bea0-d4e6-11ee-0acc-639605e39050
# ╟─4b35056f-00aa-459d-93f0-dd8422288187
# ╟─ea71085e-4f4e-4ae7-a492-3852bf266d14
# ╟─498da509-562e-44cf-aaa0-77563d85e525
# ╟─4d34509f-b36b-42e4-9e68-8f2540eddeae
# ╟─675bd95b-b86a-4479-81f2-b3d3e79f4dc6
# ╟─1861ac13-6087-4938-9705-f23ce940d7b4
# ╟─199f62a9-2b3c-43c5-a4f9-46d82eb38a9c
# ╟─8ece44ee-a8c9-4a9a-a7cc-cbcf0fbe741b
# ╟─37edc9a2-e864-4963-959b-d7760184d0ba
