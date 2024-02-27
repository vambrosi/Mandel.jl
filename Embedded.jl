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
end

# ╔═╡ 675bd95b-b86a-4479-81f2-b3d3e79f4dc6
begin
to_range(raw_range, n) = range(minimum(raw_range), maximum(raw_range); length=n)
function mandel(z:: Complex) 
    c = z 
    maxiter = 512
    for n = 1:maxiter 
        if abs2(z) > 100.0 
            return 2*(n-1) % 256
        end 
        z = z^2 + c 
    end 
    return 128
end
function mandel_slice!(grid, ix, x, yrange)
	@inbounds for (iy, y) in enumerate(yrange)
		grid[iy, ix] = mandel(Complex(x,y))
	end
end
function mandel_grid(xrange, yrange)
	grid = Matrix{UInt8}(undef, length(xrange), length(yrange))
	futures = Vector{Task}(undef, length(xrange))
	@inbounds for (ix, x) in enumerate(xrange)
		futures[ix] = Threads.@spawn mandel_slice!(grid, ix, x, yrange)
	end
	wait.(futures)
	grid
end
end

# ╔═╡ 199f62a9-2b3c-43c5-a4f9-46d82eb38a9c
let
	n = 500
	x = to_range([-2.5, 1.5], n)
	y = to_range([-2, 2], n)
	z = mandel_grid(x,y)
	data = heatmap(;x,y,z,
		colorscale=PlutoPlotly.colors.cyclical[:twilight],
		zmin = 0, zmax = 255,
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
				tickmode = "array",
				tickvals = [],
			),
			yaxis = attr(;
				range = [-2,2],
				scaleanchor = "x",
				scaleratio = 1,
				showgrid = false,
				showticklabels = false,
				zeroline = false,
				tickmode = "array",
				tickvals = [],
			),
			showlegend = false,
			height = 550,
			width = 550,
		)
	)
	add_plotly_listener!(p, "plotly_relayout", htl_js("async (e) => {
		let xrange = PLOT._fullLayout.xaxis.range
		let yrange = PLOT._fullLayout.yaxis.range
		if (PLOT.updated_data == undefined) {return}
		let update = await PLOT.updated_data({x: xrange, y: yrange})
		Plotly.restyle(PLOT, update)
	}"))
	f(d) = let
		x = range(d["x"]..., n) |> collect
		y = range(d["y"]..., n) |> collect
		z = mandel_grid(x,y) |> PlutoPlotly._preprocess
		out = (;
			x = [x],
			y = [y],
			z = [z],
		)
		out
	end
	@htl("""
	$p
	<script>
		const plt = currentScript.closest('pluto-cell').querySelector('.js-plotly-plot')
		plt.updated_data = $(AbstractPlutoDingetjes.Display.with_js_link(f))
	</script>
	""")
end


# ╔═╡ Cell order:
# ╠═7ca0bea0-d4e6-11ee-0acc-639605e39050
# ╠═675bd95b-b86a-4479-81f2-b3d3e79f4dc6
# ╠═199f62a9-2b3c-43c5-a4f9-46d82eb38a9c
