module Mandel

export mandel, julia

using StaticArraysCore, LinearAlgebra, Symbolics, PlutoPlotly, HypertextLiteral, AbstractPlutoDingetjes
using PrecompileTools

# --------------------------------------------------------------------------- #
# Complex Projective Line
# --------------------------------------------------------------------------- #

const Point = MVector{2, ComplexF64}

function distance(pt1::Point, pt2::Point)
	return norm(pt1[1] * pt2[2] - pt1[2] * pt2[1])
end

function projectivize(f::Function)
	@variables u, v, c
	frac = simplify(f(u / v, c))

	n, d = Symbolics.arguments(Symbolics.value(frac))
	fu = eval(build_function(n, u, v, c, expression=Val{false}))
	fv = eval(build_function(d, u, v, c, expression=Val{false}))

	return (pt, param) ->
		normalize!(Point(fu(pt[1], pt[2], param), fv(pt[1], pt[2], param)))
end

# --------------------------------------------------------------------------- #
# Coloring Algorithms
# --------------------------------------------------------------------------- #

const DEFAULT_COLOR::Float64 = 0.5

function continuous_iter(n::Integer, z::Complex)
	return mod((n + 1.0 - log2(log(abs(z)))) / 64.0, 1.0)
end

function escape_time(f::Function, z::Complex, c::Complex)
    for n in 1:512
		z = f(z, c)
		if abs2(z) > 1000.0
			return continuous_iter(n, z)
		end
	end

	return DEFAULT_COLOR
end

function repeated_iterate(f::Function, pt::Point, c::Number, ε::Real, max_iter::Integer)
	slow = copy(pt)
	fast = copy(pt)

	iteration = 1
	while iteration <= max_iter
		slow = f(slow, c)
		fast = f(f(fast, c), c)

		if distance(fast, slow) <= ε
			break
		end

		iteration += 1
	end

	return iteration, slow, fast
end

function near_again(f::Function, pt1::Point, pt2::Point, c::Number, ε::Real, max_period::Integer)
	pt = copy(pt1)

	iteration = 1
	while iteration <= max_period
		pt = f(pt, c)

		if distance(pt, pt2) <= ε
			break
		end

		iteration += 1
	end

	return iteration, pt1, pt2
end

function iterate_until_close(f::Function, pt1::Point, pt2::Point, c::Number, ε::Real, max_iter::Integer)
	iteration = 0
	while iteration <= max_iter
		if distance(pt1, pt2) <= ε
			break
		end

		pt1 = f(pt1, c)
		pt2 = f(pt2, c)
		iteration += 1
	end

	return iteration
end

function projective_convergence(f_proj::Function, z::Complex, c::Complex)
	ε = 1e-4
	max_iter = 200
	pt = Point(z, 1)

	time_diff, slow, fast = repeated_iterate(f_proj, pt, c, ε, max_iter)
	time_diff > max_iter && return DEFAULT_COLOR

	period, slow, fast = near_again(f_proj, slow, fast, c, ε, time_diff)
	period > time_diff && return DEFAULT_COLOR

	preperiod = iterate_until_close(f_proj, pt, fast, c, ε, max_iter)
	preperiod > max_iter && return DEFAULT_COLOR

	return mod((preperiod / period) / 64.0, 1.0)
end

algorithms = Dict(:escape_time => escape_time, :projective_convergence => projective_convergence)

# --------------------------------------------------------------------------- #
# Grid calculations
# --------------------------------------------------------------------------- #

mutable struct Viewer
	pluto_code::HypertextLiteral.Result
	tracked_point::ComplexF64
end

to_range(raw_range, n) = range(minimum(raw_range), maximum(raw_range); length=n)

function grid(f, c, xrange, yrange, is_mandel, algorithm)
	grid = Matrix{Float64}(undef, length(xrange), length(yrange))
	futures = Vector{Task}(undef, length(xrange))

	if is_mandel
		@inbounds for (ix, x) in enumerate(xrange)
			futures[ix] = Threads.@spawn mandel_slice!(grid, f, ix, x, yrange, algorithm)
		end
	else
		@inbounds for (ix, x) in enumerate(xrange)
			futures[ix] = Threads.@spawn julia_slice!(grid, f, ix, x, c, yrange, algorithm)
		end
	end
	wait.(futures)
	return grid
end

function mandel_slice!(grid, f, ix, x, yrange, algorithm)
	@inbounds for (iy, y) in enumerate(yrange)
		grid[iy, ix] = algorithm(f, 0.0im, Complex(x,y))
	end
end

function julia_slice!(grid, f, ix, x, c, yrange, algorithm)
	@inbounds for (iy, y) in enumerate(yrange)
		grid[iy, ix] = algorithm(f, Complex(x,y), c)
	end
end

function view(
	f::Function, center::ComplexF64, diam::Real, c::ComplexF64, is_mandel::Bool, coloring_alg::Symbol
)
	n = 650
	xlo, xhi = real(center) - diam / 2, real(center) + diam / 2
	ylo, yhi = imag(center) - diam / 2, imag(center) + diam / 2

	algorithm = algorithms[coloring_alg]

	# Uses projective version of the function if the algorithm tests convergence on the sphere.
	g = coloring_alg == :projective_convergence ? projectivize(f) : f

	x = to_range([xlo, xhi], n)
	y = to_range([ylo, yhi], n)
	z = grid(g, c, x, y, is_mandel, algorithm)

	data = heatmap(;x,y,z,
		colorscale=PlutoPlotly.colors.cyclical[:twilight],
		zmin = 0.0, zmax = 1.0,
		showscale = false,
		name = is_mandel ? "Mandel" : "Julia",
		hoverinfo = "none",
	)

	points = scatter(x=[0.0], y=[0.0], mode="markers+lines", name="Tracker")

	p = plot(
		[data, points],
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
		config = PlotConfig(
			responsive=false,
			modeBarButtonsToRemove=["autoscale", "select2d", "lasso2d"],
		),
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
		z = grid(g, c, x, y, is_mandel, algorithm) |> PlutoPlotly._preprocess
		out = (;
			x = [x],
			y = [y],
			z = [z],
		)
		out
	end

	add_plotly_listener!(p, "plotly_click", """e => {
		let point_x = e.points[0].x
		let point_y = e.points[0].y
		Plotly.restyle(PLOT, {x: [[point_x]], y: [[point_y]]}, 1)

		let span = PLOT.closest("span")
		span.value = [point_x, point_y]
		span.dispatchEvent(new CustomEvent("input"))
	}""")

	orbit_fixed(d) = let
		v = orbit(f, complex(d["x"], d["y"]), c, is_mandel ? 1 : 1)
		(; x = real(v), y = imag(v))
	end

	return @htl("""
	<span>
		$p
		<script>
			const plt = currentScript.closest('pluto-cell').querySelector('.js-plotly-plot')
			plt.updated_data = $(AbstractPlutoDingetjes.Display.with_js_link(new_grid))
			plt.orbit = $(AbstractPlutoDingetjes.Display.with_js_link(orbit_fixed))

			const span = currentScript.parentElement
			span.value = [0, 0]
			span.dispatchEvent(new CustomEvent("input"))
		</script>
	</span>
	""")
end

function julia(f; center=0.0im, diam=4.0, c=0.0im, coloring_alg=:escape_time)
	if c === missing
		c = 0.0im
	end

	return view(
		f, convert(ComplexF64, center), convert(Float64, diam), convert(ComplexF64, complex(c...)), false, coloring_alg
	)
end

mandel(f; center=0.0im, diam=4.0, c=0.0im, coloring_alg=:escape_time) = view(
	f, convert(ComplexF64, center), convert(Float64, diam), convert(ComplexF64, c), true, coloring_alg
)

# --------------------------------------------------------------------------- #
# Other Functions
# --------------------------------------------------------------------------- #

function orbit(f::Function, z::Complex, c::Complex, n::Integer)
	v = Vector{ComplexF64}(undef, n)
	for i in eachindex(v)
		v[i] = z
		z = f(z, c)
	end
	return v
end

# --------------------------------------------------------------------------- #
# Precompiling
# --------------------------------------------------------------------------- #
@compile_workload begin
	f(z, c) = z^2 + c
	f_proj(pt, c) = normalize!(Point(pt[1]^2 + c * pt[2]^2, pt[2]^2))

	grid(f, 0.0im, range(-2.5, 1.5, 10), range(-2.0, 2.0, 10), true, escape_time)
	grid(f, -0.1145 + 0.7428im, range(-2.0, 2.0, 10), range(-2.0, 2.0, 10), false, escape_time)
	grid(f_proj, 0.0im, range(-2.5, 1.5, 10), range(-2.0, 2.0, 10), true, projective_convergence)
	grid(f_proj, -0.1145 + 0.7428im, range(-2.0, 2.0, 10), range(-2.0, 2.0, 10), false, projective_convergence)
end

end