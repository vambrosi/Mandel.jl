module MandelMakie

export Viewer, Viewer3D

using GLMakie, Symbolics, Parameters

import Base: show

include("Dynamics.jl")

# --------------------------------------------------------------------------------------- #
# Parsing Functions
# --------------------------------------------------------------------------------------- #

struct DynamicalSystem
	f::Function
	df_dz::Function
	df_dc::Function
	crit::Function

	function DynamicalSystem(f, crit)
        @variables z, c
		expr = f(z, c)

		df_dz_expr = expand_derivatives(Differential(z)(expr))
		df_dz = build_function(df_dz_expr, z, c, expression=Val{false})

		df_dc_expr = expand_derivatives(Differential(c)(expr))
		df_dc = build_function(df_dc_expr, z, c, expression=Val{false})

        if crit isa Function
            crit_function = crit
        else
            critical_point = convert(ComplexF64, crit)
            crit_function = (_) -> critical_point
        end

		return new(f, df_dz, df_dc, crit_function)
	end
end

function to_point_family(f)
	@variables u, v, c, d
	frac = simplify(f(u / v, c / d))
	value = Symbolics.value(frac)

	num, den = Symbolics.arguments(value)
	fu = build_function(num, u, v, c, d, expression=Val{false})
	fv = build_function(den, u, v, c, d, expression=Val{false})

	return (pt_z, pt_c) -> Point(fu(pt_z.u, pt_z.v, pt_c.u, pt_c.v), fv(pt_z.u, pt_z.v, pt_c.u, pt_c.v))
end

function to_point_map(f)
	@variables c, d
	frac = simplify(f(c / d))
	value = Symbolics.value(frac)

	num, den = Symbolics.arguments(value)
	fu = build_function(num, c, d, expression=Val{false})
	fv = build_function(den, c, d, expression=Val{false})

	return (pt) -> Point(fu(pt.u, pt.v), fv(pt.u, pt.v))
end

struct RationalMap
	f::Function
	crit::Function

	function RationalMap(f, crit)
		if crit isa Function
            crit_function = to_point_map(crit)
        else
            crit_function = (_) -> Point(crit, 1.0)
        end

		return new(to_point_family(f), crit_function)
	end
end

# --------------------------------------------------------------------------------------- #
# 2D Viewer for Complex Analytic Maps
# --------------------------------------------------------------------------------------- #

abstract type View end

mutable struct ViewerOptions
	esc_radius::Float64
	max_iter::Int
	orbit_len::Int
	crit_len::Int
	coloring_algorithm::Function
end

mutable struct ViewerState
	last_button::Symbol
end

function to_complex(view::View, point)
	x,  y = point[1], point[2]

	ppu = view.pixels / view.diameter
	a = ((x - 0.5) - view.pixels / 2) / ppu
	b = ((y - 0.5) - view.pixels / 2) / ppu

	return view.center + complex(a, b)
end

function to_image_coords(view::View, z::Number)
	ppu = view.pixels / view.diameter
	w = z - view.center

	x = real(w) * ppu + 0.5 + view.pixels / 2
	y = imag(w) * ppu + 0.5 + view.pixels / 2

	return [x, y]
end

function to_image_coords(view::View, zs::Vector{<:Number})
	ppu = view.pixels / view.diameter
	xs = Float64[]
	ys = Float64[]

	for z in zs
		w = z - view.center
		push!(xs, real(w) * ppu + 0.5 + view.pixels / 2)
		push!(ys, imag(w) * ppu + 0.5 + view.pixels / 2)
	end

	return xs, ys
end

function prepare!(view::View)
	hidedecorations!(view.axis)
	hidespines!(view.axis)
	deactivate_interaction!(view.axis, :rectanglezoom)
	deactivate_interaction!(view.axis, :scrollzoom)
	deactivate_interaction!(view.axis, :dragpan)

	plt = heatmap!(
		view.axis,
		view.array,
		colormap = :twilight,
		colorrange = (0.0, 1.0),
		inspectable = false,
		inspector_label = (self, i, p) -> let
			z = to_complex(view, p)
			"x: $(real(z))\ny: $(imag(z))"
		end,
	)

	plt.inspectable = view.inspectable

	points = lift(view.points) do zs
		xs, ys = to_image_coords(view, zs)
		return Point2f.(xs, ys)
	end

	lines!(
		view.axis,
		points,
		color = (:red, 0.5),
		inspectable = false,
	)

	scatter!(
		view.axis,
		points,
		color = (:red, 1.0),
		inspector_label = (self, i, p) -> let
			z = to_complex(view, p)
			"x: $(real(z))\ny: $(imag(z))"
		end,
	)

	marks = lift(view.marks) do zs
		xs, ys = to_image_coords(view, zs)
		return Point2f.(xs, ys)
	end

	lines!(
		view.axis,
		marks,
		color = (:blue, 0.5),
		inspectable = false,
	)

	scatter!(
		view.axis,
		marks,
		color = (:blue, 1.0),
		inspector_label = (self, i, p) -> let
			z = to_complex(view, p)
			"x: $(real(z))\ny: $(imag(z))"
		end,
	)

	DataInspector(view.axis)
end

@with_kw mutable struct MandelView <: View
	axis::Axis
	center::ComplexF64 = 0.0im
	diameter::Float64 = 4.0
	pixels::Int = 1000

	init_center::ComplexF64 = center
	init_diameter::Float64 = diameter

	array::Observable{Matrix{Float64}} =
		Observable(zeros(Float64, pixels, pixels))
	points::Observable{Vector{ComplexF64}} =
		Observable([init_center])
	marks::Observable{Vector{ComplexF64}} =
		Observable(ComplexF64[])

	inspectable::Observable{Any} = Observable{Any}(false)

	function MandelView(args...)
		view = new(args...)
		prepare!(view)
		return view
	end
end

@with_kw mutable struct JuliaView <: View
	axis::Axis
	center::ComplexF64 = 0.0im
	diameter::Float64 = 4.0
	parameter::ComplexF64 = 0.0im
	pixels::Int = 1000

	init_center::ComplexF64 = center
	init_diameter::Float64 = diameter

	array::Observable{Matrix{Float64}} =
		Observable(zeros(Float64, pixels, pixels))
	points::Observable{Vector{ComplexF64}} =
		Observable([init_center])
	marks::Observable{Vector{ComplexF64}} =
		Observable(ComplexF64[])

	inspectable::Observable{Any} = Observable{Any}(false)

	function JuliaView(args...)
		view = new(args...)
		prepare!(view)
		return view
	end
end

function orbit(f::Function, z::Number, c::Number, iterates::Integer)
	orbit = [z]
	for i in 1:iterates
		push!(orbit, f(orbit[end], c))
	end
	return orbit
end

function pick_orbit!(
	julia::JuliaView,
	d_system::DynamicalSystem,
	options::ViewerOptions,
	z::Number,
)
	julia.points[] = orbit(
		d_system.f,
		z,
		julia.parameter,
		options.orbit_len - 1,
	)
	return julia
end

function corner_and_step(view::View)
	step = view.diameter / view.pixels

	corner_real = real(view.center) - view.diameter / 2 + 0.5 * step
	corner_imag = imag(view.center) - view.diameter / 2 + 0.5 * step
	corner = complex(corner_real, corner_imag)

	return corner, step
end

function mandel_slice!(array, j, f, df_dz, crit, corner, step, pxs, options)
	@inbounds for i in 1:pxs
		c = corner + step * complex(i, j)
		array[i, j] = options.coloring_algorithm(
			f, df_dz, crit(c), c, options.esc_radius, options.max_iter
		)
	end
	return nothing
end

function update_grid!(
	view::MandelView,
	f::Function,
	df_dz::Function,
	crit::Function,
	corner::ComplexF64,
	step::Float64,
	options::ViewerOptions,
)
	futures = Vector{Task}(undef, view.pixels)

	@inbounds for j in 1:view.pixels
		futures[j] = Threads.@spawn mandel_slice!(
			view.array[], j, f, df_dz, crit, corner, step, view.pixels, options
		)
	end

	wait.(futures)
	return nothing
end

function julia_slice!(array, j, f, df_dz, c, corner, step, pxs, options)
	@inbounds for i in 1:pxs
		z = corner + step * complex(i, j)
		array[i, j] = options.coloring_algorithm(
			f, df_dz, z, c,	options.esc_radius,	options.max_iter,
		)
	end
	return nothing
end

function update_grid!(
	view::JuliaView,
	f::Function,
	df_dz::Function,
	crit::Function,
	corner::ComplexF64,
	step::Float64,
	options::ViewerOptions,
)
	futures = Vector{Task}(undef, view.pixels)

	@inbounds for j in 1:view.pixels
		futures[j] = Threads.@spawn julia_slice!(
			view.array[], j, f, df_dz, view.parameter, corner, step, view.pixels, options
		)
	end

	wait.(futures)
	return nothing
end

function update!(view::View, d_system::DynamicalSystem, options::ViewerOptions)
	corner, step = corner_and_step(view)
	update_grid!(
		view,
		d_system.f,
		d_system.df_dz,
		d_system.crit,
		corner,
		step,
		options,
	)
	notify(view.array)
	notify(view.points)
	notify(view.marks)

	return view
end

function centering!(
	view::View,
	d_system::DynamicalSystem,
	options::ViewerOptions,
	point
)
	view.center = to_complex(view, point)
	update!(view, d_system, options)
	return view
end

function zoom_in!(
	view::View,
	d_system::DynamicalSystem,
	options::ViewerOptions,
	point,
)
	z = to_complex(view, point)
	view.center = 0.5 * view.center + 0.5 * z
	view.diameter /= 2.0
	update!(view, d_system, options)
	return view
end

function zoom_out!(
	view::View,
	d_system::DynamicalSystem,
	options::ViewerOptions,
	point,
)
	z = to_complex(view, point)
	view.center = 2.0 * view.center - z
	view.diameter *= 2.0
	update!(view, d_system, options)
	return view
end

function pick_parameter!(
	julia::JuliaView,
	mandel::MandelView,
	d_system::DynamicalSystem,
	options::ViewerOptions,
	point,
)
	julia.parameter = to_complex(mandel, point)
	julia.init_center = d_system.crit(julia.parameter)
	if isapprox(julia.diameter, julia.init_diameter)
		julia.center = julia.init_center
	end

	julia.marks[] = orbit(
		d_system.f,
		julia.init_center,
		julia.parameter,
		options.crit_len - 1,
	)

	mandel.points[] = [julia.parameter]

	update!(julia, d_system, options)
	pick_orbit!(julia, d_system, options, julia.points[][begin])
	return julia
end

import Dates

function save_view(filename::String, view::View)
	pxs = view.pixels
	fig = Figure(figure_padding=0, size=(pxs, pxs))
	ax = Axis(fig[1, 1], aspect=AxisAspect(1))

	hidedecorations!(ax)
	hidespines!(ax)

	heatmap!(
		ax,
		view.array[],
		colormap = :twilight,
		colorrange = (0.0, 1.0)
	)

	Makie.save(filename, fig)
end

function reset!(
	view::View,
	d_system::DynamicalSystem,
	options::ViewerOptions
)
	view.center = view.init_center
	view.diameter = view.init_diameter
	update!(view, d_system, options)
	return view
end

function add_view_buttons(
	position,
	view::View,
	d_system::DynamicalSystem,
	options::ViewerOptions,
	pressed::Bool,
)
	buttons = GridLayout(
		position,
		width=Relative(0.95),
		valign=0.97,
		tellheight = false,
		tellwidth = false,
		default_colgap = 8,
	)

	save_button = Button(buttons[1, 2], label="Save")

	on(save_button.clicks, priority=200) do event
		format = Dates.dateformat"yyyy-mm-ddTHH.MM.SS"
		time = string(Dates.format(Dates.now(), format))
		save_view(time * ".png", view)
		pressed = true
		return Consume(true)
	end

	reset_button = Button(buttons[1, 3], label="Reset")

	on(reset_button.clicks, priority=200) do event
		reset!(view, d_system, options)
		pressed = true
		return Consume(true)
	end

end

struct Viewer
	d_system::DynamicalSystem
	options::ViewerOptions

	figure::Figure
	mandel::MandelView
	julia::JuliaView

	function Viewer(
			f;
            crit=0.0im,
			mandel_center=0.0im,
			mandel_diam=4.0,
			julia_diam=4.0,
			coloring_algorithm=escape_time,
		)
        d_system = DynamicalSystem(f, crit)
		options = ViewerOptions(100.0, 200, 1, 1, coloring_algorithm)
		state = ViewerState(:pick)
		figure = Figure(figure_padding=10, size=(900, 510))

		mandel = MandelView(
			axis = Axis(figure[1, 1][1, 1], aspect=AxisAspect(1)),
			center = mandel_center,
			diameter = mandel_diam,
			pixels = 1000,
		)

		julia = JuliaView(
			axis = Axis(figure[1, 1][1, 2], aspect=AxisAspect(1)),
			center = d_system.crit(mandel_center),
			diameter = julia_diam,
			parameter = mandel_center,
			pixels = 1000,
		)

		mandel.inspectable[] = state.last_button == :pick
		julia.inspectable[] = state.last_button == :pick

		mandel.points[] = [julia.parameter]
		julia.marks[] = [d_system.crit(julia.parameter)]
		pick_orbit!(
			julia,
			d_system,
			options,
			julia.points[][begin],
		)

		update!(mandel, d_system, options)
		update!(julia, d_system, options)

		colsize!(figure.layout, 1, Relative(0.97))
		rowsize!(figure.layout, 1, Aspect(1, 0.5))
		colgap!(content(figure[1,1]), 10)

		axis_pressed = false
		plots = figure[1, 1]
		add_view_buttons(plots[1, 1], mandel, d_system, options, axis_pressed)
		add_view_buttons(plots[1, 2], julia, d_system, options, axis_pressed)

		buttons = Dict(
			:centering => Button(figure[2,1][1,1], label="Choose\nCenter"),
			:zoom => Button(figure[2,1][1,2], label="Click\nZoom"),
			:pick => Button(figure[2,1][1,3], label="Choose\nPoint"),
		)

		labels = Dict(
			:max_iter => Label(figure[2,1][1,4], "Maximum\nIterations:"),
			:orbit_len => Label(figure[2,1][1,6], "Orbit\nLength:"),
			:crit_len => Label(figure[2,1][1,8], "Critical Point\nOrbit Length:"),
			:esc_radius => Label(figure[2,1][1,10], "Escape\nRadius:"),
		)

		input_fields = Dict(
			:max_iter => Textbox(
				figure[2,1][1,5],
				width = 60,
				placeholder = string(options.max_iter),
				validator = Int,
			),
			:orbit_len => Textbox(
				figure[2,1][1,7],
				width = 60,
				placeholder = string(options.orbit_len),
				validator = Int,
			),
			:crit_len => Textbox(
				figure[2,1][1,9],
				width = 60,
				placeholder = string(options.orbit_len),
				validator = Int,
			),
			:esc_radius => Textbox(
				figure[2,1][1,11],
				width = 60,
				placeholder = string(options.esc_radius),
				validator = Float64,
			),
		)

		on(input_fields[:max_iter].stored_string) do s
		    options.max_iter = parse(Int, s)
			update!(julia, d_system, options)
			update!(mandel, d_system, options)
		end

		on(input_fields[:orbit_len].stored_string) do s
		    options.orbit_len = parse(Int, s)
			pick_orbit!(
				julia,
				d_system,
				options,
				julia.points[][begin],
			)
		end

		on(input_fields[:crit_len].stored_string) do s
		    options.crit_len = parse(Int, s)
			julia.marks[] = orbit(
				d_system.f,
				julia.init_center,
				julia.parameter,
				options.orbit_len - 1,
			)
		end

		on(input_fields[:esc_radius].stored_string) do s
		    options.esc_radius = parse(Float64, s)
			update!(julia, d_system, options)
			update!(mandel, d_system, options)
		end

		button_non_active_color = buttons[:centering].buttoncolor[]
		label_non_active_color = buttons[:centering].labelcolor[]
		button_active_color = buttons[:centering].buttoncolor_active[]
		label_active_color = buttons[:centering].labelcolor_active[]

		buttons[state.last_button].buttoncolor[] = button_active_color
		buttons[state.last_button].labelcolor[] = label_active_color

		for (symbol, button1) in buttons
			on(button1.clicks) do n
				state.last_button = symbol

				mandel.inspectable[] = state.last_button == :pick
				julia.inspectable[] = state.last_button == :pick

				for (_, button2) in buttons
					button2.buttoncolor[] = button_non_active_color
					button2.labelcolor[] = label_non_active_color
				end
				button1.buttoncolor[] = button_active_color
				button1.labelcolor[] = label_active_color
			end
		end

		for view in [mandel, julia]
			scene = view.axis.scene
			axis = view.axis

			on(events(scene).mousebutton) do event
				if axis_pressed
					axis_pressed = false
					return Consume(false)
				end

				point = mouseposition(scene)
		        if event.button == Mouse.left && is_mouseinside(axis)
		            if event.action == Mouse.press
						if state.last_button == :zoom
							zoom_in!(view, d_system, options, point)
						elseif state.last_button == :centering
							centering!(view, d_system, options, point)
						end
					end
				elseif event.button == Mouse.right && is_mouseinside(axis)
					if event.action == Mouse.press && state.last_button == :zoom
						zoom_out!(view, d_system, options, point)
					end
				end
				return Consume(false)
			end
		end

		on(events(mandel.axis.scene).mousebutton) do event
			if axis_pressed
				axis_pressed = false
				return Consume(false)
			end

			point = mouseposition(mandel.axis.scene)
			if event.button == Mouse.left && is_mouseinside(mandel.axis)
				if event.action == Mouse.press && state.last_button == :pick
					pick_parameter!(julia, mandel, d_system, options, point)
				end
			end

			return Consume(false)
		end

		on(events(julia.axis.scene).mousebutton) do event
			if axis_pressed
				axis_pressed = false
				return Consume(false)
			end

			point = mouseposition(julia.axis.scene)
			if event.button == Mouse.left && is_mouseinside(julia.axis)
				if event.action == Mouse.press && state.last_button == :pick
					pick_orbit!(julia, d_system, options, to_complex(julia, point))
				end
			end

			return Consume(false)
		end

		return new(d_system, options, figure, mandel, julia)
	end
end

Base.show(io::IO, viewer::Viewer) = display(GLMakie.Screen(), viewer.figure)

# --------------------------------------------------------------------------------------- #
# 3D Viewer for Rational Maps
# --------------------------------------------------------------------------------------- #

mutable struct Viewer3DOptions
	ε ::Float64
	max_iter::Int
end

abstract type View3D end

mutable struct MandelView3D <: View3D
	axis::LScene
	grid::Matrix{Point}

	texture::Observable{Matrix{Float64}}
	camera::Camera3D

	crit::Function
	focus::Observable{Point}
	point::Observable{Point}
end

function mandel_slice3D!(grid, texture, w, j, δ, f, crit, options)
	@inbounds for i in 0:(w-1)
		φ = i * δ
		θ = j * δ

		pt = Point(sin(φ) * exp(im * θ), 1 - cos(φ))
		grid[i + 1, j + 1] = pt

		texture[i + 1, j + 1] = multiplier(
			f,
			crit(pt),
			pt,
			options.ε,
			options.max_iter,
		)
	end
	return nothing
end

function initialize!(
	view::MandelView3D,
	f::Function,
	crit::Function,
	options::Viewer3DOptions,
)
	w, h = size(view.grid)
	δ = pi / (w - 1)

	futures = Vector{Task}(undef, h)

	@inbounds for j in 0:(h-1)
		futures[j+1] = Threads.@spawn mandel_slice3D!(
			view.grid, view.texture[], w, j, δ, f, crit, options
		)
	end

	wait.(futures)

	return nothing
end

function plot_setup!(view::View3D)
	mesh!(
		Sphere(Point3f(0), 1.0),
		color = view.texture,
		colormap=:twilight,
		colorrange = (0.0, 1.0),
		shading = FastShading,
	)

	focus = lift(view.camera.eyeposition) do eye
		x, y, z = eye
		return eye / sqrt(x^2 + y^2 + z^2)
	end

	complex_focus = let
		x, y, z = convert(Vector{Float64}, focus[])
		u0 = complex(x, y)
		v0 = 1 - z
		Observable(Point(u0, v0))
	end

	point = lift(view.point) do point
		x = real(point.u)
		y = imag(point.u)
		z = 1 - real(point.v)
		return Vec3f(x,y,z)
	end

	on(focus) do f
		event = events(view.axis.scene).mousebutton[]
		if event.button == Mouse.left
			if event.action == Mouse.release
				x, y, z = convert(Vector{Float64}, f)
				u0 = complex(x, y)
				v0 = 1 - z

				complex_focus[] = Point(u0, v0)
			end
		end
	end

	scatter!(view.axis, focus, color = :red)
	scatter!(view.axis, point, color = :blue)

	return complex_focus
end

function MandelView3D(w::Int, figure::Figure, f::Function, crit::Function, options::Viewer3DOptions)
	h = 2 * w - 1

	axis = LScene(figure[1, 1], show_axis = false)
	grid = Matrix{Point}(undef, w, h)
	texture = Observable(Array{Float64}(undef, w, h))
	camera = Camera3D(axis.scene, eyeposition = Vec3f(0, 0, -3), upvector = Vec3f(0, -1, 0), fixed_axis = false, center = false)

	mandel = MandelView3D(axis, grid, texture, camera, crit, Observable(Point(0, 1)), Observable(Point(0, 1)))
	initialize!(mandel, f, crit, options)
	mandel.focus = plot_setup!(mandel)
	mandel.point[] = mandel.focus[]

	return mandel
end

mutable struct JuliaView3D <: View3D
	axis::LScene
	grid::Matrix{Point}

	texture::Observable{Matrix{Float64}}
	camera::Camera3D

	parameter::Observable{Point}
	focus::Observable{Point}
	point::Observable{Point}
end

function julia_slice3D!(grid, texture, w, j, δ, f, c, options)
	@inbounds for i in 0:(w-1)
		φ = i * δ
		θ = j * δ

		pt = Point(sin(φ) * exp(im * θ), 1 - cos(φ))
		grid[i + 1, j + 1] = pt

		texture[i + 1, j + 1] = multiplier(
			f, pt, c, options.ε, options.max_iter,
		)
	end

	return nothing
end

function initialize!(
	view::JuliaView3D,
	f::Function,
	parameter::Observable{Point},
	options::Viewer3DOptions,
)
	w, h = size(view.grid)
	δ = pi / (w - 1)

	on(parameter) do c
		futures = Vector{Task}(undef, h)

		@inbounds for j in 0:(h-1)
			futures[j+1] = Threads.@spawn julia_slice3D!(
				view.grid, view.texture[], w, j, δ, f, c, options
			)
		end

		wait.(futures)
		notify(view.texture)
	end
	return nothing
end

function JuliaView3D(w::Int, figure::Figure, f::Function, c::Observable{Point}, options::Viewer3DOptions)
	h = 2 * w - 1

	axis = LScene(figure[1, 2], show_axis = false)
	grid = Matrix{Point}(undef, w, h)
	texture = Observable(Array{Float64}(undef, w, h))
	camera = Camera3D(axis.scene, eyeposition = Vec3f(0, 0, -3), upvector = Vec3f(0, -1, 0), fixed_axis = false, center = false)

	julia = JuliaView3D(axis, grid, texture, camera, c, Observable(Point(0, 1)), Observable(Point(0, 1)))
	initialize!(julia, f, c, options)
	julia.focus = plot_setup!(julia)

	return julia
end

# function update!(
# 	texture::Observable{Matrix{Float64}},
# 	grid::Matrix{Point},
# 	f::Function,
# 	parameter::Number,
# 	options::Viewer3DOptions,
# 	u0::Number,
# 	v0::Number,
# 	r::Number,
# )

# 	Threads.@threads for i in eachindex(grid)
# 		grid[i] = Point(
# 			r * v0 * grid[i].u + (1 - r) * u0 * grid[i].v,
# 			v0 * grid[i].v
# 		)

# 		divide!(grid[i], sqrt(norm2(grid[i])))

# 		texture[][i] = multiplier(
# 			f,
# 			grid[i],
# 			parameter,
# 			options.ε,
# 			options.max_iter,
# 		)
# 	end

# 	notify(texture)
# 	return
# end

struct Viewer3D
	rational_map::RationalMap
	figure::Figure
	options::Viewer3DOptions

	mandel::MandelView3D
	julia::JuliaView3D
end

function Viewer3D(f::Function; crit=0.0im, c=0.0im)
    rational_map = RationalMap(f, crit)
    options = Viewer3DOptions(1e-4, 200)
    figure = Figure(size = (1000, 560))

	# Initialize Mandel and Julia Views
	w = 501
    mandel = MandelView3D(w, figure, rational_map.f, rational_map.crit, options)
	julia = JuliaView3D(w, figure, rational_map.f, mandel.point, options)
	notify(mandel.point)

	rowsize!(figure.layout, 1, Aspect(1, 1))
	colgap!(figure.layout, 5)

	menu_mandel = GridLayout(
		figure[2,1],
		width=Relative(0.95),
		valign=0.97,
		tellheight = false,
		tellwidth = false,
		default_colgap = 8,
	)

	point_button = Button(menu_mandel[1,1], label="Pick Parameter")

	on(point_button.clicks) do event
		mandel.point[] = mandel.focus[]
	end

    # zoom_in = Button(figure[2,1][1,1], label="Zoom In")
    # zoom_out = Button(figure[2,1][1,2], label="Zoom Out")

    # on(zoom_in.clicks) do event
    #     x, y, z = convert(Vector{Float64}, focus[])
    #     u0 = complex(x, y)
    #     v0 = 1 - z
    #     update!(texture, grid, rational_map.f, parameter, options, u0, v0, 0.5)
    # end

    # on(zoom_out.clicks) do event
    #     x, y, z = convert(Vector{Float64}, focus[])
    #     u0 = complex(x, y)
    #     v0 = 1 - z
    #     update!(texture, grid, rational_map.f, parameter, options, u0, v0, 2.0)
    # end

    return Viewer3D(rational_map, figure, options, mandel, julia)
end

Base.show(io::IO, viewer::Viewer3D) = display(GLMakie.Screen(), viewer.figure)

end