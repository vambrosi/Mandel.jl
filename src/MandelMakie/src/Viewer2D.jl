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
	a = (x + 0.5 - view.pixels / 2) / ppu
	b = (y + 0.5 - view.pixels / 2) / ppu

	return view.center + complex(a, b)
end

function to_image_coords(view::View, z::Number)
	ppu = view.pixels / view.diameter
	w = z - view.center

	x = real(w) * ppu - 0.5 + view.pixels / 2
	y = imag(w) * ppu - 0.5 + view.pixels / 2

	return [x, y]
end

function to_image_coords(view::View, zs::Vector{<:Number})
	ppu = view.pixels / view.diameter
	xs = Float64[]
	ys = Float64[]

	for z in zs
		w = z - view.center
		push!(xs, real(w) * ppu - 0.5 + view.pixels / 2)
		push!(ys, imag(w) * ppu - 0.5 + view.pixels / 2)
	end

	return xs, ys
end

function prepare!(view::View)
	hidedecorations!(view.axis)
	hidespines!(view.axis)
	deregister_interaction!(view.axis, :rectanglezoom)
	deregister_interaction!(view.axis, :dragpan)

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

	is_zooming::Bool = false
	zooming::Timer = Timer(identity, 0.1)

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

	is_zooming::Bool = false
	zooming::Timer = Timer(identity, 0.1)

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

function update_view!(view::View, d_system::DynamicalSystem, options::ViewerOptions)
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

function zoom!(view::View, d_system::DynamicalSystem, options::ViewerOptions)
	view.is_zooming = true

	x_min, x_max, _, _ = view.axis.limits[]
	original_size = x_max - x_min
	x1_min, x1_max = view.axis.xaxis.attributes.limits[]
	current_size = x1_max - x1_min

	scale = current_size / original_size
	point = mouseposition(view.axis.scene)
	z = to_complex(view, point)

	view.diameter *= scale
	view.center = scale * view.center + (1 - scale) * z
	update_view!(view, d_system, options)
	reset_limits!(view.axis)
	view.is_zooming = false

	return view
end

function pick_parameter!(
	julia::JuliaView,
	mandel::MandelView,
	d_system::DynamicalSystem,
	options::ViewerOptions,
	point,
)
	julia.parameter = point

	julia.marks[] = orbit(
		d_system.f,
		d_system.crit(julia.parameter),
		julia.parameter,
		options.crit_len - 1,
	)

	mandel.points[] = [julia.parameter]

	update_view!(julia, d_system, options)
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
	update_view!(view, d_system, options)
	reset_limits!(view.axis)
	return view
end

function add_view_buttons(
	position,
	view::View,
	d_system::DynamicalSystem,
	options::ViewerOptions,
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

	on(save_button.clicks, priority=300) do event
		format = Dates.dateformat"yyyy-mm-ddTHH.MM.SS"
		time = string(Dates.format(Dates.now(), format))
		!isdir("imgs") && mkdir("imgs")
		save_view("imgs/" * time * ".png", view)
		return Consume(false)
	end

	reset_button = Button(buttons[1, 3], label="Reset")

	on(reset_button.clicks, priority=300) do event
		if !view.is_zooming
			reset!(view, d_system, options)
		end

		return Consume(false)
	end
end

function add_view_events(view::View, other::View, d_system::DynamicalSystem, options::ViewerOptions)
	dragging = Ref(false)
	dragstart = Ref(Point2f(0.0))
	dragend = Ref(Point2f(0.0))
	to_world_at_start = z -> to_world(view.axis.scene, z)

	on(events(view.axis.scene).mousebutton) do event
		mp = mouseposition_px(view.axis.scene)

		if event.button == Mouse.right
			if event.action == Mouse.press && is_mouseinside(view.axis)
				dragging[] = true
				to_world_at_start = z -> to_world(view.axis.scene, z)
				dragstart[] = to_world_at_start(mp)

			elseif event.action == Mouse.release && dragging[]
				dragend[] = to_world_at_start(mp)
				view.center += to_complex(view, dragstart[]) - to_complex(view, dragend[])
				translate!(view.axis.scene, Point2(0))
				update_view!(view, d_system, options)
				reset_limits!(view.axis)
				dragging[] = false
			end
		end

		if event.button == Mouse.left
			if event.action == Mouse.press && is_mouseinside(view.axis)
				point = to_complex(view, to_world_at_start(mp))
				if view isa MandelView
					pick_parameter!(other, view, d_system, options, point)
				elseif view isa JuliaView
					pick_orbit!(view, d_system, options, point)
				end
			end
		end
	end

	on(events(view.axis.scene).mouseposition) do event
		if dragging[]
			mp = mouseposition_px(view.axis.scene)
			translate!(view.axis.scene, to_world_at_start(mp) - dragstart[])
		end
	end

	on(events(view.axis.scene).scroll, priority=100) do event
		if is_mouseinside(view.axis) && !view.is_zooming
			close(view.zooming)
			view.zooming = Timer(_ -> zoom!(view, d_system, options), 0.1)
		end
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
			c=0.0im,
			mandel_center=0.0im,
			mandel_diam=4.0,
			julia_center=0.0im,
			julia_diam=4.0,
			coloring_algorithm=:escape_time,
		)

		algs = Dict(
			:escape_time => escape_time,
			:stop_time => stop_time,
			:escape_preperiod => escape_preperiod
		)

        d_system = DynamicalSystem(f, crit)
		options = ViewerOptions(100.0, 200, 1, 1, algs[coloring_algorithm])
		state = ViewerState(:zoom)
		figure = Figure(figure_padding=10, size=(900, 510))

		mandel = MandelView(
			axis = Axis(figure[1, 1][1, 1], aspect=AxisAspect(1), limits = (0.5, 1000.5, 0.5, 1000.5)),
			center = mandel_center,
			diameter = mandel_diam,
			pixels = 1000,
		)

		julia = JuliaView(
			axis = Axis(figure[1, 1][1, 2], aspect=AxisAspect(1), limits = (0.5, 1000.5, 0.5, 1000.5)),
			center = julia_center,
			diameter = julia_diam,
			parameter = c,
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

		update_view!(mandel, d_system, options)
		update_view!(julia, d_system, options)

		colsize!(figure.layout, 1, Relative(0.97))
		rowsize!(figure.layout, 1, Aspect(1, 0.5))
		colgap!(content(figure[1,1]), 10)

		plots = figure[1, 1]
		add_view_buttons(plots[1, 1], mandel, d_system, options)
		add_view_buttons(plots[1, 2], julia, d_system, options)

		labels = Dict(
			:max_iter => Label(figure[2,1][1,1], "Maximum\nIterations:"),
			:orbit_len => Label(figure[2,1][1,3], "Orbit\nLength:"),
			:crit_len => Label(figure[2,1][1,5], "Critical Point\nOrbit Length:"),
			:esc_radius => Label(figure[2,1][1,7], "Escape\nRadius:"),
		)

		input_fields = Dict(
			:max_iter => Textbox(
				figure[2,1][1,2],
				width = 60,
				placeholder = string(options.max_iter),
				validator = Int,
			),
			:orbit_len => Textbox(
				figure[2,1][1,4],
				width = 60,
				placeholder = string(options.orbit_len),
				validator = Int,
			),
			:crit_len => Textbox(
				figure[2,1][1,6],
				width = 60,
				placeholder = string(options.orbit_len),
				validator = Int,
			),
			:esc_radius => Textbox(
				figure[2,1][1,8],
				width = 60,
				placeholder = string(options.esc_radius),
				validator = Float64,
			),
		)

		on(input_fields[:max_iter].stored_string) do s
		    options.max_iter = parse(Int, s)
			update_view!(julia, d_system, options)
			update_view!(mandel, d_system, options)
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
				d_system.crit(julia.parameter),
				julia.parameter,
				options.crit_len - 1,
			)
		end

		on(input_fields[:esc_radius].stored_string) do s
		    options.esc_radius = parse(Float64, s)
			update_view!(julia, d_system, options)
			update_view!(mandel, d_system, options)
		end

		add_view_events(mandel, julia, d_system, options)
		add_view_events(julia, mandel, d_system, options)

		return new(d_system, options, figure, mandel, julia)
	end
end

Base.show(io::IO, viewer::Viewer) = display(GLMakie.Screen(), viewer.figure)