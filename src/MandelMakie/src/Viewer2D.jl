# --------------------------------------------------------------------------------------- #
# Views of Mandelbrot and Julia Sets
# --------------------------------------------------------------------------------------- #

abstract type View end

mutable struct Options
	escape_radius::Float64
	max_iterations::Int
	orbit_length::Int
	critical_length::Int
	compact_view::Bool
	coloring_algorithm::Function
end

mutable struct MandelView <: View
	center::ComplexF64
	diameter::Float64
	pixels::Int

	init_center::ComplexF64
	init_diameter::Float64

	color_levels::Observable{Matrix{Float64}}
	points::Observable{Vector{ComplexF64}}
	marks::Observable{Vector{ComplexF64}}

	function MandelView(; center = 0.0im, diameter = 4.0, pixels = 1000)
		diameter > 0.0 || throw("diameter must be a positive real")
		pixels > 0 || throw("pixels must be a positive integer")

		color_levels = zeros(Float64, pixels, pixels)
		points = ComplexF64[center]
		marks = ComplexF64[]

		return new(
			center, diameter, pixels, center,
			diameter, color_levels, points, marks
		)
	end
end

mutable struct JuliaView <: View
	center::ComplexF64
	diameter::Float64
	parameter::ComplexF64
	pixels::Int

	init_center::ComplexF64
	init_diameter::Float64

	color_levels::Observable{Matrix{Float64}}
	points::Observable{Vector{ComplexF64}}
	marks::Observable{Vector{ComplexF64}}

	function JuliaView(; center = 0.0im, diameter = 4.0, parameter = 0.0im, pixels = 1000)
		diameter > 0.0 || throw("diameter must be a positive real")
		pixels > 0 || throw("pixels must be a positive integer")

		color_levels = zeros(Float64, pixels, pixels)
		points = ComplexF64[center]
		marks = ComplexF64[]

		return new(
			center, diameter, parameter, pixels, center,
			diameter, color_levels, points, marks
		)
	end
end

# --------------------------------------------------------------------------------------- #
# Change of Coordinates between Complex Plane and Pixel Space
# --------------------------------------------------------------------------------------- #

function to_complex_plane(view::View, pixel_vector)
	x,  y = pixel_vector[1], pixel_vector[2]

	upp = view.diameter / view.pixels
	a = (x + 0.5 - view.pixels / 2) * upp
	b = (y + 0.5 - view.pixels / 2) * upp

	return view.center + complex(a, b)
end

function to_pixel_space(view::View, z::Number)
	ppu = view.pixels / view.diameter
	w = z - view.center

	x = real(w) * ppu - 0.5 + view.pixels / 2
	y = imag(w) * ppu - 0.5 + view.pixels / 2

	return [x, y]
end

function to_pixel_space(view::View, zs::Vector{<:Number})
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

# --------------------------------------------------------------------------------------- #
# Updating Plots
# --------------------------------------------------------------------------------------- #

function pick_orbit!(
	julia::JuliaView,
	d_system::DynamicalSystem,
	options::Options,
	z::Number,
)
	julia.points[] = orbit(
		d_system.f,
		z,
		julia.parameter,
		options.orbit_length - 1,
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
			f, df_dz, crit(c), c, options.escape_radius, options.max_iterations
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
	options::Options,
)
	futures = Vector{Task}(undef, view.pixels)

	@inbounds for j in 1:view.pixels
		futures[j] = Threads.@spawn mandel_slice!(
			view.color_levels[], j, f, df_dz, crit, corner, step, view.pixels, options
		)
	end

	wait.(futures)
	return nothing
end

function julia_slice!(array, j, f, df_dz, c, corner, step, pxs, options)
	@inbounds for i in 1:pxs
		z = corner + step * complex(i, j)
		array[i, j] = options.coloring_algorithm(
			f, df_dz, z, c,	options.escape_radius, options.max_iterations,
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
	options::Options,
)
	futures = Vector{Task}(undef, view.pixels)

	@inbounds for j in 1:view.pixels
		futures[j] = Threads.@spawn julia_slice!(
			view.color_levels[], j, f, df_dz, view.parameter, corner, step, view.pixels, options
		)
	end

	wait.(futures)
	return nothing
end

function update_view!(view::View, d_system::DynamicalSystem, options::Options)
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
	notify(view.color_levels)
	notify(view.points)
	notify(view.marks)

	return view
end

function pick_parameter!(
	julia::JuliaView,
	mandel::MandelView,
	d_system::DynamicalSystem,
	options::Options,
	point,
)
	julia.parameter = point

	julia.marks[] = orbit(
		d_system.f,
		d_system.crit(julia.parameter),
		julia.parameter,
		options.critical_length - 1,
	)

	mandel.points[] = [julia.parameter]

	update_view!(julia, d_system, options)
	pick_orbit!(julia, d_system, options, julia.points[][begin])
	return julia
end

# --------------------------------------------------------------------------------------- #
# User Interface
# --------------------------------------------------------------------------------------- #

struct Frame
	axis::Axis
	view::Ref{Any}
	events::Dict{Symbol, Any}
end

function translate_axis!(axis, z)
	translate!(axis.scene, 0, 0, z)
	translate!(axis.elements[:background], 0, 0, z-1)
	translate!(axis.elements[:xoppositeline], 0, 0, z+1)
	translate!(axis.elements[:yoppositeline], 0, 0, z+1)
	translate!(axis.xaxis.elements[:axisline], 0, 0, z+1)
	translate!(axis.yaxis.elements[:axisline], 0, 0, z+1)
end

function create_frames!(figure, options, mandel, julia)
	mandel_limits = (0.5, 0.5 + mandel.pixels, 0.5, 0.5 + mandel.pixels)
	julia_limits = (0.5, 0.5 + julia.pixels, 0.5, 0.5 + julia.pixels)

	left_axis = Axis(figure[1,1][1,1], aspect=1, limits=mandel_limits)
	right_axis = Axis(figure[1,1][1,1], aspect=1, limits=julia_limits)
	translate_axis!(left_axis, 10)
	translate_axis!(right_axis, 0)

	if options.compact_view
		left_axis.width = Relative(0.3)
		left_axis.height = Relative(0.3)
		left_axis.valign = 0.03
		left_axis.halign = 0.03
	end

	for axis in [left_axis, right_axis]
		hidedecorations!(axis)
		deregister_interaction!(axis, :rectanglezoom)
		deregister_interaction!(axis, :dragpan)
	end

	return Frame(left_axis, mandel, Dict()), Frame(right_axis, julia, Dict())
end

function create_plot!(frame::Frame)
	empty!(frame.axis)
	view = frame.view[]

	heatmap!(
		frame.axis,
		view.color_levels,
		colormap = :twilight,
		colorrange = (0.0, 1.0),
		inspectable = false,
		inspector_label = (self, i, p) -> let
			z = to_complex_plane(view, p)
			"x: $(real(z))\ny: $(imag(z))"
		end,
	)

	point_vectors = lift(view.points) do zs
		xs, ys = to_pixel_space(view, zs)
		return Point2f.(xs, ys)
	end

	lines!(
		frame.axis,
		point_vectors,
		color = (:red, 0.5),
		inspectable = false,
	)

	scatter!(
		frame.axis,
		point_vectors,
		color = (:red, 1.0),
		inspector_label = (self, i, p) -> let
			z = to_complex_plane(view, p)
			"x: $(real(z))\ny: $(imag(z))"
		end,
	)

	mark_vectors = lift(view.marks) do zs
		xs, ys = to_pixel_space(view, zs)
		return Point2f.(xs, ys)
	end

	lines!(
		frame.axis,
		mark_vectors,
		color = (:blue, 0.5),
		inspectable = false,
	)

	scatter!(
		frame.axis,
		mark_vectors,
		color = (:blue, 1.0),
		inspector_label = (self, i, p) -> let
			z = to_complex_plane(view, p)
			"x: $(real(z))\ny: $(imag(z))"
		end,
	)

	return frame
end

function zoom!(frame::Frame, d_system::DynamicalSystem, options::Options)
	frame.events[:is_zooming] = true
	view = frame.view[]

	x_min, x_max, _, _ = frame.axis.limits[]
	original_size = x_max - x_min
	x1_min, x1_max = frame.axis.xaxis.attributes.limits[]
	current_size = x1_max - x1_min

	scale = current_size / original_size
	vector = mouseposition(frame.axis.scene)
	z = to_complex_plane(view, vector)

	view.diameter *= scale
	view.center = scale * view.center + (1 - scale) * z
	update_view!(view, d_system, options)
	reset_limits!(frame.axis)
	frame.events[:is_zooming] = false

	return frame
end

# import Dates

# function save_view(filename::String, view::View)
# 	pxs = view.pixels
# 	fig = Figure(figure_padding=0, size=(pxs, pxs))
# 	ax = Axis(fig[1, 1], aspect=AxisAspect(1))

# 	hidedecorations!(ax)
# 	hidespines!(ax)

# 	heatmap!(
# 		ax,
# 		view.color_levels[],
# 		colormap = :twilight,
# 		colorrange = (0.0, 1.0)
# 	)

# 	Makie.save(filename, fig)
# end

# function reset!(
# 	view::View,
# 	d_system::DynamicalSystem,
# 	options::Options
# )
# 	view.center = view.init_center
# 	view.diameter = view.init_diameter
# 	update_view!(view, d_system, options)
# 	reset_limits!(view.axis)
# 	return view
# end

# function add_view_buttons(
# 	position,
# 	view::View,
# 	d_system::DynamicalSystem,
# 	options::Options,
# )
# 	buttons = GridLayout(
# 		position,
# 		width=Relative(0.95),
# 		valign=0.97,
# 		tellheight = false,
# 		tellwidth = false,
# 		default_colgap = 8,
# 	)

# 	save_button = Button(buttons[1, 2], label="Save")

# 	on(save_button.clicks, priority=300) do event
# 		format = Dates.dateformat"yyyy-mm-ddTHH.MM.SS"
# 		time = string(Dates.format(Dates.now(), format))
# 		!isdir("imgs") && mkdir("imgs")
# 		save_view("imgs/" * time * ".png", view)
# 		return Consume(false)
# 	end

# 	reset_button = Button(buttons[1, 3], label="Reset")

# 	on(reset_button.clicks, priority=300) do event
# 		if !view.is_zooming
# 			reset!(view, d_system, options)
# 		end

# 		return Consume(false)
# 	end
# end

function add_frame_events!(
	frame::Frame,
	topframe::Frame,
	d_system::DynamicalSystem,
	options::Options,
	julia::JuliaView,
)
	axis = frame.axis
	scene = axis.scene

	dragging = false
	dragstart = Point2f(0.0)
	dragend = Point2f(0.0)

	to_world_at_start = z -> to_world(scene, z)

	is_zooming = false
	zooming = Timer(identity, 0.1)

	is_topframe = frame == topframe
	z_level = is_topframe ? 10 : 0

	on(events(scene).mousebutton) do event
		mp = mouseposition_px(scene)
		view = frame.view[]

		if event.button == Mouse.right
			if event.action == Mouse.press && is_mouseinside(axis) &&
					(is_topframe || !is_mouseinside(topframe.axis))
				dragging = true
				to_world_at_start = z -> to_world(scene, z)
				dragstart = to_world_at_start(mp)

			elseif event.action == Mouse.release && dragging
				dragend = to_world_at_start(mp)
				view.center += to_complex_plane(view, dragstart) - to_complex_plane(view, dragend)
				translate!(scene, 0, 0, z_level)
				update_view!(view, d_system, options)
				reset_limits!(axis)
				dragging = false
			end
		end

		if event.button == Mouse.left
			if event.action == Mouse.press && is_mouseinside(axis) &&
					(is_topframe || !is_mouseinside(topframe.axis))
				point = to_complex_plane(view, to_world_at_start(mp))
				if view isa MandelView
					pick_parameter!(julia, view, d_system, options, point)
				elseif view isa JuliaView
					pick_orbit!(view, d_system, options, point)
				end
			end
		end
	end

	on(events(scene).mouseposition) do event
		if dragging
			mp = mouseposition(scene)
			translate!(scene, mp - dragstart..., z_level)
		end
	end

	on(events(scene).scroll, priority=100) do event
		if is_mouseinside(axis) && !is_zooming &&
				(is_topframe || is_mouseinside(topframe.axis))
			close(zooming)
			zooming = Timer(_ -> zoom!(frame, d_system, options), 0.1)
		end
	end

	frame.events[:dragging] = dragging
	frame.events[:dragstart] = dragstart
	frame.events[:dragend] = dragstart
	frame.events[:is_zooming] = is_zooming
	frame.events[:zooming] = zooming
end

function add_buttons!(figure, left_frame, right_frame, mandel, julia, d_system, options)
	layout = GridLayout(figure[2,1])

	labels = Dict(
		:max_iter => Label(layout[1,3], "Maximum\nIterations:"),
		:orbit_len => Label(layout[1,5], "Orbit\nLength:"),
		:critical_length => Label(layout[1,7], "Critical Point\nOrbit Length:"),
		:escape_radius => Label(layout[1,9], "Escape\nRadius:"),
	)

	inputs = Dict(
		:switch_layout => Button(layout[1,1], label="↰", halign=:left),
		:switch_positions => Button(layout[1,2], label="↔", halign=:left),
		:max_iter => Textbox(
			layout[1,4],
			width = 60,
			placeholder = string(options.max_iterations),
			validator = Int,
		),
		:orbit_len => Textbox(
			layout[1,6],
			width = 60,
			placeholder = string(options.orbit_length),
			validator = Int,
		),
		:critical_length => Textbox(
			layout[1,8],
			width = 60,
			placeholder = string(options.orbit_length),
			validator = Int,
		),
		:escape_radius => Textbox(
			layout[1,10],
			width = 60,
			placeholder = string(options.escape_radius),
			validator = Float64,
		),
	)

	on(inputs[:switch_layout].clicks, priority=300) do event
		if options.compact_view
			left_frame.axis.width = nothing
			left_frame.axis.height = nothing
			left_frame.axis.valign = :center
			left_frame.axis.halign = :center

			figure[1, 1][1, 1] = left_frame.axis
			figure[1, 1][1, 2] = right_frame.axis

			inputs[:switch_layout].label = "↳"
		else
			left_frame.axis.width = Relative(0.3)
			left_frame.axis.height = Relative(0.3)
			left_frame.axis.valign = 0.01
			left_frame.axis.halign = 0.01

			figure[1, 1][1, 1] = left_frame.axis
			figure[1, 1][1, 1] = right_frame.axis

			trim!(contents(figure[1,1])...)
			inputs[:switch_layout].label = "↰"
		end

		options.compact_view = !options.compact_view
	end

	on(inputs[:switch_positions].clicks, priority=300) do event
		left_frame.view[], right_frame.view[] =
			right_frame.view[], left_frame.view[]
		create_plot!(left_frame)
		create_plot!(right_frame)
	end

	on(inputs[:max_iter].stored_string) do s
		options.max_iterations = parse(Int, s)
		update_view!(julia, d_system, options)
		update_view!(mandel, d_system, options)
	end

	on(inputs[:orbit_len].stored_string) do s
		options.orbit_length = parse(Int, s)
		pick_orbit!(
			julia,
			d_system,
			options,
			julia.points[][begin],
		)
	end

	on(inputs[:critical_length].stored_string) do s
		options.critical_length = parse(Int, s)
		julia.marks[] = orbit(
			d_system.f,
			d_system.crit(julia.parameter),
			julia.parameter,
			options.critical_length - 1,
		)
	end

	on(inputs[:escape_radius].stored_string) do s
		options.escape_radius = parse(Float64, s)
		update_view!(julia, d_system, options)
		update_view!(mandel, d_system, options)
	end

	return inputs
end

struct Viewer
	d_system::DynamicalSystem
	options::Options

	figure::Figure
	left_frame::Frame
	right_frame::Frame

	mandel::MandelView
	julia::JuliaView
	inputs::Dict{Symbol, Any}

	function Viewer(
			f;
            crit=0.0im,
			c=0.0im,
			mandel_center=0.0im,
			mandel_diameter=4.0,
			julia_center=0.0im,
			julia_diameter=4.0,
			compact_view=true,
			coloring_algorithm=:escape_time,
		)

		algs = Dict(
			:escape_time => escape_time,
			:stop_time => stop_time,
			:escape_preperiod => escape_preperiod
		)

        d_system = DynamicalSystem(f, crit)
		options = Options(100.0, 200, 1, 1, compact_view, algs[coloring_algorithm])
		figure = Figure(figure_padding=10, size=(800, 750))

		mandel = MandelView(
			center = mandel_center, diameter = mandel_diameter, pixels = 1000
		)

		julia = JuliaView(
			center = julia_center, diameter = julia_diameter, parameter = c, pixels = 1000
		)

		mandel.points[] = [julia.parameter]
		julia.marks[] = [d_system.crit(julia.parameter)]
		pick_orbit!(
			julia,
			d_system,
			options,
			julia.points[][begin],
		)

		left_frame, right_frame = create_frames!(figure, options, mandel, julia)
		create_plot!(left_frame)
		create_plot!(right_frame)

		add_frame_events!(left_frame, left_frame, d_system, options, julia)
		add_frame_events!(right_frame, left_frame, d_system, options, julia)

		update_view!(mandel, d_system, options)
		update_view!(julia, d_system, options)

		# colsize!(figure.layout, 1, Relative(0.97))
		# rowsize!(figure.layout, 1, Aspect(1, 0.5))
		colgap!(content(figure[1,1]), 10)

		# plots = figure[1, 1]
		# add_view_buttons(plots[1, 1], mandel, d_system, options)
		# add_view_buttons(plots[1, 2], julia, d_system, options)

		inputs = add_buttons!(figure, left_frame, right_frame, mandel, julia, d_system, options)

		return new(d_system, options, figure, left_frame, right_frame, mandel, julia, inputs)
	end
end

Base.show(io::IO, viewer::Viewer) = display(GLMakie.Screen(), viewer.figure)