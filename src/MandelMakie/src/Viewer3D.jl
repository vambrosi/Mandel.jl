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
		x = real(point[1])
		y = imag(point[1])
		z = 1 - real(point[2])
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