mutable struct State
    ε ::Float64
	max_iter::Int
end

mutable struct View3D
    scene::LScene
    camera::Camera3D

    points::Matrix{Point}
    texture::Observable{Matrix{Float64}}

    sphere_focus::Observable{Vec3f}
    point_focus::Observable{Point}

    sphere_mark::Observable{Vec3f}
    point_mark::Observable{Point}
end

function points_grid(meridian, longitude)
    points = Matrix{Point}(undef, meridian, longitude)
    δ = pi / (meridian - 1)

    for c in CartesianIndices(points)
        φ = (c[1] - 1) * δ
        θ = (c[2] - 1) * δ

        points[c] = Point(sin(φ) * exp(im * θ), 1 - cos(φ))
    end

    return points
end

function update_slice!(texture, height, column, points, f, crit::Function, state)
    @inbounds for row in 1:height
        pt = points[row, column]
        texture[row, column] = multiplier(
            f, crit(pt), pt, state.ε, state.max_iter
        )
    end
end

function update_slice!(texture, height, column, points, f, c::Observable{Point}, state)
    @inbounds for row in 1:height
        pt = points[row, column]
        texture[row, column] = multiplier(
            f, pt, c[], state.ε, state.max_iter,
        )
    end
end

function update_view!(texture, points, f, c, state)
    height, width = size(points)
    futures = Vector{Task}(undef, width)

    @inbounds for column in 1:width
        futures[column] = Threads.@spawn update_slice!(
            texture[], height, column, points, f, c, state
        )
    end

    wait.(futures)
    notify(texture)
end

function plot_setup!(scene, points, texture)
    mesh!(
        scene,
        Sphere(Point3f(0), 1.0),
        color = texture,
        colormap=:twilight,
        colorrange = (0.0, 1.0),
        shading = FastShading,
    )

    camera = Camera3D(scene.scene, eyeposition = Vec3f(0, 0, -3), upvector = Vec3f(0, -1, 0), fixed_axis = false, center = false)

    sphere_focus = lift(camera.eyeposition) do eye
        return normalize(eye)
    end

    point_focus = lift(sphere_focus) do f
        x, y, z = convert(Vector{Float64}, f)
        return Point(complex(x, y), 1 - z)
    end

    sphere_mark = Observable(Vec3f(0, 0, -1))

    point_mark = lift(sphere_mark) do m
        x, y, z = convert(Vector{Float64}, m)
        return Point(complex(x, y), 1 - z)
    end

    scatter!(scene, sphere_focus, color = :red)
    scatter!(scene, sphere_mark, color = :blue)

    return View3D(scene, camera, points, texture, sphere_focus, point_focus, sphere_mark, point_mark)
end

function view3D(meridian, figure, f, c, state)
    longitude = 2 * meridian - 1

	scene = LScene(figure, show_axis = false)
	points = points_grid(meridian, longitude)

    texture = Observable(Matrix{Float64}(undef, meridian, longitude))
    update_view!(texture, points, f, c, state)
    view = plot_setup!(scene, points, texture)

	return view
end

struct Viewer3D
	rational_map::RationalMap
	figure::Figure
	state::State

	mandel::View3D
	julia::View3D
end

function update_points!(points, texture, mobius, f, c, state)
    for point in points
        point .= normalize(mobius * point)
    end
    update_view!(texture, points, f, c, state)
end

function put_menu!(figure, view, rational_map, c, state)
    menu = GridLayout(
		figure,
		width=Relative(0.95),
		valign=0.97,
		tellheight = false,
		tellwidth = false,
		default_colgap = 8,
	)

    zoominbutton = Button(menu[1,1], label="↑", width=30)

    on(zoominbutton.clicks) do event
        mobius = antipodal_hyperbolic(view.point_focus[], 0.5)
        update_points!(view.points, view.texture, mobius, rational_map.f, c, state)
	end

    zoomoutbutton = Button(menu[1,2], label="↓", width=30)

    on(zoomoutbutton.clicks) do event
        mobius = antipodal_hyperbolic(view.point_focus[], 2)
        update_points!(view.points, view.texture, mobius, rational_map.f, c, state)
	end

    return menu
end

function Viewer3D(f::Function; crit=0.0im, c=0.0im)
    rational_map = RationalMap(f, crit)
    state = State(1e-4, 200)
    figure = Figure(size = (1000, 560))

	# Initialize Mandel and Julia Views
	meridian = 501
    mandel = view3D(meridian, figure[1, 1], rational_map.f, rational_map.crit, state)
	julia = view3D(meridian, figure[1, 2], rational_map.f, mandel.point_mark, state)

	rowsize!(figure.layout, 1, Aspect(1, 1))
	colgap!(figure.layout, 5)

	mandel_menu = put_menu!(figure[2,1], mandel, rational_map, rational_map.crit, state)
    julia_menu = put_menu!(figure[2,2], julia, rational_map, mandel.point_mark, state)

    point_button = Button(mandel_menu[1,3], label="·", width=30)
	on(point_button.clicks) do event
		mandel.sphere_mark[] = mandel.sphere_focus[]
        update_view!(julia.texture, julia.points, rational_map.f, mandel.point_mark, state)
	end

    return Viewer3D(rational_map, figure, state, mandel, julia)
end

Base.show(io::IO, viewer::Viewer3D) = display(GLMakie.Screen(), viewer.figure)