mutable struct State
    ε::Float64
    max_iter::Int
end

mutable struct View3D
    scene::LScene
    camera::Camera3D

    points::Observable{Matrix{Point}}
    texture::Observable{Matrix{Float64}}

    focus_vector::Observable{Vec3f}
    focus_point::Observable{Point}

    mark_vector::Observable{Vec3f}
    mark_point::Observable{Point}
end

function reset_points!(points)
    longitudes, _ = size(points[])
    δ = pi / (longitudes - 1)

    for c in CartesianIndices(points[])
        φ = (c[1] - 1) * δ
        θ = (c[2] - 1) * δ

        if φ == 0
            points[][c] = Point(1, 0)
        else
            points[][c] = normalize(Point(sin(φ) * exp(im * θ), 1 - cos(φ)))
        end
    end

    notify(points)
    return points
end

function update_texture!(texture, points, f, crit::Function, state)
    height, width = size(points)

    Threads.@threads for column in 1:width
        for row in 1:height
            pt = points[row, column]
            @inbounds texture[][row, column] = multiplier(f, crit(pt), pt, state.ε, state.max_iter)
        end
    end

    notify(texture)
    return texture
end

function update_texture!(texture, points, f, c::Observable{Point}, state)
    height, width = size(points)

    Threads.@threads for column in 1:width
        for row in 1:height
            pt = points[row, column]
            @inbounds texture[][row, column] = multiplier(f, pt, c[], state.ε, state.max_iter)
        end
    end

    notify(texture)
    return texture
end

function indices_to_vector(indices, longitudes)
    δ = pi / (longitudes - 1)

    φ = (indices[1] - 1) * δ
    θ = (indices[2] - 1) * δ

    x = 1.001 * sin(φ) * cos(θ)
    y = 1.001 * sin(φ) * sin(θ)
    z = 1.001 * cos(φ)

    return Vec3f(x, y, z)
end

function closest_point(point, points)
    x, y, z = normalize(point)
    θ, φ = mod2pi(angle(complex(x, y))), acos(z)

    longitudes, meridians = size(points)
    row = round(Int, (longitudes - 1) * φ / pi + 1)
    column = round(Int, (meridians - 1) * θ / (2 * pi) + 1)

    return points[row, column]
end

function closest_indices(point, points)
    d = 4.0
    indices = CartesianIndex(1, 1)

    for c in CartesianIndices(points)
        current_d = distance(point, points[c])
        if current_d < d
            d = current_d
            indices = c
        end
    end

    return indices
end

function plot_setup!(scene, points, texture)
    mesh!(
        scene,
        Sphere(Point3f(0), 1.0),
        color=texture,
        colormap=:twilight,
        colorrange=(0.0, 1.0),
        shading=FastShading,
    )

    camera = Camera3D(
        scene.scene,
        eyeposition=Vec3f(0, 0, -3),
        upvector=Vec3f(0, -1, 0),
        fixed_axis=false,
        center=false,
    )

    focus_vector = @lift(normalize($(camera.eyeposition)))

    focus_point = @lift(closest_point($focus_vector, $points))

    mark_point = Observable(Point(0, 1))

    mark_vector = @lift(indices_to_vector(closest_indices($mark_point, $points), size($points)[1]))

    scatter!(scene, focus_vector, color=:red)
    scatter!(scene, mark_vector, color=:blue)

    return View3D(scene, camera, points, texture, focus_vector, focus_point, mark_vector, mark_point)
end

function init_view3D(longitudes, figure, f, c, state)
    meridians = 2 * longitudes - 1

    scene = LScene(figure, show_axis=false)
    points = Observable(Matrix{Point}(undef, longitudes, meridians))

    texture = Observable(Matrix{Float64}(undef, longitudes, meridians))
    on(points) do points
        update_texture!(texture, points, f, c, state)
    end

    reset_points!(points)
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

function update_points!(points, mobius)
    for c in CartesianIndices(points[])
        points[][c] .= normalize(mobius * points[][c])
    end

    notify(points)
end

function put_menu!(figure, view)
    menu = GridLayout(
        figure,
        width=Relative(0.95),
        valign=0.97,
        tellheight=false,
        tellwidth=false,
        default_colgap=8,
    )

    zoominbutton = Button(menu[1, 1], label="⌕₊", width=30)

    on(zoominbutton.clicks) do _
        mobius = antipodal_hyperbolic(view.focus_point[], 0.5)
        update_points!(view.points, mobius)
    end

    zoomoutbutton = Button(menu[1, 2], label="⌕₋", width=30)

    on(zoomoutbutton.clicks) do _
        mobius = antipodal_hyperbolic(view.focus_point[], 2)
        update_points!(view.points, mobius)
    end

    resetbutton = Button(menu[1, 3], label="↺", width=30)

    on(resetbutton.clicks) do _
        reset_points!(view.points)
        update_cam!(view.scene.scene, Vec3f(0, 0, -3), Vec3f(0, 0, 0), Vec3f(0, -1, 0))
    end

    markbutton = Button(menu[1, 4], label="⋅", width=30)

    on(markbutton.clicks) do _
        view.mark_point[] .= closest_point(view.focus_vector[], view.points[])
        notify(view.mark_point)
    end

    return markbutton
end

function Viewer3D(f::Function; crit=0.0im, c=0.0im)
    rational_map = RationalMap(f, crit)
    state = State(1e-4, 200)
    figure = Figure(size=(1000, 560))

    # Initialize Mandel and Julia Views
    meridian = 501
    mandel = init_view3D(meridian, figure[1, 1], rational_map.f, rational_map.crit, state)
    julia = init_view3D(meridian, figure[1, 2], rational_map.f, mandel.mark_point, state)

    rowsize!(figure.layout, 1, Aspect(1, 1))
    colgap!(figure.layout, 5)

    markbutton = put_menu!(figure[2, 1], mandel)
    put_menu!(figure[2, 2], julia)

    on(markbutton.clicks) do _
        notify(julia.points)
    end

    return Viewer3D(rational_map, figure, state, mandel, julia)
end

Base.show(io::IO, viewer::Viewer3D) = display(GLMakie.Screen(), viewer.figure)