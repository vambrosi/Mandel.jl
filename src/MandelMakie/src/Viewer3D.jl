using GeometryBasics: Tesselation, uv_normal_mesh

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

    # mark_vector::Observable{Vec3f}
    mark_point::Observable{Point}

    path_points::Observable{Vector{Point}}
    path_vectors::Observable{Vector{Vec3f}}
    path_length::Observable{Int}
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

    x = sin(φ) * cos(θ)
    y = sin(φ) * sin(θ)
    z = cos(φ)

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

function arc(v0, v1, segments)
    v_start, v_end = normalize(v0), normalize(v1)
    arc_steps = Vector{Vec3f}(undef, segments)
    arc_steps[1] = v_start

    # In case vectors are antipodal, pick a random great circle.
    if norm(v_start + v_end) < 0.01
        # Pick random vector
        u = v_start
        while norm(u - v_start) < 0.01 || norm(u - v_end) < 0.01
            u = rand(0.0:0.1:1.0, 3)
        end

        # Make it orthonormal to v_start
        u = normalize(u - (u ⋅ v_start) * v_start)

        # Pick first vector to be 1 / segments away along the great circle
        t = pi / segments
        arc_steps[2] = cos(t) * v_start + sin(t) * u

    # Otherwise pick the first point along the arc containing v_start and v_end
    else
        # First vector after v_start (slerp formula)
        ϕ = acos(clamp(v_start ⋅ v_end/(norm(v_start) * norm(v_end)), -1, 1))
        t = 1 / segments # how far to walk along the arc
        arc_steps[2] = 1 / sin(ϕ) * (sin((1 - t) * ϕ) * v_start + sin(t * ϕ) * v_end)
    end

    # Other vectors (by reflecting pt_{i-2} across pt_{i-1} along the arc)
    dot_product = 2 * arc_steps[1] ⋅ arc_steps[2]
    for i in 3:segments
        arc_steps[i] = dot_product * arc_steps[i-1] - arc_steps[i-2]
    end

    return arc_steps
end

function plot_setup!(scene, points, texture, is_mandel, f, c)
    msh = uv_normal_mesh(Tesselation(Sphere(Point3f(0), 1f0), 50))
    mesh!(
        scene,
        msh,
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
    # mark_vector = @lift(indices_to_vector(closest_indices($mark_point, $points), size($points)[1]))

    path_length = is_mandel ? Observable(1) : Observable(1)
    path_points = is_mandel ? @lift([$mark_point]) : @lift(orbit(f, $mark_point, $c, $path_length))
    path_vectors = @lift begin
        vectors = Vector{Vec3f}(undef, $path_length)
        for i in 1:$path_length
            @inbounds vectors[i] = 1.01 .* indices_to_vector(closest_indices($path_points[i], $points), size($points)[1])
        end
        return vectors
    end

    trace = @lift begin
        segments = 20
        trace = Vector{Vec3f}(undef, ($path_length - 1) * segments + 1)
        for i in 2:$path_length
            v1 = $path_vectors[i - 1]
            v2 = $path_vectors[i]
            trace[(i-2) * segments + 1 : (i-1) * segments] = 1.01 .* arc(v1, v2, segments)
        end

        trace[end] = $path_vectors[end]
        return trace
    end

    scatter!(scene, focus_vector, color=:red)
    scatter!(scene, path_vectors, color=:blue, markersize = 12)
    lines!(scene, trace, color=:blue, linewidth = 2)

    return View3D(
        scene, camera, points, texture, focus_vector, focus_point,
        mark_point, path_points, path_vectors, path_length
    )
end

function init_view3D(longitudes, figure, f, c, state, is_mandel)
    meridians = 2 * longitudes - 1

    scene = LScene(figure, show_axis=false)
    points = Observable(Matrix{Point}(undef, longitudes, meridians))

    texture = Observable(Matrix{Float64}(undef, longitudes, meridians))
    on(points) do points
        update_texture!(texture, points, f, c, state)
    end

    reset_points!(points)
    view = plot_setup!(scene, points, texture, is_mandel, f, c)

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

    return markbutton, menu
end

function Viewer3D(f::Function; crit=0.0im, c=0.0im, longitudes=501)
    rational_map = RationalMap(f, crit)
    state = State(1e-4, 200)
    figure = Figure(size=(1000, 560))

    # Initialize Mandel and Julia Views
    longitudes = longitudes % 2 == 1 ? longitudes : longitudes + 1 # has to be an odd number
    mandel = init_view3D(longitudes, figure[1, 1], rational_map.f, rational_map.crit, state, true)
    julia = init_view3D(longitudes, figure[1, 2], rational_map.f, mandel.mark_point, state, false)

    rowsize!(figure.layout, 1, Aspect(1, 1))
    colgap!(figure.layout, 5)

    markbutton, _ = put_menu!(figure[2, 1], mandel)
    _, menu = put_menu!(figure[2, 2], julia)

    on(markbutton.clicks) do _
        notify(julia.points)
    end

    Label(menu[1, 5], "Iterates:")

    orbit_length_input = Textbox(
        menu[1, 6],
        width = 50,
        placeholder = string(julia.path_length[] - 1),
        validator = Int,
    )

    on(orbit_length_input.stored_string) do s
        julia.path_length[] = max(parse(Int, s), 0) + 1
    end

    return Viewer3D(rational_map, figure, state, mandel, julia)
end

Base.show(io::IO, viewer::Viewer3D) = display(GLMakie.Screen(), viewer.figure)