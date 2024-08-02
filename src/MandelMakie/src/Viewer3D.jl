using GeometryBasics: Tesselation, uv_normal_mesh, coordinates
using GeometryBasics: mesh as meshGB

const id_mobius = Mobius(1, 0, 0, 1)

mutable struct State
    ε::Float64
    max_iter::Int
end

struct View3D
    scene::LScene
    camera::Camera3D

    vertex_colors::Observable{Vector{Float64}}
    mobius::Observable{Mobius}

    focus::Observable{Point}
    focus_antipode::Observable{Point}
    mark::Observable{Point}

    path::Observable{Vector{Point}}
    path_length::Observable{Int}
end

function vector_to_point(v; mobius::Mobius = id_mobius)
    # Assumes v is in the 2-sphere
    isapprox(v, Vec3f(0, 0, 1)) && return Point(1, 0)

    θ, φ = mod2pi(angle(complex(v[1], v[2]))), acos(v[3])
    pt = mobius * Point(sin(φ) * exp(im * θ), 1 - cos(φ))

    # Returns a point in the 3-sphere
    return normalize(pt)
end

function point_to_vector(pt::Point; mobius::Mobius = id_mobius)
    new_pt = mobius * pt
    w = 2 * new_pt[1] * conj(new_pt[2])
    z = abs2(new_pt[1]) - abs2(new_pt[2])
    v = Vec3f(real(w), imag(w), z)

    # Returns a point in the 2-sphere
    return normalize(v)
end

function orthogonal_vector(v::Vec3f)
    w = v + Vec3f(norm(v) * sign(v[1]), 0, 0)
    B = I - 2 * w * w' / (w' * w)
    return B[:, 2]
end

function update_colors!(vertex_colors, vertices, f_proj, c::Point, mobius, state)
    Threads.@threads for i in eachindex(vertices)
        pt = vector_to_point(normalize(vertices[i]); mobius = mobius)
        @inbounds vertex_colors[][i] = multiplier(f_proj, pt, c, state.ε, state.max_iter)
    end
    notify(vertex_colors)
end

function update_colors!(vertex_colors, vertices, f_proj, crit::Function, mobius, state)
    Threads.@threads for i in eachindex(vertices)
        c = vector_to_point(normalize(vertices[i]); mobius = mobius)
        @inbounds vertex_colors[][i] =
            multiplier(f_proj, crit(c), c, state.ε, state.max_iter)
    end
    notify(vertex_colors)
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
        ϕ = acos(clamp(v_start ⋅ v_end / (norm(v_start) * norm(v_end)), -1, 1))
        t = 1 / segments # how far to walk along the arc
        arc_steps[2] = 1 / sin(ϕ) * (sin((1 - t) * ϕ) * v_start + sin(t * ϕ) * v_end)
    end

    # Other vectors (by reflecting pt_{i-2} across pt_{i-1} along the arc)
    dot_product = 2 * arc_steps[1] ⋅ arc_steps[2]
    for i = 3:segments
        arc_steps[i] = dot_product * arc_steps[i-1] - arc_steps[i-2]
    end

    return arc_steps
end

function View3D(figure, f_proj, c, state; is_mandel, mark)
    scene = LScene(figure, show_axis = false)

    _mesh = meshGB(Tesselation(Sphere(Point3f(0), 1), 800))
    vertices = coordinates(_mesh)

    mobius = Observable{Mobius}(id_mobius)
    inv_mobius = @lift(inv($mobius))
    vertex_colors = Observable(similar(vertices, Float64))

    onany(mobius, c) do M, c
        update_colors!(vertex_colors, vertices, f_proj, c, M, state)
    end
    notify(mobius)

    mesh!(
        scene,
        _mesh,
        color = vertex_colors,
        colormap = :twilight,
        colorrange = (0.0, 1.0),
        shading = FastShading,
        inspectable = false,
    )

    camera = Camera3D(
        scene.scene,
        eyeposition = Vec3f(0, 0, -3),
        upvector = Vec3f(0, -1, 0),
        fixed_axis = false,
        center = false,
    )

    focus_vector = @lift(normalize($(camera.eyeposition)))
    focus_point = @lift(vector_to_point($focus_vector; mobius = $mobius))
    focus_antipode_point = @lift(vector_to_point(-$focus_vector; mobius = $mobius))

    mark_point = isfinite(mark) ? Observable(Point(mark, 1)) : Observable(Point(1, 0))

    path_length = Observable{Int}(1)
    path_points =
        is_mandel ? @lift([$mark_point]) :
        @lift(orbit(f_proj, $mark_point, $c, $path_length))
    path_vectors = @lift(1.001 .* point_to_vector.($path_points; mobius = $inv_mobius))

    trace = @lift begin
        segments = 20
        trace = Vector{Vec3f}(undef, ($path_length - 1) * segments + 1)
        for i = 2:$path_length
            v1 = $path_vectors[i-1]
            v2 = $path_vectors[i]
            trace[(i-2)*segments+1:(i-1)*segments] = 1.005 .* arc(v1, v2, segments)
        end

        trace[end] = $path_vectors[end]
        return trace
    end

    scatter!(
        scene,
        focus_vector,
        color = :red,
        overdraw = true,
        marker = '+',
        markersize = 15,
        glowwidth = 1,
        glowcolor = :white,
        inspector_label = (self, i, p) ->
            string(to_complex_plane(vector_to_point(normalize(p); mobius = mobius[]))),
    )
    scatter!(
        scene,
        path_vectors,
        color = :blue,
        glowwidth = 2,
        glowcolor = :white,
        inspector_label = (self, i, p) ->
            string(to_complex_plane(vector_to_point(normalize(p); mobius = mobius[]))),
    )
    lines!(scene, trace, color = :blue, linewidth = 2, inspectable = false)

    return View3D(
        scene,
        camera,
        vertex_colors,
        mobius,
        focus_point,
        focus_antipode_point,
        mark_point,
        path_points,
        path_length,
    )
end

abstract type AbstractViewer3D end

struct Viewer3D <: AbstractViewer3D
    rational_map::RationalMap
    figure::Figure
    state::State

    mandel::View3D
    julia::View3D
end

function put_menu!(figure, view)
    menu = GridLayout(
        figure,
        width = Relative(0.95),
        valign = 0.97,
        tellheight = false,
        tellwidth = false,
        default_colgap = 8,
    )

    zoominbutton = Button(menu[1, 1], label = "⌕₊", width = 30)

    on(zoominbutton.clicks) do _
        view.mobius[] = hyperbolic(view.focus[], view.focus_antipode[], 0.5) * view.mobius[]
    end

    zoomoutbutton = Button(menu[1, 2], label = "⌕₋", width = 30)

    on(zoomoutbutton.clicks) do _
        view.mobius[] = hyperbolic(view.focus[], view.focus_antipode[], 2.0) * view.mobius[]
    end

    resetbutton = Button(menu[1, 3], label = "↺", width = 30)

    on(resetbutton.clicks) do _
        focus = point_to_vector(view.focus[])
        view.mobius[] = id_mobius

        v = orthogonal_vector(focus)
        update_cam!(view.scene.scene, 3 * focus, Vec3f(0, 0, 0), v)
    end

    markbutton = Button(menu[1, 4], label = "⋅", width = 30)

    on(markbutton.clicks) do _
        view.mark[] = view.focus[]
    end

    return markbutton, menu
end

function put_julia_menu!(menu, julia)
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

    return menu
end

function Viewer3D(f::Function; crit = 0.0im, c = 0.0im)
    rational_map = RationalMap(f, crit)
    state = State(1e-4, 200)
    figure = Figure(size = (1000, 560))

    DataInspector(figure)

    # Initialize Mandel and Julia Views
    mandel = View3D(
        figure[1, 1],
        rational_map.f_proj,
        Observable(rational_map.crit),
        state;
        is_mandel = true,
        mark = c,
    )
    julia = View3D(
        figure[1, 2],
        rational_map.f_proj,
        mandel.mark,
        state;
        is_mandel = false,
        mark = 0.0im,
    )

    rowsize!(figure.layout, 1, Aspect(1, 1))
    colgap!(figure.layout, 5)

    markbutton, _ = put_menu!(figure[2, 1], mandel)
    _, menu = put_menu!(figure[2, 2], julia)

    on(markbutton.clicks) do _
        notify(julia.mobius)
    end

    put_julia_menu!(menu, julia)

    return Viewer3D(rational_map, figure, state, mandel, julia)
end

struct Julia3D <: AbstractViewer3D
    rational_map::RationalMap
    figure::Figure
    state::State
    julia::View3D
    reset::Any
end

function Julia3D(f::Function, parameter::Number = 0.0im)
    if hasmethod(f, Tuple{ComplexF64,ComplexF64})
        g = (z, c) -> f(z, c)
    elseif hasmethod(f, ComplexF64)
        g = (z, _) -> f(z)
    else
        throw("Function not defined!")
    end

    rational_map = RationalMap(g, 0.0im)
    state = State(1e-4, 200)
    figure = Figure(size = (500, 560))

    DataInspector(figure)

    julia = View3D(
        figure[1, 1],
        rational_map.f_proj,
        Observable(Point(parameter, 1)),
        state;
        is_mandel = false,
        mark = 0.0im,
    )

    rowsize!(figure.layout, 1, Aspect(1, 1))
    colgap!(figure.layout, 5)

    _, menu = put_menu!(figure[2, 1], julia)
    put_julia_menu!(menu, julia)

    return Julia3D(rational_map, figure, state, julia, contents(menu[1, 3])[1])
end

Base.show(io::IO, viewer::AbstractViewer3D) = display(GLMakie.Screen(), viewer.figure)

function set_parameter!(viewer::Viewer3D, c)
    viewer.mandel.mark[] = isfinite(c) ? Point(c, 1) : Point(1, 0)
    notify(viewer.julia.mobius)
    return nothing
end
