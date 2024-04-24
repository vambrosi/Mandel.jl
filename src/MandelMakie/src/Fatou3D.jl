using GLMakie.Colors: LCHab

struct Fatou3D <: AbstractViewer3D
    f::Function
    f_projective::Function

    figure::Figure
    scene::LScene
    camera::Camera3D

    mobius::Observable{Mobius}
    focus::Observable{Point}

    attractors::Vector{Vector{Point}}
end

function hue_from_fatou(component, n)
    if n == 2
        return component == 1 ? 30 : 270
    else
       return mod(360 * (component - 1) / n + 30, 360)
    end
end

function color_from_fatou(fatou::FatouIterationDistance, n_fatou::Int, max_iter::Int)
    fatou.component_index == 0 && return LCHab(30, 0, 0)

    t = mod(fatou.preperiod / (fatou.period * 16), 2.0)
    hue = hue_from_fatou(fatou.component_index, n_fatou)

    return LCHab(20 * cospi(t) + 50, 40, hue)
end

function update_fatou!(vertex_colors, vertices, f, c, mobius, post_critical_orbit_points, n_fatou)
    Threads.@threads for i in eachindex(vertices)
        pt = vector_to_point(normalize(vertices[i]); mobius=mobius)
        fatou = convergence_time(f, pt, c, post_critical_orbit_points, 1e-4, 200)
        @inbounds vertex_colors[][i] = color_from_fatou(fatou, n_fatou, 200)
    end
    notify(vertex_colors)
end

function nearby_up_to_shift(orbit1, orbit2, ε)
    # Since both are orbits, only need to check for one point
    pt1 = orbit1[1]

    for pt2 in orbit2
        distance(pt1, pt2) < ε && return true
    end

    return false
end

function get_attractors(f_proj, crit_pts, c)
    attractors = []
    for crit_pt in crit_pts
        orbit = attracting_orbit(f_proj, crit_pt, Point(c, 1), 1e-4, 500)
        all([!nearby_up_to_shift(orbit, other_orbit, 1e-4) for other_orbit in attractors]) &&
            push!(attractors, orbit)
    end
    return attractors
end

function Fatou3D(f::Function, c::Number=0.0im; show_critical_points=false)
    if hasmethod(f, Tuple{ComplexF64, ComplexF64})
        g = (z, c) -> f(z, c)
    elseif hasmethod(f, ComplexF64)
        g = (z, _) -> f(z)
    else
        throw("Function not defined!")
    end

    f_proj = to_point_family(g)
    crit_pts = critical_points(g, c)

    figure = Figure(size=(600, 650))
    scene = LScene(figure[1, 1], show_axis=false)

    rowsize!(figure.layout, 1, Aspect(1, 1))
    colgap!(figure.layout, 5)

    mobius = Observable(id_mobius)
    inv_mobius = @lift(inv($mobius))

    fatou_mesh = meshGB(Tesselation(Sphere(Point3f(0), 1), 800))
    fatou_vertices = coordinates(fatou_mesh)
    fatou_vertex_colors = Observable(similar(fatou_vertices, LCHab))

    mesh!(
        scene,
        fatou_mesh,
        color=fatou_vertex_colors,
        shading=NoShading,
    )

    camera = Camera3D(
        scene.scene,
        eyeposition=Vec3f(0, 0, -3),
        upvector=Vec3f(0, -1, 0),
        fixed_axis=false,
        center=false,
    )

    focus_vector = @lift(normalize($(camera.eyeposition)))
    focus_point = @lift(vector_to_point($focus_vector; mobius=$mobius))
    focus_antipode_point = @lift(vector_to_point(- $focus_vector; mobius=$mobius))
    scatter!(scene, focus_vector, color=:red)

    if show_critical_points
        critical_vectors = @lift(1.001 .* point_to_vector.(crit_pts; mobius=$inv_mobius))
        scatter!(scene, critical_vectors, color=:green)
    end

    attractor_points = get_attractors(f_proj, crit_pts, c)
    attractor_vectors = []
    attractor_traces = []

    for attractor in attractor_points
        vectors = @lift(1.002 .* point_to_vector.(attractor; mobius=$inv_mobius))
        push!(attractor_vectors, vectors)

        trace = @lift begin
            segments = 20
            path_length = length(attractor)
            trace = Vector{Vec3f}(undef, path_length * segments + 1)
            for i in 1:path_length
                start_index, end_index = i, mod1(i+1, path_length)
                v1 = $vectors[start_index]
                v2 = $vectors[end_index]
                trace[(start_index-1) * segments + 1 : (start_index) * segments] = 1.005 .* arc(v1, v2, segments)
            end

            trace[end] = $vectors[1]
            return trace
        end
        push!(attractor_traces, trace)
    end

    n_fatou = length(attractor_vectors)
    for (i, orbit) in enumerate(attractor_vectors)
        scatter!(scene, orbit, color = LCHab(40, 40, hue_from_fatou(i, n_fatou)), glowwidth=2, glowcolor=:white)
    end

    for (i, trace) in enumerate(attractor_traces)
        lines!(scene, trace, color = LCHab(40, 40, hue_from_fatou(i, n_fatou)), linewidth=2)
    end

    on(mobius) do M
        update_fatou!(fatou_vertex_colors, fatou_vertices, f_proj, Point(c, 1), M, attractor_points, n_fatou)
    end
    notify(mobius)

    menu = GridLayout(
        figure[2, 1],
        width=Relative(0.95),
        valign=0.97,
        tellheight=false,
        tellwidth=false,
        default_colgap=8,
    )

    zoominbutton = Button(menu[1, 1], label="⌕₊", width=30)

    on(zoominbutton.clicks) do _
        mobius[] = hyperbolic(focus_point[], focus_antipode_point[], 0.5) * mobius[]
    end

    zoomoutbutton = Button(menu[1, 2], label="⌕₋", width=30)

    on(zoomoutbutton.clicks) do _
        mobius[] = hyperbolic(focus_point[], focus_antipode_point[], 2.0) * mobius[]
    end

    resetbutton = Button(menu[1, 3], label="↺", width=30)

    on(resetbutton.clicks) do _
        focus_vector = point_to_vector(focus_point[])
        mobius[] = id_mobius
        update_cam!(scene.scene, 3 * focus_vector, Vec3f(0, 0, 0))
    end

    return Fatou3D(f, f_proj, figure, scene, camera, mobius, focus_point, attractor_points)
end