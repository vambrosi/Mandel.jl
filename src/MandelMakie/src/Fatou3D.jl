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
    reset::Any
end

grad = convert.(LCHab, cgrad(:twilight))

function get_color(index, n)
    n == 1 && return grad[255].c, grad[255].h
    n > 8 && return 40, mod(360 * (index - 1) / n + 30, 360)

    #RGB CA736C 202 115 108 red
    index == 1 && return 38.35791098655309, 29.49442824523494
    #RGB 5794D0 87 148 208 blue
    index == 2 && return 36.822489073365745, 266.89189631383687
    #RGB 47A477 71 164 119 green
    index == 3 && return 41.31205868565174, 158.3817555910295
    #RGB 8D9741 141 151 65 yellow
    index == 4 && return 46.2028504263141, 110.51854840517339
    #RGB 00A2AF 0 162 175 cyan
    index == 5 && return 34.439876883846594, 209.09162086774987
    #RGB BC73A4 188 115 164 magenta
    index == 6 && return 38.49553456026353, 338.47504267487597
    #RGB BA823A 186 130 58 orange
    index == 7 && return 48.58788683488929, 72.5301368954557
    #RGB 9481CC 148 129 204 purple
    index == 8 && return 43.17999997459473, 302.8931595387337
end

function color_from_fatou(fatou::FatouIterationDistance, n_fatou::Int)
    if n_fatou == 1
        fatou.component_index == 0 && return grad[255]
        depth = mod1(16 * fatou.preperiod ÷ fatou.period, 510)
        return grad[depth]
    end

    fatou.component_index == 0 && return LCHab(30, 0, 0)

    t = mod(fatou.preperiod / (fatou.period * 16), 2.0)
    chroma, hue = get_color(fatou.component_index, n_fatou)

    return LCHab(20 * cospi(t) + 50, chroma, hue)
end

function update_fatou!(
    vertex_colors,
    vertices,
    f,
    c,
    mobius,
    post_critical_orbit_points,
    n_fatou,
)
    Threads.@threads for i in eachindex(vertices)
        pt = vector_to_point(normalize(vertices[i]); mobius = mobius)
        fatou = convergence_time(f, pt, c, post_critical_orbit_points, 1e-4, 200)
        @inbounds vertex_colors[][i] = color_from_fatou(fatou, n_fatou)
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

function find_attractors(f_proj, crit_pts, c)
    attractors = []
    for crit_pt in crit_pts
        orbit = attracting_orbit(f_proj, crit_pt, Point(c, 1), 1e-4, 500)
        all([
            !nearby_up_to_shift(orbit, other_orbit, 1e-4) for other_orbit in attractors
        ]) && push!(attractors, orbit)
    end
    return attractors
end

function Fatou3D(f::Function, c::Number = 0.0im; show_critical_points = false)
    if hasmethod(f, Tuple{ComplexF64,ComplexF64})
        g = (z, c) -> f(z, c)
    elseif hasmethod(f, ComplexF64)
        g = (z, _) -> f(z)
    else
        throw("Function not defined!")
    end

    f_proj = to_point_family(g)
    crit_pts = get_critical_points(g, c)

    figure = Figure(size = (600, 650))
    scene = LScene(figure[1, 1], show_axis = false)

    DataInspector(figure)

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
        color = fatou_vertex_colors,
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

    if show_critical_points
        critical_vectors = @lift(1.001 .* point_to_vector.(crit_pts; mobius = $inv_mobius))
        scatter!(
            scene,
            critical_vectors,
            color = :green,
            glowwidth = 2,
            glowcolor = :white,
            inspector_label = (self, i, p) ->
                string(to_complex_plane(vector_to_point(normalize(p); mobius = mobius[]))),
        )
    end

    attractor_points = find_attractors(f_proj, crit_pts, c)
    attractor_vectors = []
    attractor_traces = []

    for attractor in attractor_points
        vectors = @lift(1.002 .* point_to_vector.(attractor; mobius = $inv_mobius))
        push!(attractor_vectors, vectors)

        trace = @lift begin
            segments = 20
            path_length = length(attractor)
            trace = Vector{Vec3f}(undef, path_length * segments + 1)
            for i = 1:path_length
                start_index, end_index = i, mod1(i + 1, path_length)
                v1 = $vectors[start_index]
                v2 = $vectors[end_index]
                trace[(start_index-1)*segments+1:(start_index)*segments] =
                    1.005 .* arc(v1, v2, segments)
            end

            !isempty($vectors) && (trace[end] = $vectors[1])
            return trace
        end
        push!(attractor_traces, trace)
    end

    n_fatou = length(attractor_vectors)
    for (i, orbit) in enumerate(attractor_vectors)
        chroma, hue = get_color(i, n_fatou)
        scatter!(
            scene,
            orbit,
            color = LCHab(40, chroma, hue),
            glowwidth = 2,
            glowcolor = :white,
            inspector_label = (self, i, p) ->
                string(to_complex_plane(vector_to_point(normalize(p); mobius = mobius[]))),
        )
    end

    for (i, trace) in enumerate(attractor_traces)
        chroma, hue = get_color(i, n_fatou)
        lines!(
            scene,
            trace,
            color = LCHab(40, chroma, hue),
            linewidth = 2,
            inspectable = false,
        )
    end

    on(mobius) do M
        update_fatou!(
            fatou_vertex_colors,
            fatou_vertices,
            f_proj,
            Point(c, 1),
            M,
            attractor_points,
            n_fatou,
        )
    end
    notify(mobius)

    menu = GridLayout(
        figure[2, 1],
        width = Relative(0.95),
        valign = 0.97,
        tellheight = false,
        tellwidth = false,
        default_colgap = 8,
    )

    zoominbutton = Button(menu[1, 1], label = "⌕₊", width = 30)

    on(zoominbutton.clicks) do _
        mobius[] = hyperbolic(focus_point[], focus_antipode_point[], 0.5) * mobius[]
    end

    zoomoutbutton = Button(menu[1, 2], label = "⌕₋", width = 30)

    on(zoomoutbutton.clicks) do _
        mobius[] = hyperbolic(focus_point[], focus_antipode_point[], 2.0) * mobius[]
    end

    resetbutton = Button(menu[1, 3], label = "↺", width = 30)

    on(resetbutton.clicks) do _
        f = 3 * point_to_vector(focus_point[])
        mobius[] = id_mobius

        v = orthogonal_vector(f)
        update_cam!(scene.scene, f, Vec3f(0, 0, 0), v)
    end

    return Fatou3D(
        f,
        f_proj,
        figure,
        scene,
        camera,
        mobius,
        focus_point,
        attractor_points,
        resetbutton,
    )
end
