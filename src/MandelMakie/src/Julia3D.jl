using GeometryBasics: mesh as meshGB
using GeometryBasics: coordinates
using ChunkSplitters

id_mobius = Mobius(1, 0, 0, 1)

struct Julia3D
    f::Function
    f_projective::Function

    figure::Figure
    scene::LScene
    camera::Camera3D

    mobius::Observable{Mobius}
    focus::Observable{Point}
    vertex_colors::Observable{Vector{Float64}}
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

function update_colors!(vertex_colors, vertices, f, c, mobius)
    Threads.@threads for i in eachindex(vertices)
        pt = vector_to_point(normalize(vertices[i]); mobius=mobius)
        @inbounds vertex_colors[][i] = multiplier(f, pt, c, 1e-4, 200)
    end

    # Threads.@threads for indexes in chunks(vertices; n=Threads.nthreads())
    #     for i in indexes
    #         pt = vector_to_point(normalize(vertices[i]), mobius=mobius)
    #         @inbounds vertex_colors[][i] = multiplier(f, pt, c, 1e-4, 200)
    #     end
    # end
    notify(vertex_colors)
end

function Julia3D(f::Function, c::Number; show_critical_points=false)
    f_proj = to_point_family(f)
    crit_pts = critical_points(f, c)

    figure = Figure(size=(600, 650))
    scene = LScene(figure[1, 1], show_axis=false)
    _mesh = meshGB(Tesselation(Sphere(Point3f(0), 1), 800))

    rowsize!(figure.layout, 1, Aspect(1, 1))
    colgap!(figure.layout, 5)

    vertices = coordinates(_mesh)
    vertex_colors = Observable(similar(vertices, Float64))
    mobius = Observable(id_mobius)
    inv_mobius = @lift(inv($mobius))

    on(mobius) do M
        update_colors!(vertex_colors, vertices, f_proj, Point(c, 1), M)
    end
    notify(mobius)

    mesh!(
        scene,
        _mesh,
        color=vertex_colors,
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
    focus_point = @lift(vector_to_point($focus_vector; mobius=$mobius))
    scatter!(scene, focus_vector, color=1, colormap=:tab10, colorrange = (1, 10))

    if show_critical_points
        critical_vectors = @lift(1.005 .* point_to_vector.(crit_pts; mobius=$inv_mobius))
        scatter!(scene, critical_vectors, color=2, colormap=:tab10, colorrange = (1, 10))
    end

    post_critical_orbit_points = []
    post_critical_orbit_vectors = []
    post_critical_orbit_traces = []
    for crit_pt in crit_pts
        orbit_points = attracting_orbit(f_proj, crit_pt, Point(c, 1), 1e-4, 500)
        push!(post_critical_orbit_points, orbit_points)

        orbit_vectors = @lift(1.005 .* point_to_vector.(orbit_points; mobius=$inv_mobius))
        push!(post_critical_orbit_vectors, orbit_vectors)

        trace = @lift begin
            segments = 20
            path_length = length(orbit_points)
            trace = Vector{Vec3f}(undef, (path_length - 1) * segments + 1)
            for i in 2:path_length
                v1 = $orbit_vectors[i - 1]
                v2 = $orbit_vectors[i]
                trace[(i-2) * segments + 1 : (i-1) * segments] = 1.01 .* arc(v1, v2, segments)
            end

            trace[end] = $orbit_vectors[end]
            return trace
        end
        push!(post_critical_orbit_traces, trace)
    end

    for (i, orbit) in enumerate(post_critical_orbit_vectors)
        scatter!(scene, orbit, color = 2 + mod1(2*i, 8), colormap=:tab10, colorrange = (1, 10))
    end

    for (i, trace) in enumerate(post_critical_orbit_traces)
        lines!(scene, trace, color = 2 + mod1(2*i, 8), colormap=:tab10, colorrange = (1, 10), linewidth = 2)
    end

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
        mobius[] = antipodal_hyperbolic(focus_point[], 0.5) * mobius[]
    end

    zoomoutbutton = Button(menu[1, 2], label="⌕₋", width=30)

    on(zoomoutbutton.clicks) do _
        mobius[] = antipodal_hyperbolic(focus_point[], 2.0) * mobius[]
    end

    resetbutton = Button(menu[1, 3], label="↺", width=30)

    on(resetbutton.clicks) do _
        mobius[] = id_mobius
        update_cam!(scene.scene, Vec3f(0, 0, -3), Vec3f(0, 0, 0), Vec3f(0, -1, 0))
    end

    return Julia3D(f, f_proj, figure, scene, camera, mobius, focus_point, vertex_colors)
end

Base.show(io::IO, viewer::Julia3D) = display(GLMakie.Screen(), viewer.figure)