module MandelGTK

export viewer

using Gtk4, Graphics, Cairo, Colors, FixedPointNumbers

include("ColorSchemes.jl")

function iterate(f, z, c, max_iter)
    for i in 1:max_iter
        z = f(z, c)

        if abs(z) > 10.0
            return i, abs(z)
        end
    end

    return nothing, abs(z)
end

color(::Nothing, ::Any) = twilight_RGB24[500]

function color(iterations::Integer, r::Real)
    height = (iterations - log(2, log(r))) / 50.0
    index = mod1(round(Int, 999 * mod(height, 1.0)), 1000)
    return twilight_RGB24[index]
end


function viewer(f::Function, c::Number)
    win = GtkWindow("Mandel", 500, 500)
    
    canvas = GtkCanvas()
    push!(win, canvas)

    @guarded draw(canvas) do widget
        ctx = getgc(canvas)
        h = height(canvas)
        w = width(canvas)

        center = 0.0im
        diameter = 4.0
        upp = diameter / min(h, w)

        m2 = Matrix{RGB24}(undef, w, h)

        for i in 1:w
            for j in 1:h
                x = (i + 0.5 - w / 2) * upp
                y = - (j + 0.5 - h / 2) * upp
                z = center + complex(x, y)

                m2[i, j] = color(iterate(f, z, c, 200)...)
            end
        end

        plt = CairoImageSurface(m2)

        set_source_surface(ctx, plt, 0, 0)
        paint(ctx)
    end

    show(win)
    
    if !isinteractive()
        signal_connect(win, :close_request) do widget
            GC.gc()
        end
        @async Gtk4.GLib.glib_main()
    end

    return nothing
end

end