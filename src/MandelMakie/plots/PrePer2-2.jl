using GLMakie

function height(c, ε, N)
    R = 1 / ε
    t = 1
    z = (-c^2 - 1) / (2 * c)

    for i in 1:N
        w = z * (c + 1)^2 / (z^2 + c)^2
        t = t * w
        x = -(c^2 - 2 * z^2 + 1) / (2 * z^2 + 2 * c)
        z = x

        if abs(t) > R
            return 0.0
        elseif abs(t) < ε
            return 0.25
        end
    end

    return 0.5
end

function iteration_matrix(center, radius, pixels)
    pixel_length = 2 * radius / pixels
    corner = center - (radius + pixel_length / 2) * complex(1, 1)

    heights = Matrix{Float64}(undef, pixels, pixels)
    for j in 1:pixels
        for i in 1:pixels
            c = corner + pixel_length * complex(i, j)
            heights[i, j] = height(c, 1e-5, 100000)
        end
    end

    return heights
end

function plot_matrix(m)
    d = size(m, 1)

    fig = Figure(size = (d, d))
    ax = Axis(fig[1,1], aspect = AxisAspect(1))

    heatmap!(ax, m, colormap = :twilight, colorrange = (0.0, 1.0))

    hidedecorations!(ax)
    hidespines!(ax)

    return fig
end

fig1 = plot_matrix(iteration_matrix(-3.0, 10.0, 3000));
save(joinpath("plots", "preper_2-2_parameter.png"), fig1)
