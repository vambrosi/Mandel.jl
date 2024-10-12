using GLMakie

function height(c, R, N)
    z = 0

    for i in 1:N
        if abs(z) > R
            return i / 3 - log(2, log(abs(z)))
        end

        z = (z^2 - 1 - c + c^3) / (z^2 - c^2)
    end

    return NaN
end

function iteration_matrix(center, radius, pixels)
    pixel_length = 2 * radius / pixels
    corner = center - (radius + pixel_length / 2) * complex(1, 1)

    heights = Matrix{Float64}(undef, pixels, pixels)
    for j in 1:pixels
        for i in 1:pixels
            c = corner + pixel_length * complex(i, j)
            heights[i, j] = height(c, 1000.0, 1000) / 20.0

            isfinite(heights[i, j]) || (heights[i, j] = 0.5)
        end
    end

    return heights
end

function Δ(f)
    Δf = zero(f) 
    for y = 2:size(f,2)-1, x = 2:size(f,1)-1
        Δf[x,y] = f[x-1,y] + f[x+1,y] + f[x,y-1] + f[x,y+1] - 4f[x,y]
    end
    return Δf
end

function pool(A, a)
    B = zero(A)
    for y = 2:size(A,2)-1
        for x = 2:size(A,1)-1
            if maximum(abs, A[x-1:x+1, y-1:y+1] .- A[x, y]) > a
                B[x,y] = 0.5
            end
        end
    end
    return B
end

threshold(x, a) = x > a ? 0.5 : 0.0

function plot_matrix(m)
    d = size(m, 1)

    fig = Figure(size = (d, d))
    ax = Axis(fig[1,1], aspect = AxisAspect(1))

    heatmap!(ax, m, colormap = :twilight, colorrange = (0.0, 1.0))

    hidedecorations!(ax)
    hidespines!(ax)

    return fig
end

fig1 = plot_matrix(mod.(iteration_matrix(0.25, 3.0, 3000), 1.0));
save(joinpath("plots", "continuous.png"), fig1)

levels1 = threshold.(mod.(Δ(iteration_matrix(0.25, 3.0, 3000)), 1.0), 0.001);
fig2 = plot_matrix(levels1);
save(joinpath("plots", "laplacian.png"), fig2)

levels2 = pool(iteration_matrix(0.25, 3.0, 3000), 0.05);
fig3 = plot_matrix(levels2);
save(joinpath("plots", "pooling.png"), fig3)