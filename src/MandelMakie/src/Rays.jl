module Rays
using LinearAlgebra
using PolynomialRoots

export PolynomialWrapper, compute_ray!
# Constants
eps_land = 1e-6
eps_coland = 1e-3
iterations = 1000
outer_radius = 10
quality = 1

struct PolynomialWrapper
    polynomial::Vector{ComplexF64}
    stored_psi_values::Dict{Rational{Int64},Vector{ComplexF64}}
    degree::Int

    function PolynomialWrapper(polynomial)
        degree = length(polynomial) - 1
        if degree == -1
            degree = 0
        end
        new(polynomial, Dict{Rational{Int64},Vector{ComplexF64}}(), degree)
    end
end

wrapper = PolynomialWrapper([im, 0, 1])

function σ(wrapper::PolynomialWrapper, q::Rational{Int64})
    ret = q * wrapper.degree
    return ret - floor(Int, ret)
end

@assert σ(wrapper, 1 // 7) == 2 // 7

function newton(wrapper::PolynomialWrapper, z::ComplexF64, guess::ComplexF64)
    p = copy(wrapper.polynomial)
    p[1] -= z
    # return minimum([w for w in roots(p) | abs(guess-w[1])])
    lroots = roots(p)
    return lroots[argmin(abs.(guess .- lroots))]
end

@assert newton(wrapper, 7 + 0.0im, 7 + 0.0im) == 2.6524580874978474 - 0.18850439234335525im

function compute_ray!(wrapper::PolynomialWrapper, θout::Rational{Int64})
    if haskey(wrapper.stored_psi_values, θout)
        return
    end
    θi = θout
    angle_list = Rational{Int64}[]
    while !in(θi, angle_list)
        push!(angle_list, θi)
        θi = σ(wrapper, θi)
    end

    if θi != θout  # pre-periodic
        compute_ray!(wrapper, σ(wrapper, θout))
        wrapper.stored_psi_values[θout] = ComplexF64[]
        for j in reverse(1:quality)
            r = 100 * wrapper.degree^(j / quality)
            num = r * exp(2π * θout * im)
            push!(wrapper.stored_psi_values[θout], num)
        end
        while true
            guess = last(wrapper.stored_psi_values[θout])
            depth = length(wrapper.stored_psi_values[θout]) - quality + 1
            if depth ≥ length(wrapper.stored_psi_values[σ(wrapper, θout)])
                return
            end
            image = wrapper.stored_psi_values[σ(wrapper, θout)][depth]
            new = newton(wrapper, image, guess)
            push!(wrapper.stored_psi_values[θout], new)
        end
    end

    # assume periodic and not already in thing
    for θ in angle_list
        wrapper.stored_psi_values[θ] = ComplexF64[]
    end
    for j in reverse(1:quality)
        r = 100 * wrapper.degree^(j / quality)
        for θ in angle_list
            num = r * exp(2π * θ * im)
            push!(wrapper.stored_psi_values[θ], num)
        end
    end
    while true
        for θ in angle_list
            guess = last(wrapper.stored_psi_values[θ])
            image = wrapper.stored_psi_values[σ(wrapper, θ)][end-quality+1]
            new = newton(wrapper, image, guess)
            push!(wrapper.stored_psi_values[θ], new)
            if abs(new - guess) < eps_land
                return
            end
        end
    end
end

function landing_point(wrapper::PolynomialWrapper, θ::Rational{Int64})
    compute_ray!(wrapper, θ)
    return last(wrapper.stored_psi_values[θ])
end

function co_land(
    wrapper::PolynomialWrapper,
    laminations::Vector{Vector{Rational{Int64}}},
)::Bool
    for eqclass in laminations
        if length(eqclass) < 2
            continue
        end
        p = landing_point(wrapper, eqclass[1])
        for θ in eqclass
            if abs(p - landing_point(wrapper, θ)) > eps_coland
                return false
            end
        end
    end
    return true
end

@assert co_land(wrapper, [[1 // 7, 2 // 7, 4 // 7]])
@assert !co_land(wrapper, [[1 // 7, 1 // 14]])

function auto_rays(coefficients::Vector{ComplexF64}, periods)
    w = PolynomialWrapper(coefficients)
    lamination = []
    d = w.degree
    for period in 1:maximum(periods)
        denom = d^period - 1
        classes = []
        for numeraitor in 0:(denom-1)
            if denom > 1 && numeraitor == 0
                continue
            end
            if gcd(denom,numeraitor) > 1
                continue
            end

            angle = numeraitor//denom
            added = false
            compute_ray!(w,angle)
            for (i, class) in enumerate(classes)
                if co_land(w,[[class[1], angle]]) 
                    push!(classes[i],angle)
                    added = true
                    break
                end
            end
            if !added 
                push!(classes,[angle])
            end
        end
        for class in classes 
            if length(class) > 1
                push!(lamination, class)
            end
        end
    end
    angles = []
    println(lamination)
    relivant_rays = reduce(vcat, lamination, init=[])

    return [w.stored_psi_values[angle] for angle in relivant_rays]
end

function rays(coefficients::Vector{ComplexF64}, angles)
    w = PolynomialWrapper(coefficients)
    for a in angles
        compute_ray!(w, a)
    end

    return collect(values(w.stored_psi_values))
end
end
