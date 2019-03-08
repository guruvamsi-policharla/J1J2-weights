using Distributions
using LinearAlgebra
include("skyrm_aux.jl")

function sample_gauss_var(v,var)
    #Input is the central vector around which we flip
    v_new = [rand(Normal(v[1],var)),rand(Normal(v[2],var)),rand(Normal(v[3],var))]
    v_new = v_new/sqrt(dot(v_new,v_new))
    return convert(Vector,v_new)
end

function initialise_nail(M::Int, N::Int)
    lat = Array{Vector{Float64},2}(undef,M, N);
    var = 0.000000002
    for i = 1:M
        for j = 1:N
            if (mod(i+j,2) == 0)
                lat[i,j] = sample_gauss_var([0,0,1],var)
            elseif (mod(i+j,2) == 1)
                lat[i,j] = sample_gauss_var([0,0,-1],var)
            end
        end
    end
    return lat
end

function initialise_stripe(M::Int, N::Int)
    lat = Array{Vector{Float64},2}(undef,M, N);
    var = 0.000000002
    for i = 1:M
        for j = 1:N
            if (mod(i,2) == 0)
                lat[i,j] = sample_gauss_var([0,0,1],var)
            elseif (mod(i,2) == 1)
                lat[i,j] = sample_gauss_var([0,0,-1],var)
            end
        end
    end
    return lat
end

function skyrm_hist_nail(size)
    skyrm_arr = Array{Int64}(undef, size)
        for i in 1:size
            lat = initialise_nail(4,4)
            skyrm_arr[i] = abs(round(skyrmion_number(lat)))
        end
return skyrm_arr
end

function skyrm_hist_stripe(size)
    skyrm_arr = Array{Int64}(undef, size)
        for i in 1:size
            lat = initialise_stripe(4,4)
            skyrm_arr[i] = abs(round(skyrmion_number(lat)))
        end
return skyrm_arr
end
