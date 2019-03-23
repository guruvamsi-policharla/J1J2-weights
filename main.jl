#Includes
using Distributed
addprocs(Sys.CPU_THREADS)
#addprocs(4)
println(nprocs())

#@everywhere include("/home/guru/repos/AntiFerro-J1J2Weights/skyrm_aux.jl")
#@everywhere include("/home/guru/repos/AntiFerro-J1J2Weights/error_aux.jl")
#@everywhere include("/home/guru/repos/AntiFerro-J1J2Weights/energy_aux.jl")
#@everywhere include("/home/guru/repos/AntiFerro-J1J2Weights/lat_aux.jl")

@everywhere include("skyrm_aux.jl")
@everywhere include("error_aux.jl")
@everywhere include("energy_aux.jl")
@everywhere include("lat_aux.jl")

@everywhere using SharedArrays
@everywhere using Distributions
@everywhere using StatsBase
@everywhere using LinearAlgebra
@everywhere using JLD2
using Dates
#driver
Tmin = 0.1
Tchange = 0.1
Tmax = 1
N = 4
Temperature = Tmin:Tchange:Tmax
#J_space = [0,0.25,0.5,0.75,1.0]
J_space = [0.0:0.1:0.3;0.35:0.05:0.65;0.7:0.1:1]

E_temp = SharedArray{Float64,6}(length(Temperature),length(J_space),4,3,2,nprocs()-1)
mag_temp = SharedArray{Float64,6}(length(Temperature),length(J_space),4,3,2,nprocs()-1)
skyrm_temp = SharedArray{Float64,6}(length(Temperature),length(J_space),4,3,2,nprocs()-1)
magbind_temp = SharedArray{Float64,5}(length(Temperature),length(J_space),4,2,nprocs()-1)
skyrmbind_temp = SharedArray{Float64,5}(length(Temperature),length(J_space),4,2,nprocs()-1)

proc_complete = SharedArray{Int,1}(nprocs())

for i in 1:nprocs()
    proc_complete[i] = 0
end
@distributed for i in 2:nprocs()
	E_temp[:,:,:,:,:,i-1],
    skyrm_temp[:,:,:,:,:,i-1],
	mag_temp[:,:,:,:,:,i-1],
	magbind_temp[:,:,:,:,i-1],
	skyrmbind_temp[:,:,:,:,i-1] = fetch(@spawnat i montecarlo(Temperature,N,J_space))
    proc_complete[i] = 1
end

proc_complete[1] = 1

for i in 1:5000
    if(mean(proc_complete) == 1)
        println(proc_complete)
        @save "J1J2data"*string(N)*"x"*string(N)*"fullresbind"*string(Dates.now())*".jld2" E_temp skyrm_temp mag_temp magbind_temp skyrmbind_temp Temperature N J_space
	break
    end
    println(proc_complete)
    sleep(35)
end

println("ending")
