"""
This file analyzes a given matrix which is loded from jld2 file regarding 
community measurements as well as node measurements
"""

using Distributions
using UUIDs
using LinearAlgebra
using StatsBase
using Statistics, StatsPlots
using JLD2 # for storing and loading data
using Plots
pgfplotsx()
# using PGFPlotsX
using Statistics

include("types.jl")
include("backgroundFunctions.jl") 
include("templateplots.jl")

@load "stableSettings_arctic.jld2" stableSettings s stab_trials trials comm


N = comm.N
println("Compareing $N species.")

mymeasure = zeros(stab_trials)
leading = zeros(stab_trials)

resilience = zeros(stab_trials)
svmaxV = zeros(stab_trials)
svminV = zeros(stab_trials)

returntime = zeros(3, stab_trials)
returntime_nodes = zeros(N, stab_trials)

resistance = zeros(3, stab_trials)
resistance_nodes = zeros(N, stab_trials)

degree = zeros(N)
for n in 1:N
    degree[n]=sum(comm.A[:,n])+sum(comm.A[n,:])
end
println(degree)


Sens = zeros(N, stab_trials)
Infl = zeros(N, stab_trials)
vulnerability = zeros(3, stab_trials)
vulnerability_node = zeros(N, stab_trials)

reactivity = zeros(3, stab_trials)
reactive_nodes = zeros(N, stab_trials)


for i in 1:stab_trials
    J = stableSettings[i].J
    
    #set up eigendecomposition of J
    E  = eigen(J)              # E.values :: Vector{Complex}, E.vectors :: Matrix{Complex}
    EVals,  EVecs  = E.values,  E.vectors
    EVecsT = inv(EVecs)

    #Resilience
    resilience[i] = maximum(real(EVals))
    vals = svdvals(EVecs)
    svmaxV[i] = maximum(vals)
    svminV[i] = minimum(vals)

    #return time
    obsGram = obs_Gram(EVals, EVecs, EVecsT)
    S_eigen = eigvals(obsGram)
    returntime[1, i] = minimum(S_eigen)
    returntime[2, i] = sum(S_eigen) / comm.N
    returntime[3, i] = maximum(S_eigen)
    returntime_nodes[:, i] = diag(obsGram)

    #resistance
    contrGram = contr_Gram(EVals, EVecs, EVecsT)
    invcG = inv(contrGram)
    S_eigen = eigvals(invcG)
    resistance[1, i] = minimum(S_eigen)
    resistance[2, i] = sum(S_eigen) / comm.N
    resistance[3, i] = maximum(S_eigen)
    resistance_nodes[:, i] = diag(invcG)

    # #Hankel singular values
    # W = eigvals(contrGram * obsGram)
    # hankelsv[:,i] = sqrt.(W)

    Sens[:,i], Infl[:,i] = calculateSensandInf(N, J, ComplexF64.(EVals), ComplexF64.(EVecs), ComplexF64.(EVecsT))
    D = svdvals(J)
    vulnerability[1,i] = 1/(D[1]^2) 
    vulnerability[2,i] =  (1 / comm.N) * sum(1 ./ (D .^2))
    vulnerability[3,i] = 1/(D[comm.N]^2) 
    J_inv = inv(J)
    # println(J_inv)
    for n in 1:comm.N
        vulnerability_node[n,i] = (comm.taxa_list[n].params.alpha * norm(J_inv[:,n]))^2 #shall take column
    end

    # max initial change of the distance to equilibria when abundancies are shifted by a vector of norm 1
    S = (J + J') / 2       # symmetric part
    S_eigen = eigvals(S)
    reactivity[1, i] = minimum(S_eigen)
    reactivity[2, i] = sum(S_eigen) / comm.N
    reactivity[3, i] = maximum(S_eigen)
    reactive_nodes[:, i] = real.(diag(J))
end
println(resilience)
resilience = signed_geomean(resilience)
println("Resilience $(resilience)")
svmaxV = signed_geomean(svmaxV)
svminV = signed_geomean(svminV)
println("Singular values max: $(svmaxV), min: $svminV, condition $(svmaxV/svminV)")

resistance = [geomean(row) for row in eachrow(resistance)]
resistance_nodes = [geomean(row) for row in eachrow(resistance_nodes)]

returntime = [geomean(row) for row in eachrow(returntime)]
returntime_nodes = [geomean(row) for row in eachrow(returntime_nodes)]

# hankelsv = [geomean(row) for row in eachrow(hankelsv)]

Sens = [geomean(row) for row in eachrow(Sens)]
Infl = [geomean(row) for row in eachrow(Infl)]
vulnerability = [signed_geomean(row) for row in eachrow(vulnerability)]
vulnerability_node = [signed_geomean(row) for row in eachrow(vulnerability_node)]

reactivity = [signed_geomean(row) for row in eachrow(reactivity)]
# println(reactivity)
# reactive_nodes = [signed_geomean(row) for row in eachrow(reactive_nodes)]
# println(reactive_nodes)

# # Plot the results
make_barplot(resistance_nodes, getfield.(comm.taxa_list, :name), "resistance", "resistance_arctic.tex", resistance)
make_barplot(returntime_nodes, getfield.(comm.taxa_list, :name), "return time", "returntime_arctic.tex", returntime)
# make_barplot(reactive_nodes, getfield.(comm.taxa_list, :name), "initial amplification rate", "reactive_arctic.png", reactivity)
make_barplot(Sens, getfield.(comm.taxa_list, :name), "sensitivity", "sensitivity_arctic.tex")
make_barplot(Infl, getfield.(comm.taxa_list, :name), "influence", "influence_arctic.tex")
make_barplot(vulnerability_node, getfield.(comm.taxa_list, :name), "absolute impact (normed pert)", "impact_arctic_norm.tex")

# # make_barplot(reactive_nodes, getfield.(comm.taxa_list, :name), "reactivity", "reactivity_arctic.png", reactivity)
# plt = scatter(1:N, degree)
# savefig(plt,"degree_arctic.png")
# println(getfield.(comm.taxa_list, :name))