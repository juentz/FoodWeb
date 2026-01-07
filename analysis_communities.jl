"""
This file analyzes different communities loaded from seperate files
"""

using Distributions
using UUIDs
using LinearAlgebra
using StatsBase
using Statistics, StatsPlots, DataFrames
using JLD2 # for storing and loading data
using Plots
pgfplotsx()
# using PGFPlotsX
using Statistics

include("types.jl")
include("backgroundFunctions.jl") 
include("templateplots.jl")

#amount of communities
L = 2
reac_resil_total = zeros(2,L)
returntime_total = zeros(3,L)
resistance_total = zeros(3,L)
psw = zeros(L)

#Load different communities
@load "stableSettings_Arctic_all.jld2" stableSettings s stab_trials trials comm
reac_resil = zeros(2, stab_trials)
returntime = zeros(3, stab_trials)
resistance = zeros(3, stab_trials)
for i in 1:stab_trials
    J = stableSettings[i].J
    
    #set up eigendecomposition of J
    E  = eigen(J)              # E.values :: Vector{Complex}, E.vectors :: Matrix{Complex}
    EVals,  EVecs  = E.values,  E.vectors
    EVecsT = inv(EVecs)

    #Resilience
    reac_resil[1, i] = maximum(real(EVals))
    # max initial change of the distance to equilibria when abundancies are shifted by a vector of norm 1
    S = (J + J') / 2       # symmetric part
    S_eigen = eigvals(S)
    reac_resil[2, i] = maximum(S_eigen)

    #return time
    obsGram = obs_Gram(EVals, EVecs, EVecsT)
    S_eigen = eigvals(obsGram)
    returntime[1, i] = minimum(S_eigen)
    returntime[2, i] = sum(S_eigen) / comm.N
    returntime[3, i] = maximum(S_eigen)

    #resistance
    contrGram = contr_Gram(EVals, EVecs, EVecsT)
    invcG = inv(contrGram)
    S_eigen = eigvals(invcG)
    resistance[1, i] = minimum(S_eigen)
    resistance[2, i] = sum(S_eigen) / comm.N
    resistance[3, i] = maximum(S_eigen)
end
psw[1] = stab_trials/trials
reac_resil_total[:,1] = [signed_geomean(row) for row in eachrow(reac_resil)]
returntime_total[:,1] = [geomean(row) for row in eachrow(returntime)]
resistance_total[:,1] = [geomean(row) for row in eachrow(resistance)]

#------------------------------------------------------------------------------------------------------------------
@load "stableSettings_arctic_noMetridinidae.jld2" stableSettings s stab_trials trials comm
reac_resil = zeros(2,stab_trials)
returntime = zeros(3, stab_trials)
resistance = zeros(3, stab_trials)
for i in 1:stab_trials
    J = stableSettings[i].J
    
    #set up eigendecomposition of J
    E  = eigen(J)              # E.values :: Vector{Complex}, E.vectors :: Matrix{Complex}
    EVals,  EVecs  = E.values,  E.vectors
    EVecsT = inv(EVecs)

    #Resilience
    reac_resil[1, i] = maximum(real(EVals))
    # max initial change of the distance to equilibria when abundancies are shifted by a vector of norm 1
    S = (J + J') / 2       # symmetric part
    S_eigen = eigvals(S)
    reac_resil[2, i] = maximum(S_eigen)

    #return time
    obsGram = obs_Gram(EVals, EVecs, EVecsT)
    S_eigen = eigvals(obsGram)
    returntime[1, i] = minimum(S_eigen)
    returntime[2, i] = sum(S_eigen) / comm.N
    returntime[3, i] = maximum(S_eigen)

    #resistance
    contrGram = contr_Gram(EVals, EVecs, EVecsT)
    invcG = inv(contrGram)
    S_eigen = eigvals(invcG)
    resistance[1, i] = minimum(S_eigen)
    resistance[2, i] = sum(S_eigen) / comm.N
    resistance[3, i] = maximum(S_eigen)
end
psw[2] = stab_trials/trials
reac_resil_total[:,2] = [signed_geomean(row) for row in eachrow(reac_resil)]
returntime_total[:,2] = [geomean(row) for row in eachrow(returntime)]
resistance_total[:,2] = [geomean(row) for row in eachrow(resistance)]
println(psw)
println(reac_resil_total)
println(returntime_total)
println(resistance_total)

x= [1.0, 2.0]

commlabels = ["Community $i" for i in 1:length(x)]

# compute ymin / ymax vectors
min = returntime_total[1, :]
max = returntime_total[3, :]
med = returntime_total[2, :]
compare_comm_plot(min, max, commlabels,"mean deviation triad", "commcompare_rtt.tex", med)

compare_comm_plot(resistance_total[1,:], resistance_total[3,:], commlabels,"resistance triad", "commcompare_resist.tex", resistance_total[2,:])
# compare_comm_plot(reac_resil_total[1,:], reac_resil_total[2,:], commlabels,"resilience-reactivity", "commcompare_reac_resil.tex")

# choose cap width in x-data units (10% of minimum x spacing)
# dx = minimum(diff(sort(x)))
# cap = 0.03    # tweak 0.08..0.2 as needed

# plt = plot(xlim=(minimum(x)-1, maximum(x)+1), ylim=(minimum(ymin), maximum(ymax)))
# for (xi, ylo, yhi) in zip(x, ymin, ymax)
#     # vertical line
#     plot!(plt, [xi, xi], [ylo, yhi], color=:black, linewidth=2, label=false)
#     # bottom cap
#     plot!(plt, [xi - cap, xi + cap], [ylo, ylo], color=:black, linewidth=2, label=false)
#     # top cap
#     plot!(plt, [xi - cap, xi + cap], [yhi, yhi], color=:black, linewidth=2, label=false)
# end

# # optional: add marker for median/mean

# scatter!(plt, x, med, color=:black, marker=:circle, markersize = 6, label=false)
# # rotated labels

# plot!(
#     plt,
#     xticks = (x, commlabels),
#     xrotation = 45,     # optional: rotate labels if long
#     yscale = :log10,       # << log10 scale
#     ylabel = "return time"
# )

# savefig(plt, "compare_arctic.png")
