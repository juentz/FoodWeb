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
# Load communities from 'communities.jld2' file



function run_stability(trials, s, e)
    J = zeros(N, N)
    stableSettings = stableJacPara[]   # vector of your struct
    stab_trials = 0

    for i in 1:trials
        generalised_jacobian!(J, s, e[i])
        # println(J)
        EVals = eigvals(J)

        if maximum(real(EVals)) < -1e-9
            stab_trials += 1
            print(stab_trials)
            sJP = stableJacPara(copy(J), e[i])
            push!(stableSettings, sJP)
        end
    end

    return stableSettings, stab_trials
end
# #---------------------------------------------------------------------------------
# """
# Define final community for usage
# """
# #---------------------------------------------------------------------------------
@load "communities_arctic.jld2" communities
metacomm = communities[1]
group_comm = communities[2]
sort_group_comm = communities[3]
noDetri_comm = communities[4]
noToxic_comm = communities[5]
#------------------
# Setting used
comm = noToxic_comm
trials = 10000000
filename = "stableSettings_arctic.jld2"
# -----------------

# """
# BELOW: TEST SETTING
# """
# # prey, pred small, pred big
#  A = [0.0 0.0 0.0 0.0;
#      1.0 0.0 0.0 0.0;
#      0.0 0.0 0.0 0.0;
#      1.0 1.0 1.0 0.0]

# # amount of taxa
# N = size(A, 1)
# # print("Amount of species = $N.\n")

# taxa_list = [Taxa(name="prey", WoRMS=1, producer=1.0), 
#             Taxa(name="small predator",  WoRMS=2, producer=0.5),
#             #  Taxa(name="small predator red",  WoRMS=22, producer=0.0), 
#             Taxa(name="big predator I",  WoRMS=3, producer=0.0),
#             Taxa(name="big predator II",  WoRMS=4, producer=0.0)]
# metacomm = Community(name = "test", N= N, A=A, taxa_list=taxa_list)
# comm = Community(name = "test", N= N, A=A, taxa_list=taxa_list)
# alpha =  ([0.1, 0.3, 0.3, 0.0005])
# [t.params.alpha = a for (t,a) in zip(taxa_list, alpha)] # set the parameter alpha for each taxa
# # producer = getfield.(comm.taxa_list, :producer)
# ex_taxa_id = 4 # index of taxa
# inv_taxon_id = 3
# # if Invasion == 1
# #     #make new community without inv_taxon_id
# #     comm = ex_community(inv_taxon_id, comm)
# # end
# K = [0.0, -0.5, -0.5, 0.0]
# """
# END: TEST SETTING
# """



println("$(comm.name) is set up - lets go!")

#---------------------------------------------------------------------------------
"""
Save stable parameter and jacobian for this community
"""
#---------------------------------------------------------------------------------

N = comm.N
println("Amount of species = $N.")

# define (random) parameters for each of the N taxa
A = comm.A
producer = getfield.(comm.taxa_list, :producer)
#turnover rates per taxa the higher the trophie the smaller
alpha = [t.params.alpha for t in comm.taxa_list] #takes alpha from each taxa in the community

s = structural_parameters(N, A, alpha, producer)
e = random_parameters(N, trials)

# find stable settings
stableSettings, stab_trials = run_stability(trials, s, e)
println("\n$stab_trials of $trials trials were stable.")

#save stable settings
@save filename comm trials stab_trials stableSettings s 
println("Information stored in $filename")