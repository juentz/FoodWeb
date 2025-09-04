using Distributions
using UUIDs
using LinearAlgebra
using StatsBase
using Statistics, StatsPlots
using JLD2 # for storing and loading data
using Plots

include("types.jl")
include("backgroundFunctions.jl") 
# Load communities from 'communities.jld2' file
@load "communities.jld2" communities



"""
Stability analysis is always included.
What other of analysis shall be conducted? 
y=1, n=0
"""

"A. Detection of Sensibile and Influential Taxa"
SensandInfAna = 1
"B. Impact of Perturbation"
perturbation = 1
K = [0.0, 0.0, -0.5] #perturbation vector. What taxa experiences what change.

# amount of trials for random variables
trials = 100

"Set the communy for the analyis"
comm = communities[3] 

N = comm.N
A = comm.A
producer = getfield.(comm.taxa_list, :producer)
#turnover rates per taxa the higher the trophie the smaller
alpha = [t.params.alpha for t in comm.taxa_list] #takes alpha from each taxa in the community
print("\nAmount of species = $N.\n")

"""
BELOW: TEST SETTING
"""
# prey, pred small, pred big
#  A = [0.0 0.0 0.0;
#      1.0 1.0 0.0;
#      0.0 1.0 1.0]

# # amount of taxa
# N = size(A, 1)
# print("Amount of species = $N.\n")

# # primary producer y = 1; n = 0
# # producer = vec(sum(A, dims=2) .== 0) #now: 1 if primary producer
# producer = [1, 0.5, 0.0]
# alpha = 1 ./ ([1.0, 7.0, 40.0])


"DEFINE PARAMETERS"
# define (random) parameters for each of the N taxa
s = structural_parameters(N, A, alpha, producer)
e = random_parameters(N, trials)
 

# TODO: check parameter entries of species in the web
# if parameters are given, exchange random parameters with known ones


# predefined variables for later use
J = zeros(N, N)
stability = zeros(trials)
Sens = zeros(N, trials)
Infl = zeros(N, trials)
Impact = zeros(N, trials)


for i = 1:trials
    #compute Jacobean for every set of parameters
    generalised_jacobian!(J,s,e[i])
    # get left and right eigenvectors and corresponding values
    E  = eigen(J)              # E.values :: Vector{Complex}, E.vectors :: Matrix{Complex}
    EVals,  EVecs  = E.values,  E.vectors

    # Sort both spectra consistently (by real part, then imaginary part)
    p1 = sortperm(EVals,  by = x -> (real(x), imag(x)))
    EVals,  EVecs  = EVals[p1],  EVecs[:, p1]

    # check asymptotically stable
    if real(EVals[end]) < 0
        # count stable conditions
        stability[i] = 1# λ < 0 ? 1 : 0

        # Compute value for sensitivity and influence for every taxa and trial
        if SensandInfAna == 1
            Sens[:,i], Infl[:,i] = calculateSensandInf(N, J, ComplexF64.(EVals), ComplexF64.(EVecs), p1)      
        end

        #pertubation analysis
        if perturbation == 1
                @assert size(J,1) == size(J,2) "J must be square"
                @assert size(J,1) == length(K) "J and K dimensions are incompatible"
                Impact[:, i] = -J \ K # solves Jx=K for x
        end
    end
end

# proportion stable webs 
num_stab = sum(stability)
s_rate = num_stab/trials
print("\n$(s_rate*100) % of the conducted $trials trials are stable.\n")


# VISUALIZATION --------------------------------------------------------------------------------
# names for each taxa
taxa_names = [comm.taxa_list[n].name for n in 1:N]

# Visualization Sens and Inf analysis
if SensandInfAna == 1
    # Extract only entries where stability == 1
    Sens = Sens[:, stability .== 1]
    Infl = Infl[:, stability .== 1]

    # sensitive and influential taxa
    #makes the geometric mean over all the entries where the jacobi matrix was stable.
    #sensmean = mean(Sens[:, stability .== 1], dims=2) # geometric mean without exp
    sensmean = vec(mean(Sens, dims=2))
    inflmean = vec(mean(Infl, dims=2))

    if num_stab > 0
        print("\n The mean sensitivity per taxa out of $num_stab trials:\n $sensmean\n")
        print("\n The mean influence per taxa out of $num_stab trials:\n $inflmean.\n")
    end

    # print("Sensibility\n max $(maximum(Sens, dims=2))\n min $(minimum(Sens, dims=2))\n")
    # print("Influence\n max $(maximum(Infl, dims=2))\n min $(minimum(Infl, dims=2))\n")

    # Plot: Sens and Inf for each trial -------
    # First plot (Sensitivity)
    p1 = plot(xlabel="trials", ylabel="sensitivity", title="Sensitivity of Taxa")
    for n in 1:N
        plot!(p1, 1:num_stab, Sens[n, :], label=taxa_names[n], lw=1)
    end
    # Second plot (Influence)
    p2 = plot(xlabel="trials", ylabel="influence", title="Influence of Taxa")
    for n in 1:N
        plot!(p2, 1:num_stab, Infl[n, :], label=taxa_names[n], lw=1)
    end
    # Combine them vertically and save
    plt = plot(p1, p2, layout=(2,1))
    savefig(plt, "Sens_Infl_pertrial.png")

    # Plot: Mean sens and infl per taxa -------
    plt = groupedbar(taxa_names, [sensmean inflmean], 
    label=["Sensitivity" "Influence"], 
    ylabel="mean value", 
    title="Mean Sensitivity and Influence per Taxa",
    bar_width=0.4,
    grouped=true)
    savefig(plt, "Sens_Infl_mean.png")
end


# Visualization of perturbation/impact analysis
if perturbation == 1
    # Extract only entries where stability == 1
    Impact = Impact[:, stability .== 1]

    # Plot: Impact of perturbation on taxa per trial
    plot(xlabel="trials", 
        ylabel="impact", 
        title="Impact of Perturbation $K")
    for n in 1:N
        plot!(1:num_stab, Impact[n, :], label=taxa_names[n], lw = 1)
    end
    savefig("Impact_plot.png")

    # Plot: Mean impact per taxa
    meanImp = mean(Impact, dims =2)
    plot( bar(taxa_names, meanImp))
    savefig("Impact_mean.png")
end
