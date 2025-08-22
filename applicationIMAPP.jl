using Distributions
using UUIDs
using LinearAlgebra
using StatsBase
using Statistics

using Plots

include("types.jl")
include("backgroundFunctions.jl") 

# prey, pred small, pred big
 A = [0.0 0.0 0.0;
     1.0 1.0 0.0;
     0.0 1.0 1.0]
# # pred big, pred small, prey
# A = [1.0 1.0 0.0;
#      0.0 1.0 1.0;
#      0.0 0.0 0.0]

# amount of taxa
N = size(A, 1)
print("Amount of species = $N.\n")


# primary producer y = 1; n = 0
# producer = vec(sum(A, dims=2) .== 0) #now: 1 if primary producer
producer = [0.0, 0.0, 0.5]
producer = [0.5, 0.0, 0.0]

# iterations for random variables
iteration = 10000

"speeding up turnover rates - is irrelevant "
#turnover rates per taxa the higher the trophie the smaller
alpha = 1 ./ ([20; 7.0; 1.0])
alpha = 1 ./ ([1.0; 7.0; 20.0]) 


"DEFINE PARAMETERS"
# define (random) parameters for each of the N taxa
s = structural_parameters(N, A, alpha, producer)
e = random_parameters(N, iteration)
 
# TODO: check parameter entries of species in the web
# if parameters are given, exchange random parameters with known ones

"STABILITY RATIO"
J = zeros(N, N)
stability = zeros(iteration)
Sens = zeros(N, iteration)
Infl = zeros(N, iteration)
Sensmat = Array{Float64}(undef, N, 0)
# Inflmat = zeros(N)
#  print("EVNorms: $(size(Sens)) should be $N times $iteration")

for i = 1:iteration
    #compute Jacobean for every set of parameters
    generalised_jacobian!(J,s,e[i])
    # get left and right eigenvectors and corresponding values
    E  = eigen(J)              # E.values :: Vector{Complex}, E.vectors :: Matrix{Complex}
    

    EVals,  EVecs  = E.values,  E.vectors
    # print(size(EVals), size(EVecs))

    # Sort both spectra consistently (by real part, then imaginary part)
    p1 = sortperm(EVals,  by = x -> (real(x), imag(x)))
    
    EVals,  EVecs  = EVals[p1],  EVecs[:, p1]
    #  print(size(EVals), size(EVecs))


    # for n = 1:N
    # λ = real(eigvals!(J)[end])
    if real(EVals[end]) < 0
        stability[i] = 1# λ < 0 ? 1 : 0

        "here function, that takes J, EVals and EVecs, p1 - computes EValsT and EVecsT and returns Sensitivity and Influence vector"
        ET = eigen(transpose(J))   # transpose, to find left EVal and EVec
        EValsT, EVecsT = ET.values, ET.vectors
        # Sort both spectra consistently (by real part, then imaginary part)
        # p2 = sortperm(EValsT, by = x -> (real(x), imag(x)))
        EValsT, EVecsT = EValsT[p1], EVecsT[:, p1]
        
        # print(size(EValsT), size(EVecsT))
        # Sanity check: eigenvalues(J) ≈ eigenvalues(Jᵀ)
        if any(abs.(EVals .- EValsT) .> 1e-8)
            println("WARNING: J and transpose(J) differ beyond tolerance = 1e-8.\n ",
                    "This can be roundoff or ordering.")
            # for (a, b) in zip(EVals, EValsT)
            #     println(a, "   |   ", b)
            # end
        end

        # EVNorms = sqrt.([abs(sum(EVecsT[:, n] .* EVecs[:, n])) for n in 1:N])
        EVNorms = sqrt.(abs.(sum(EVecsT .* EVecs, dims=1)))[:]
        # print("EVNorms: $(size(EVNorms)) should be $N")
        
        for n = 1:N
            Sens[n, i] = log(sum(abs(EVecs[n,k]) * abs(real(1/EVals[k]))  for k = 1:N))
            Infl[n, i] = log(sum(abs(EVecsT[n,k]) * abs(real(1/EVals[k]))  for k = 1:N))

            # Aufderheide
            # Sens[n, i] = log.(sum(abs(EVecs[n,k]) * abs(real(1/EVals[k])) / EVNorms[k]  for k = 1:N)) # Aufderheide
            # Infl[n, i] = log.(sum(abs(EVecsT[n,k]) * abs(real(1/EVals[k])) /EVNorms[k] for k = 1:N))
        end


        # TODO function that returns vector of sens and infl for every iteration step
        # and put them together to a matrix here.
        # if i == 1
        #     Sensmat, Infmat = calculateSensandInf(N, J, EVals, EVecs, p1)
       
        # else
        # sens, inf = calculateSensandInf(N, J, EVals, EVecs, p1)
        # Sensmat, Infmat = hcat(Sensmat, sens), hcat(Infmat, inf)
        # end

        # to reduce savingspace, could as well directly compute the mean
        # inflmean .+= log.(sum(abs.(EVecsT) ./ abs.(EVals)', dims=2))
        # sensmean =+ [log(sum(abs.(EVecs[n, :]) ./ abs.(EVals))) for n in 1:size(EVecs, 1)]
        # my_sens = log.(my_sens)
        
    end
    
end
# proportion stable webs
sumstab = sum(stability)
s_rate = sumstab/iteration
print("\n $(s_rate*100) % of the conducted $iteration trials are stable.\n")


# sensitive and influential taxa
#makes the geometric mean over all the entries where the jacobi matrix was stable.
sensmean = mean(Sens[:, stability .== 1], dims=2) # geometric mean without exp
inflmean = mean(Infl[:, stability .== 1], dims=2)
# print("sensmean: $(size(sensmean)) should be $N ")
if sumstab > 0
    print("\n The mean sensitivity per taxa out of $sumstab trials:\n $sensmean\n $sensman2")
    print("\n The mean influence per taxa out of $sumstab trials:\n $inflmean.\n")
end

print("Sensibility\n max $(maximum(Sens, dims=2))\n min $(minimum(Sens, dims=2))\n")
print("Influence\n max $(maximum(Infl, dims=2))\n min $(minimum(Infl, dims=2))\n")

plot(1:iteration, Sens[2, :], 
     xlabel="Index i", 
     ylabel="Sens[1,i]", 
     title="Plot of Sens[1, i]", 
     legend=false,
     marker=:circle,    # optional: adds markers
     lw=0.5,
     ylim=(-2.5, 6))              # line width
# Save as PNG
savefig("Sens_plot.png")


"IDENTIFY IMPORTANT PARAMETERS FOR STABILITY"
corr_γ, corr_μ, corr_ϕ, corr_ψ = correlation(e, stability, N)
#Plot correlations of every parameter per taxa
# p = plot_correlations(corr_γ, corr_μ, corr_ϕ, corr_ψ)
# savefig(p, "correlations.png")

"IDENTIFY SENSITIVE AND INFLUENTIAL TAXA"
