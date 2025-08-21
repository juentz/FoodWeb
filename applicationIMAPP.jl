using Distributions
using UUIDs
using LinearAlgebra
using StatsBase
using Statistics
using XLSX # for extracting the matrices out of excel files
using Plots

include("types.jl")
include("generalisedjacobian.jl") 

"GET DATA"
# Open Excel file
xf = XLSX.readxlsx("matrices/Arctic.xlsx")
# Get worksheet
sheet = xf["PredPreyMatrixArctic1"]
# Extract foodweb
# A = Float64.(coalesce.(XLSX.getdata(sheet, "D4:BS71"), 0.0))

 A = [0.0 0.0 0.0;
     1.0 0.0 1.0;
     0.0 1.0 1.0]

N = size(A, 1)


# primary producer
# producer = vec(sum(A, dims=2) .== 0) #now: 1 if primary producer
producer = [0.5, 0, 0]

# iterations for random variables
iteration = 3

"speeding up turnover rates - is irrelevant "
#turnover rates per taxa
alpha = 1 ./ ([1.0; 7.0; 20.0]) 


# define (random) parameters 
s = structural_parameters(N, A, alpha, producer)
e = random_parameters(N, iteration)

"STABILITY RATIO"
J = zeros(N, N)
stability = zeros(iteration)
Sens = zeros(N,iteration)
Infl = zeros(N, iteration)

for i = 1:iteration
    #compute Jacobean for every set of parameters
    generalised_jacobian!(J,s,e[i])
    # get left and right eigenvectors and corresponding values
    E  = eigen(J)              # E.values :: Vector{Complex}, E.vectors :: Matrix{Complex}
    

    EVals,  EVecs  = E.values,  E.vectors
    

    # Sort both spectra consistently (by real part, then imaginary part)
    p1 = sortperm(EVals,  by = x -> (real(x), imag(x)))
    
    EVals,  EVecs  = EVals[p1],  EVecs[:, p1]
    


    # for n = 1:N
    # λ = real(eigvals!(J)[end])
    if real(EVals[end]) < 0
        stability[i] = 1# λ < 0 ? 1 : 0

        "here function, that takes J, EVals and EVecs - computes EValsT and EVecsT and returns Sensitivity and Influence vector"
        ET = eigen(transpose(J))   # transpose, to find left EVal and EVec
        EValsT, EVecsT = ET.values, ET.vectors
        # Sort both spectra consistently (by real part, then imaginary part)
        p2 = sortperm(EValsT, by = x -> (real(x), imag(x)))
        EValsT, EVecsT = EValsT[p2], EVecsT[:, p2]

        # Sanity check: eigenvalues(J) ≈ eigenvalues(Jᵀ)
        if any(abs.(EVals .- EValsT) .> 1e-8)
            println("WARNING: J and transpose(J) differ beyond tolerance = 1e-8.\n ",
                    "This can be roundoff or ordering.")
            # for (a, b) in zip(EVals, EValsT)
            #     println(a, "   |   ", b)
            # end
        end

        # EVNorms = [abs(sum(EVecsT[:, n] .* EVecs[:, n])) for n in 1:N]
        for n = 1:N
            Sens[n,i] = log(sum(abs(EVecs[n,k]) / abs(EVals[k]) for k in eachindex(EVals)))
        end
        # my_sens = log.(my_sens)
        print(Sens)
    end



end
s_rate = sum(stability)/iteration
print("\n $(s_rate*100) % of the conducted $iteration trials are stable.\n")


"IDENTIFY IMPORTANT PARAMETERS FOR STABILITY"
corr_γ, corr_μ, corr_ϕ, corr_ψ = correlation(e, stability, N)
#Plot correlations of every parameter per taxa
# p = plot_correlations(corr_γ, corr_μ, corr_ϕ, corr_ψ)
# savefig(p, "correlations.png")

"IDENTIFY SENSITIVE AND INFLUENTIAL TAXA"
