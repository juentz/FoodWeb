using Distributions
using UUIDs
using LinearAlgebra
using StatsBase
using Statistics, StatsPlots
using JLD2 # for storing and loading data
using Plots

#for debugging
using Debugger
using JuliaInterpreter
using Infiltrator

include("types.jl")
include("backgroundFunctions.jl") 

"""
BELOW: TEST SETTING
"""
# prey, pred small, pred big
 A = [0.0 0.0 0.0 0.0;
     1.0 0.0 0.0 0.0;
     1.0 1.0 0.0 0.0;
     1.0 1.0 1.0 0.0]

# amount of taxa
N = size(A, 1)
print("Amount of species = $N.\n")

taxa_list = [Taxa(name="prey", WoRMS=1, producer=1.0), 
            Taxa(name="small predator",  WoRMS=2, producer=0.0), 
            Taxa(name="big predator I",  WoRMS=3, producer=0.0),
            Taxa(name="big predator II",  WoRMS=4, producer=0.0)]
metacomm = Community(name = "test", N= N, A=A, taxa_list=taxa_list)
comm = Community(name = "test", N= N, A=A, taxa_list=taxa_list)
alpha =  ([1.0, 7.0, 40.0, 40.0])
[t.params.alpha = a for (t,a) in zip(taxa_list, alpha)] # set the parameter alpha for each taxa
# producer = getfield.(comm.taxa_list, :producer)
ex_taxa_id = 4 # index of taxa
inv_taxon_id = 3

#make new community without inv_taxon_id
# comm = ex_community(inv_taxon_id, comm)

K = [0.0,-1, -1.0, 0.0]
"""
END: TEST SETTING
"""
trials = 100000
res = communityAna(comm, trials, 1)

taxa_names = [comm.taxa_list[n].name for n in 1:comm.N]
plot(xlabel="trials", 
        ylabel="impact", 
        title="Impact of Perturbation $K")
    
    plot!(1:res.sum_stab, res.reactivity, label="reactivity", lw = 1)
    plot!(1:res.sum_stab, res.vulnerability, label="vulnerability", lw = 1)
    plot!(1:res.sum_stab, res.condNrs, label="condnr", lw = 1)
    # plot!(1:sum_stab, condNr, label="condNR", lw = 1, color=:black)
    savefig("diferentmeasures.png")
println(cor(res.reactivity, res.vulnerability ))
println(cor(res.reactivity, res.condNrs ))
println(cor(res.condNrs, res.vulnerability ))
println(cor(res.leading, res.reactivity))
println(cor(res.leading, res.condNrs))
exit()



"Extinction of different species----------------------------------------------"
ex_id_list = [1, 2, 3, 4]
psw, condNr, ex_psw, ex_condNrs = compareExtinctions(comm, ex_id_list)
print(psw, ex_psw, condNr, ex_condNrs)

# Left axis: psw
scatter(1:length(ex_psw), ex_psw;
    label = "ex_psw", color=:blue, marker=:circle,
    xlabel = "Index", ylabel = "psw", ylims=[0.2,1.0], legend = false,
    yforeground_color_text = :blue)

hline!( [psw]; label = "psw", color=:blue, linewidth=1, legend = false)

# Right axis: condNrs
ax2 = twinx()
scatter!(ax2, 1:length(ex_condNrs), ex_condNrs;
    label = "ex_condNrs", color=:red, marker=:diamond,
    ylabel = "condNrs", legend = false,
    yforeground_color_text = :red)

hline!(ax2, [condNr]; label = "condNr", color=:red, linewidth=1, legend = false)

# plot(p1, p2)   # overlay
savefig("expswandcond.png")

exit()

"Plots with condition number ----------------------------------------------------"
sum_stab, Sens, Infl, Impact, condNr = communityAnaPert(comm, trials, 1, 0, K)
print("\nAmount of stable webs $(sum_stab/trials).\n")
taxa_names = [comm.taxa_list[n].name for n in 1:comm.N]
plot(xlabel="trials", 
        ylabel="impact", 
        title="Impact of Perturbation $K")
    for n in 1:comm.N
        plot!(1:sum_stab, Infl[n, :], label=taxa_names[n], lw = 1)
    end
    # plot!(1:sum_stab, condNr, label="condNR", lw = 1, color=:black)
    savefig("Senstest.png")

sensmean = vec(exp.(mean(Sens, dims=2)))
inflmean = vec(exp.(mean(Infl, dims=2)))
print(inflmean)

normalize_SensInf!(comm.N, sensmean)
normalize_SensInf!(comm.N, inflmean)

taxa_names = [comm.taxa_list[n].name for n in 1:comm.N]
plot_norm_SensInf(taxa_names, sensmean, inflmean)
savefig("testEvvsSVSens.png")

# print(Sens[5,:])
# print(condNr[5])
# print(Sens[5,:]/condNr[5])
print("sumstab $(floor(sum_stab))\n")
exit()
# Plot: Impact of perturbation on taxa per trial
    plot(xlabel="trials", 
        ylabel="impact", 
        title="Impact of Perturbation $K")
    for n in 1:comm.N
        plot!(1:sum_stab, Sens[n, :], label=taxa_names[n], lw = 1)
    end
    plot!(1:sum_stab, condNr, label="condNR", lw = 1, color=:black)
    savefig("Senstorig.png")
for n in 1:Int(sum_stab)
    Impact[:,n] = Impact[:,n] * sum(condNr)/condNr[n]
end
    # Plot: Impact of perturbation on taxa per trial
    plot(xlabel="trials", 
        ylabel="impact", 
        title="Influence")
    for n in 1:comm.N
        plot!(1:sum_stab, Impact[n, :], label=taxa_names[n], lw = 1)
    end
    savefig("Sensnorm.png")

    # Plot: Impact of perturbation on taxa per trial
    plot(xlabel="trials", 
        ylabel="impact", 
        title="Impact of Perturbation $K")
    
        plot(1:sum_stab, condNr, label=false, lw = 1)
    
    savefig("condNr.png")

    # Define the bin edges
    edges = 2 .^ (0:7)   # [1,2,4,8,16,32,64,128]

    # Count how many values fall into each bin
    counts = [sum((edges[i] .<= condNr) .& (condNr .< edges[i+1])) for i in 1:length(edges)-1]
    norm_counts = counts ./ sum(counts)

    # Plot as a bar chart
    bar(0:6, norm_counts,
    xticks=(0:6, ["2^$i - 2^$(i+1)" for i in 0:6]),
    xlabel="Range", ylabel="Normalized count",
    legend=false)

    # histogram(condNr, bins=edges, normalize=:probability, legend=false)
    savefig("condNr2.png")


exit()

"check norms of collumns of neg inv jacobean ----------------------------------------"
N = comm.N
A = comm.A
# print(A)
producer = getfield.(comm.taxa_list, :producer)
#turnover rates per taxa the higher the trophie the smaller
alpha = [t.params.alpha for t in comm.taxa_list] #takes alpha from each taxa in the community
print("\nAmount of species = $N.")

# define (random) parameters for each of the N taxa
s = structural_parameters(N, A, alpha, producer)
e = random_parameters(N, trials)

J = zeros(N, N)
stability = zeros(trials)
SensEV = zeros(N, trials)
SensEVreal = zeros(N, trials)
SensSV = zeros(N, trials)
SensJ = zeros(N, trials)
SensJreal = zeros(N, trials)
Infl = zeros(N, trials)

generalised_jacobian!(J,s,e[1])


# get left and right eigenvectors and corresponding values
E  = eigen(J)              # E.values :: Vector{Complex}, E.vectors :: Matrix{Complex}
EVals,  EVecs  = E.values,  E.vectors

 # Sort both spectra consistently (by real part, then imaginary part)
p1 = sortperm(EVals,  by = x -> (real(x), imag(x)))
EVals,  EVecs  = EVals[p1],  EVecs[:, p1]

ET = eigen(J')   # transpose, to find left EVal and EVec
EValsT, EVecsT = ET.values, ET.vectors
# Sort both spectra consistently (by real part, then imaginary part)
# p2 = sortperm(EValsT, by = x -> (real(x), imag(x)))
# EValsT, EVecsT = EValsT[p1], EVecsT[:, p1]
# print(J)

cols_abs = [norm(EVecs[:,k]) for k in 1:N]
print(cols_abs)

# print(EVecs * Diagonal(EVals) * inv(EVecs))
EVecsT = inv(EVecs)
betrag2 = 0
betrag =0
cols_abs = [norm(EVecsT[k,]) for k in 1:N]
print(cols_abs)
print(EVecs* EVecsT)

e = zeros(Float64, N)
if real(EVals[end]) < -10e-6 # not 0 to avoid numerical error
    @assert eigenvalues_are_distinct(N, ComplexF64.(EVals)) == true "trial $i: Sens and Inf Ana not possible. J is not diagonalizable."
    n =1 #for in 1:comm.N #for every taxon
        SensEVreal = sum( (EVecs[n,k]) * (1/(EVals[k])) * EVecsT[k,:] for k = 1:N)
        betrag = sum( (EVecs[n,k] / (EVals[k]))^2 for k = 1:N)
        # SensEV = sum( (EVecsT[k,n]) * ((EVals[k])) * EVecs[:,k] for k = 1:N)
        SVD = svd(inv(J))
        SensSV = sum( (SVD.U[n,k]) * ((SVD.S[k])) * SVD.Vt[k,:]  for k = 1:N)
        betrag2 = sum( (SVD.U[n,k] * SVD.S[k])^2   for k = 1:N)
        e = [i == n ? 1.0 : 0.0 for i in 1:N]
        # print(e)
        # SensJ = J*e #L1 Norm
        SensJreal = e' * inv(J) 
    # end
end
# print("$SensEV\n")

print("\n$SensSV\n")

# print("$SensJ\n")
print("\n$betrag2 vs $betrag vs  $((norm(SensEVreal, 2))^2) vs $((norm(SensJreal, 2))^2)\n")

print("$SensEVreal\n")
print("$SensJreal\n")