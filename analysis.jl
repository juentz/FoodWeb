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
@load "communities.jld2" communities
metacomm = communities[1]

function find_redundant(comm::Community, A::Matrix{Float64}, alpha::Vector)
    N = size(A, 1)
    # Make a dictionary mapping signatures → indices
    groups = Dict{String, Vector{Int}}()

    for i in 1:N
        rowkey = join(A[i, :], "")
        colkey = join(A[:, i], "")
        key = string(alpha[i], "_", rowkey, "_", colkey)
        push!(get!(groups, key, Int[]), i)
    end

    # Print redundancies
    for (key, inds) in groups
        if length(inds) > 1
            # names = getfield.(comm.taxa_list[inds], :name)
            println("Redundant indices: ", inds) #names)
        end
    end
end
function sort_matrix_by_alpha(A::Matrix{Float64}, alpha::Vector)
    @assert size(A, 1) == size(A, 2) == length(alpha) "Matrix must be square and match alpha length"

    # Sort indices by alpha descending
    sorted_inds = sortperm(alpha; rev=true)

    # Permute both rows and columns
    A_sorted = A[sorted_inds, sorted_inds]

    return A_sorted, sorted_inds
end
function find_unique_indices(A::Matrix{Float64}, alpha::Vector)
    N = size(A, 1)
    @assert N == size(A, 2) "Matrix A must be square"
    @assert N == length(alpha) "Length of alpha must match A"

    # Dictionary: key = signature → indices with same key
    groups = Dict{UInt64, Vector{Int}}()

    for i in 1:N
        # Create a hash key combining alpha, row, and column pattern
        key = hash((alpha[i], A[i, :], A[:, i]))
        push!(get!(groups, key, Int[]), i)
    end

    # For each redundant group, keep only the first index
    unique_indices = [first(v) for v in values(groups)]

    return sort(unique_indices)
end


alphas = [t.params.alpha for t in metacomm.taxa_list]

# Print groups of species
find_redundant(metacomm, metacomm.A, alphas)

# make unique groups
unique_inds = find_unique_indices(metacomm.A, alphas)
N = length(unique_inds)
# println("Keep indices: ", unique_inds)
unique_WoRMS = [metacomm.taxa_list[i].WoRMS for i in unique_inds]
unique_producer = [metacomm.taxa_list[i].producer for i in unique_inds]
unique_comm = makecommunity("uniquecomm", unique_WoRMS, metacomm.A, metacomm.taxa_list) # by alpha sorted community without redundancies
# println([unique_comm.taxa_list[i].WoRMS for i in 1:N])

group_names = String.(getfield.(metacomm.taxa_list[unique_inds], :name) )
group_names[8] = "Spinocalanidae, AetideidaeI"
group_names[9] = "AetideidaeII"
group_names[13] = "Scolecitrichidae, Tharybidae"
group_names[18] = "Augaptilidae"
group_names[22] = "Oithonidae"
group_names[40] = "Hydrozoa"
# println(length(group_names))

group_list = [Taxa(name=group_names[n], WoRMS=unique_WoRMS[n], producer=unique_producer[n]) for n in 1:N]
group_alphas = [unique_comm.taxa_list[i].params.alpha for i in 1:N]
# println(group_alphas)
[t.params.alpha = a for (t,a) in zip(group_list, (group_alphas .*365))] # set the parameter alpha for each taxa
group_comm = Community(name = "ArcticGroup", N =unique_comm.N, A = unique_comm.A, taxa_list = group_list)
# println([group_comm.taxa_list[i].WoRMS for i in 1:N])
# # unique_inds = find_unique_indices(metacomm.A, alphas)
# # println("Keep indices: ", unique_inds)

# # sort matrix from small to big (high to low turnoverrate)
A, sorted_inds = sort_matrix_by_alpha(group_comm.A, group_alphas)
# # println(sorted_inds)
sort_group_comm = Community(name = "SortedGroupCommunity", N= group_comm.N, A=A, taxa_list=group_comm.taxa_list[sorted_inds])
# #attach turnover (alpha) value
# alphas = alphas[sorted_inds]

# WoRMS 3 - Detritus
ex_taxa_id = 3
ex_ind = findfirst(t -> t.WoRMS == 3, sort_group_comm.taxa_list) #find index of extincting taxon
noDetri_comm = ex_community(ex_taxa_id, sort_group_comm)

# # make unique groups
# unique_inds = find_unique_indices(sortmetacomm.A, alphas)
# # println("Keep indices: ", unique_inds)
# unique_WoRMS = [sortmetacomm.taxa_list[i].WoRMS for i in unique_inds]
# unique_comm = makecommunity("uniquecomm", unique_WoRMS, sortmetacomm.A, sortmetacomm.taxa_list) # by alpha sorted community without redundancies

comm = noDetri_comm
N = comm.N
println("$(comm.name) is set up - lets go!")
# println([group_comm.taxa_list[i].producer for i in 1:N])

# plank_WoRMS = [comm.taxa_list[i].WoRMS for i in 1:N-7] #excludes mammals, fish, birds
# plankcomm = makecommunity("PlankArctis", plank_WoRMS, comm.A[1:N-7,1:N-7], comm.taxa_list[1:N-7])
# # println(comm.A[N-10:N, N-10:N])
# # println([t.producer for t in comm.taxa_list])


"""
Stability analysis is always included.
What other of analysis shall be conducted? 
y=1, n=0
"""

"A. Analysis of the given community.
    A.1 Detection of Sensibile and Influential Taxa"
    SensandInfAna = 0
    "A.2 Impact of Perturbation"
    Perturbation = 0
    K = [0.0, 0.0, -0.5] #perturbation vector. What taxa experiences what change in growthrate F_n.

"B. Erase taxa and detect changes in stability and important taxa
    Includes analysis of original community"
Extinction = 0
ex_taxa_id = 4 #104895 # index of taxa

"C. Include taxa and detect changes in stability and important taxa
    Includes analysis of original community."
Invasion = 0
inv_taxon_id = 1337 # index of taxa

# amount of trials for random variables
trials = 1000000

"Set the communy for the analyis"
# comm = communities[1] 

# println(comm.A)

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
# if Invasion == 1
#     #make new community without inv_taxon_id
#     comm = ex_community(inv_taxon_id, comm)
# end
# K = [0.0, -0.5, -0.5, 0.0]
# """
# END: TEST SETTING
# """
# println(signed_geomean([-5.0, 0.1, 1.0]))
# exit()

"for the original community the analysis function communityAnaPert"
res = communityAna(comm, trials, 1)
println("$(res.sum_stab) out of $trials were stable.")
# println("Resilience $(res.resilience[1])")
# println("My Measure $(res.mymeasure[1])")
# println("Vulnerability $(res.vulnerability[:,1])")
# println("Reactivity $(res.reactivity[:,1])")
# println("Return Time $(res.returntime[:,1])")
# println("Resistance $(res.resistance[:,1])")
# println("$(res.sum_stab) out of $trials equilibria were stable.")
# println("Mean resilience: $(signed_geomean(res.resilience))")
# # take the geometric mean
resilience = [signed_geomean(row) for row in eachrow(res.resilience)]
println("Resilience $(resilience)")

resistance = [geomean(row) for row in eachrow(res.resistance)]
resistance_nodes = [geomean(row) for row in eachrow(res.resistance_nodes)]

returntime = [geomean(row) for row in eachrow(res.returntime)]
returntime_nodes = [geomean(row) for row in eachrow(res.returntime_nodes)]

reactivity = [signed_geomean(row) for row in eachrow(res.reactivity)]
println(reactivity)
reactive_nodes = [signed_geomean(row) for row in eachrow(res.reactive_nodes)]
println(reactive_nodes)

vulnerability = [signed_geomean(row) for row in eachrow(res.vulnerability)]
vulnerability_node = [signed_geomean(row) for row in eachrow(res.vulnerability_node)]
# println(vulnerability_node)

# # # Plot the results
# make_barplot(resistance_nodes, getfield.(comm.taxa_list, :name), "resistance", "resistancebar2.tex", resistance)
# make_barplot(returntime_nodes, getfield.(comm.taxa_list, :name), "return time", "returntimebar2.tex", returntime)
make_barplot(reactive_nodes, getfield.(comm.taxa_list, :name), "reactivity", "reactivitybar2.tex", reactivity)
make_barplot(vulnerability_node, getfield.(comm.taxa_list, :name), "absolute impact", "impactbar2.tex", vulnerability)
# plot_measure_nodes(resistance_nodes, getfield.(comm.taxa_list, :name), resistance, "Resistance Measure", "testresistance_group.tex")
# plot_measure_nodes(returntime_nodes, getfield.(comm.taxa_list, :name), returntime, "$(res.sum_stab) Return Time Measure", "testreturntime_group.png")
# plot_measure_nodes(reactive_nodes, getfield.(comm.taxa_list, :name), reactivity, "reactivity", "testreactivity_group.png")
# plot_measure_nodes(vulnerability_node, getfield.(comm.taxa_list, :name), vulnerability, "$(res.sum_stab) Fragile Equ. Measure", "testfragequil_group.png")

# Sens = [geomean(row) for row in eachrow(res.Sens)]
# Infl = [geomean(row) for row in eachrow(res.Infl)]
# make_barplot(Sens, getfield.(comm.taxa_list, :name), "sensitivity", "sensitivitybar2.tex")
# make_barplot(Infl, getfield.(comm.taxa_list, :name), "influence", "influencebar2.tex")

# normalize_SensInf!(N, Sens)
# normalize_SensInf!(N, Infl)
# plot_norm_SensInf(getfield.(comm.taxa_list, :name), Sens, Infl, "textest.tex")



if SensandInfAna == 1 || Perturbation == 1
    sum_stab, Sens, Infl, Impact = communityAnaPert(comm, trials, SensandInfAna, Perturbation, K)
end


if Extinction == 1
    println("taxa id $ex_taxa_id")
    ex_ind = findfirst(t -> t.WoRMS == ex_taxa_id, comm.taxa_list) #find index of extincting taxon
    
    # perturbation vector - how is the effect of extinction captured
    # ex_K = comm.A[ex_ind,:] .* (1/sum(comm.A[ex_ind,:])) #if taxon does not eat n, intex is 0. If it ate n we expect positive effect on growth rate of n. Anteilig wie viele prey er hatte (geht besser)
    ex_K = zeros(comm.N)
    ex_K[ex_ind] = -0.5
    sum_stab, Sens, Infl, ex_Impact = communityAnaPert(comm, trials, 1, 1, ex_K)
    println("$sum_stab of the conducted $trials trials are stable for the original community")
    # print("Sens for orig comm: $(mean(Sens, dims=2))\n")
    # make community, where ex_taxa_id is excluded
    ex_comm = ex_community(ex_taxa_id, comm)

    #last variable 1 if Sens and Inf ana shall be conducted. Else 0
    #save stability ratio, Sens and Infl matrices
    ex_res = communityAna(ex_comm, trials, 1)
    ex_sum_stab, ex_Sens, ex_Infl = ex_res.sum_stab, ex_res.Sens, ex_res.Infl
    println("Index of extincting taxon is $ex_ind")
    println("$ex_sum_stab of the conducted $trials trials are stable, if taxon $ex_taxa_id extincts.\n")
    returntime = [geomean(row) for row in eachrow(ex_res.returntime)]
    returntime_nodes = [geomean(row) for row in eachrow(ex_res.returntime_nodes)]
    println(getfield.(ex_comm.taxa_list, :name))
    plot_measure_nodes(returntime_nodes, getfield.(ex_comm.taxa_list, :name), returntime, "Return Time Measure", "testreturntime_ex.png")

    # Plot: geometric mean sens and infl per taxa -------
    # print("Sens for orig comm: $(mean(Sens, dims=2))\n")
    sensmean = vec(exp.(mean(Sens, dims=2)))
    normalize_SensInf!(comm.N,sensmean)
    # print("sensmean $(exp.(vec(mean(Sens, dims=2)))) \n")
    inflmean = vec(exp.(mean(Infl, dims=2)))
    normalize_SensInf!(comm.N,inflmean)
    # signed geomean to exclude extremes but keep negative impact negative.
    meanImp = [signed_geomean(row) for row in eachrow(ex_Impact)]
    # print(meanImp)
    # meanImp = geomean(ex_Impact, dims=2)[:]

    ex_sensmean = vec(exp.(mean(ex_Sens, dims=2)))
    normalize_SensInf!(ex_comm.N,ex_sensmean)
    insert!(ex_sensmean, ex_ind, 0)
    ex_inflmean = vec(exp.(mean(ex_Infl, dims=2)))
    normalize_SensInf!(ex_comm.N,ex_inflmean)
    insert!(ex_inflmean, ex_ind, 0)

    println(abs.(sensmean - ex_sensmean))
    println(abs.(inflmean - ex_inflmean))

    "VISUALIZATION -------------------------------------------------------------------"
    taxa_names = [comm.taxa_list[n].name for n in 1:comm.N]
    # --- Plot 1: Sensitivity & Influence original community ---
    p1 = groupedbar(
        # taxa_names,
        [sensmean inflmean],
        label = ["sensitivity" "influence"],
        xticks=(1:comm.N, taxa_names),  # taxa names on y-axis
        # ylabel = "taxa",
        ylabel = "\ngeometric mean",
        title = "original community: $(round(sum_stab/trials * 100,digits=1))% stable trials\n",
        bar_width = 0.27,
        grouped = true,
        # xmirror=true,   # mirror x-axis to the top
        # orientation = :horizontal   # make bars horizontal
    )
    p1 = plot_norm_SensInf(taxa_names, sensmean, inflmean)
    # --- Plot 2: Impact ---
    p2 = bar(
        meanImp;
        # orientation=:horizontal,
        legend=false,
        ylabel="\nimpact",
        # xlims=(-0.5, 2),
        title="perturbation $ex_K \n",
        bar_width=0.27,
        xticks=(1:comm.N, taxa_names),  # taxa names on y-axis
        # yticks=([], []),     # remove taxa labels (shared with p1)
        # xmirror=true   # mirror x-axis to the top
    )
    hline!([0], color=:black, lw=2, label=false)
    # --- Plot 1: Sensitivity & Influence original community ---
    # taxa_names = [ex_comm.taxa_list[n].name for n in 1:ex_comm.N]
    p3 = groupedbar(
        # taxa_names,
        [ex_sensmean ex_inflmean],
        label = ["sensitivity" "influence"],
        xticks=(1:comm.N, taxa_names),  # taxa names on y-axis
        # ylabel = "taxa",
        ylabel = "\ngeometric mean",
        title = "after extinction: $(round(ex_sum_stab/trials * 100,digits=1))% stable trials\n",
        bar_width = 0.27,
        grouped = true,
        # xmirror=true,   # mirror x-axis to the top
        # orientation = :horizontal   # make bars horizontal
    )
    p3 = plot_norm_SensInf(taxa_names, ex_sensmean, ex_inflmean)

    # --- Combine plots side by side ---
    plt = plot(p1, p2, p3, layout=(3,1), size=(700,700))

    savefig(plt, "Extinction.png")
end


if Invasion == 1
    # find index of the invading taxon in the metacommunity
    inv_ind = findfirst(t -> t.WoRMS == inv_taxon_id, metacomm.taxa_list)
    
    #TODO extract perturbation vector 
    old_inds = [t.WoRMS for t in comm.taxa_list] # IDs already in community
    comm_inds = [findfirst(t -> t.WoRMS == id, metacomm.taxa_list) for id in old_inds]
    # print(comm_inds)
    inv_K = metacomm.A[inv_ind, comm_inds] .* -0.5
    # print(inv_K)
    # inv_K = zeros(N)
    # if invador eats taxon set that index of K -0.2 (something negative) otherwise 0

    # analysis of origin community with perturbation
    sum_stab, Sens, Infl, Impact = communityAnaPert(comm, trials, 1, 1, inv_K)

    
    # make community, where ex_taxa_id is excluded
    inv_comm = inv_community(inv_taxon_id, comm, metacomm)
    # print("invA $(inv_comm.A)\n invproducers $(getfield.(inv_comm.taxa_list, :producer))")
    inv_sum_stab, inv_Sens, inv_Infl = communityAna(inv_comm, trials, 1)
    
    print("\n$(inv_sum_stab/trials*100) % of the conducted $trials trials are stable, if taxon $inv_taxon_id invades.\n")

    # Plot: geometric mean sens and infl per taxa -------
    inv_ind = findfirst(t -> t.WoRMS == inv_taxon_id, inv_comm.taxa_list) #find index of inversion taxon in inversion community
    sensmean = vec(exp.(mean(Sens, dims=2)))
    inflmean = vec(exp.(mean(Infl, dims=2)))
    insert!(sensmean, inv_ind, 0.0)
    insert!(inflmean, inv_ind, 0.0)
    # signed geomean to exclude extremes but keep negative impact negative.
    meanImp = [signed_geomean(row) for row in eachrow(Impact)]
    insert!(meanImp, inv_ind, 0.0)
    # print(meanImp)
    # meanImp = geomean(ex_Impact, dims=2)[:]

    inv_sensmean = vec(exp.(mean(inv_Sens, dims=2)))
    inv_inflmean = vec(exp.(mean(inv_Infl, dims=2)))

    taxa_names = [inv_comm.taxa_list[n].name for n in 1:inv_comm.N]
    "VISUALIZATION -------------------------------------------------------------------"
    # --- Plot 1: Sensitivity & Influence original community ---
    p1 = groupedbar(
        # taxa_names,
        [sensmean inflmean],
        label = ["sensitivity" "influence"],
        xticks=(1:inv_comm.N, taxa_names),  # taxa names on y-axis
        # ylabel = "taxa",
        ylabel = "\ngeometric mean",
        title = "before invasion: $(round(sum_stab/trials * 100,digits=1))% stable trials\n",
        bar_width = 0.27,
        grouped = true,
        # xmirror=true,   # mirror x-axis to the top
        # orientation = :horizontal   # make bars horizontal
    )
    # --- Plot 2: Impact ---
    p2 = bar(
        meanImp;
        # orientation=:horizontal,
        legend=false,
        ylabel="\nimpact",
        # xlims=(-0.5, 2),
        title="Perturbation $inv_K \n",
        bar_width=0.27,
        xticks=(1:inv_comm.N, taxa_names),  # taxa names on y-axis
        # yticks=([], []),     # remove taxa labels (shared with p1)
        # xmirror=true   # mirror x-axis to the top
    )
    hline!([0], color=:black, lw=2, label=false)
    # --- Plot 1: Sensitivity & Influence original community ---
    # taxa_names = [ex_comm.taxa_list[n].name for n in 1:ex_comm.N]
    p3 = groupedbar(
        # taxa_names,
        [inv_sensmean inv_inflmean],
        label = ["sensitivity" "influence"],
        xticks=(1:inv_comm.N, taxa_names),  # taxa names on y-axis
        # ylabel = "taxa",
        ylabel = "\ngeometric mean",
        title = "after invasion: $(round(inv_sum_stab/trials * 100,digits=1))% stable trials\n",
        bar_width = 0.27,
        grouped = true,
        # xmirror=true,   # mirror x-axis to the top
        # orientation = :horizontal   # make bars horizontal
    )

    # --- Combine plots side by side ---
    plt = plot(p1, p2, p3, layout=(3,1), size=(700,700))

    savefig(plt, "Invasion.png")
end




# VISUALIZATION --------------------------------------------------------------------------------
# names for each taxa
taxa_names = [comm.taxa_list[n].name for n in 1:comm.N]

# Visualization Sens and Inf analysis
if SensandInfAna == 1
    # proportion stable webs 
    print("\n$(sum_stab/trials*100) % of the conducted $trials trials of community $(comm.name) are stable.\n")
    # sensitive and influential taxa
    #makes the geometric mean over all the entries where the jacobi matrix was stable.
    #sensmean = mean(Sens[:, stability .== 1], dims=2) # geometric mean without exp
    


    # if sum_stab > 0
    #     print("\n The mean sensitivity per taxa out of $sum_stab trials:\n $sensmean\n")
    #     print("\n The mean influence per taxa out of $sum_stab trials:\n $inflmean.\n")
    # end

    # print("Sensibility\n max $(maximum(Sens, dims=2))\n min $(minimum(Sens, dims=2))\n")
    # print("Influence\n max $(maximum(Infl, dims=2))\n min $(minimum(Infl, dims=2))\n")

    # Plot: Sens and Inf for each trial -------
    # First plot (Sensitivity)
    p1 = plot(xlabel="trials", ylabel="sensitivity", title="Sensitivity of Taxa")
    for n in 1:comm.N
        plot!(p1, 1:sum_stab, Sens[n, :], label=taxa_names[n], lw=1)
    end
    # Second plot (Influence)
    p2 = plot(xlabel="trials", ylabel="influence", title="Influence of Taxa")
    for n in 1:comm.N
        plot!(p2, 1:sum_stab, Infl[n, :], label=taxa_names[n], lw=1)
    end
    # Combine them vertically and save
    plt = plot(p1, p2, layout=(2,1))
    savefig(plt, "Sens_Infl_pertrial.png")

    # Plot: geometric mean sens and infl per taxa -------
    sensmean = vec(exp.(mean(Sens, dims=2)))
    inflmean = vec(exp.(mean(Infl, dims=2)))

    plt = groupedbar(
        # taxa_names,
        [sensmean inflmean],
        label = ["Sensitivity" "Influence"],
        xticks=(1:length(taxa_names), taxa_names),  # taxa names on y-axis
        xrotation=45,
        size=(1200, 600),
        # ylabel = "Taxa",
        ylabel = "geometric mean",
        title = "Sensitivity and Influence per Taxon",
        bar_width = 0.35,
        grouped = true,
        # xmirror=true,   # mirror x-axis to the top
        # orientation = :horizontal   # make bars horizontal
    )
    # Add thick black line at x=0
    # hline!([0], color=:black, lw=2,  label=false)
    savefig(plt, "Sens_Infl_mean.png")
    # print("$(sensmean[2]) and $(inflmean[2])")
end


# Visualization of perturbation/impact analysis
if Perturbation == 1
    if SensandInfAna == 0
        # proportion stable webs 
        print("\n$(sum_stab/trials*100) % of the conducted $trials trials in communy $(comm.name) are stable.\n")
    end
    
    # Plot: Impact of perturbation on taxa per trial
    plot(xlabel="trials", 
        ylabel="impact", 
        title="Impact of Perturbation $K")
    for n in 1:comm.N
        plot!(1:sum_stab, Impact[n, :], label=taxa_names[n], lw = 1)
    end
    savefig("Impact_plot.png")

    # TODO taking log mean seams reasonable to flatten extremes
    # Plot: Mean impact per taxa
    meanImp = mean(Impact, dims=2)[:]
    bar(
        meanImp;
        # orientation=:horizontal, # horizontal bars
        legend=false,
        title="Impact of Perturbation $K", # x-axis label
        xticks=(1:length(taxa_names), taxa_names),  # taxa names on y-axis
        # xmirror=true,   # mirror x-axis to the top
        bar_width=0.27
    )
    # Add thick black line at x=0
    hline!([0], color=:black, lw=2,label=false)
    savefig("Impact_mean.png")
end
