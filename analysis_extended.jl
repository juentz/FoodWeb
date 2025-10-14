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
metacomm = communities[1]



"""
Stability analysis is always included.
What other of analysis shall be conducted? 
y=1, n=0
"""

"A. Analysis of the given community.
    A.1 Detection of Sensibile and Influential Taxa"
    SensandInfAna = 0
    "A.2 Impact of Perturbation"
    Perturbation = 1
    K = [0.0, 0.0, -0.5] #perturbation vector. What taxa experiences what change in growthrate F_n.

"B. Erase taxa and detect changes in stability and important taxa
    Includes analysis of original community"
Extinction = 0
ex_taxa_id = 103259 # index of taxa

"C. Include taxa and detect changes in stability and important taxa
    Includes analysis of original community."
Invasion = 0
inv_taxon_id = 1337 # index of taxa

# amount of trials for random variables
trials = 100

"Set the communy for the analyis"
comm = communities[3] 


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
if Invasion == 1
    #make new community without inv_taxon_id
    comm = ex_community(inv_taxon_id, comm)
end
K = [0.0, -0.5, -0.5, 0.0]
"""
END: TEST SETTING
"""


"for the original community the analysis function communityAnaPert"

if SensandInfAna == 1 || Perturbation == 1
    sum_stab, Sens, Infl, Impact = communityAnaPert(comm, trials, SensandInfAna, Perturbation, K)
end


if Extinction == 1
    ex_ind = findfirst(t -> t.WoRMS == ex_taxa_id, comm.taxa_list) #find index of extincting taxon
    print("Index of extincting taxon is $ex_ind\n")
    # perturbation vector - how is the effect of extinction captured
    # ex_K = comm.A[ex_ind,:] .* (1/sum(comm.A[ex_ind,:])) #if taxon does not eat n, intex is 0. If it ate n we expect positive effect on growth rate of n. Anteilig wie viele prey er hatte (geht besser)
    ex_K = zeros(comm.N)
    ex_K[ex_ind] = -0.5
    sum_stab, Sens, Infl, ex_Impact = communityAnaPert(comm, trials, 1, 1, ex_K)
    # print("Sens for orig comm: $(mean(Sens, dims=2))\n")
    # make community, where ex_taxa_id is excluded
    ex_comm = ex_community(ex_taxa_id, comm)

    #last variable 1 if Sens and Inf ana shall be conducted. Else 0
    #save stability ratio, Sens and Infl matrices
    ex_sum_stab, ex_Sens, ex_Infl = communityAna(ex_comm, trials, 1)
    
    print("\n$(ex_sum_stab/trials*100) % of the conducted $trials trials are stable, if taxon $ex_taxa_id extincts.\n")

    # Plot: geometric mean sens and infl per taxa -------
    # print("Sens for orig comm: $(mean(Sens, dims=2))\n")
    sensmean = vec(exp.(mean(Sens, dims=2)))
    # print("sensmean $(exp.(vec(mean(Sens, dims=2)))) \n")
    inflmean = vec(exp.(mean(Infl, dims=2)))
    # signed geomean to exclude extremes but keep negative impact negative.
    meanImp = [signed_geomean(row) for row in eachrow(ex_Impact)]
    # print(meanImp)
    # meanImp = geomean(ex_Impact, dims=2)[:]

    ex_sensmean = vec(exp.(mean(ex_Sens, dims=2)))
    insert!(ex_sensmean, ex_ind, 0)
    ex_inflmean = vec(exp.(mean(ex_Infl, dims=2)))
    insert!(ex_inflmean, ex_ind, 0)

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

    normalize_SensInf!(comm.N, sensmean)
    normalize_SensInf!(comm.N, inflmean)
    insert!(sensmean, inv_ind, 0.0)
    insert!(inflmean, inv_ind, 0.0)
    p1 = plot_norm_SensInf(taxa_names, sensmean, inflmean)
    plot!(p1; title="community $(comm.name)")

    # p1 = groupedbar(
    #     # taxa_names,
    #     [sensmean inflmean],
    #     label = ["sensitivity" "influence"],
    #     xticks=(1:inv_comm.N, taxa_names),  # taxa names on y-axis
    #     # ylabel = "taxa",
    #     ylabel = "\ngeometric mean",
    #     title = "before invasion: $(round(sum_stab/trials * 100,digits=1))% stable trials\n",
    #     bar_width = 0.27,
    #     grouped = true,
    #     # xmirror=true,   # mirror x-axis to the top
    #     # orientation = :horizontal   # make bars horizontal
    # )

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
    # --- Plot 1: Sensitivity & Influence invasion community ---
    # taxa_names = [ex_comm.taxa_list[n].name for n in 1:ex_comm.N]

    normalize_SensInf!(inv_comm.N, inv_sensmean)
    normalize_SensInf!(inv_comm.N, inv_inflmean)
    p3 = plot_norm_SensInf(taxa_names, inv_sensmean, inv_inflmean)
    plot!(p3; title="community $(inv_comm.name)")

    # p3 = groupedbar(
    #     # taxa_names,
    #     [inv_sensmean inv_inflmean],
    #     label = ["sensitivity" "influence"],
    #     xticks=(1:inv_comm.N, taxa_names),  # taxa names on y-axis
    #     # ylabel = "taxa",
    #     ylabel = "normalized over-/under-average",
    #     ylimits = [-1.0,1.0],
    #     title = "after invasion: $(round(inv_sum_stab/trials * 100,digits=1))% stable trials\n",
    #     bar_width = 0.27,
    #     grouped = true,
    #     # xmirror=true,   # mirror x-axis to the top
    #     # orientation = :horizontal   # make bars horizontal
    # )
    # hline!([0], color=:black, lw=2, label=false)

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

    normalize_SensInf!(N, sensmean)
    normalize_SensInf!(N, inflmean)

    #plot 
    plt = plot_norm_SensInf(taxa_names, sensmean, inflmean)
    plot!(plt; title="community $(comm.name)")
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

    #showing the reality similar to above showing all impacts per taxon
    boxplot(reshape(taxa_names, 1, :), 
        transpose(Impact), 
        legend=false)
    hline!([0], color=:black, lw=1,label=false)
    savefig("Impact_box.png")


    # signed geomean to exclude extremes but keep negative impact negative.
    meanImp = [signed_geomean(row) for row in eachrow(Impact)]
    # Plot: Mean impact per taxa
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
