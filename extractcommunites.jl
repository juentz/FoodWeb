"
Takes information out of the provided excel files in the matrices folder.
Defines a list of all taxa with its (known) properties/parameters
Defines different communties that shall be analyzed with its properties

Possibility to put some extra knowledge in a single community

saves communities in extra file 'communities.jld2' which is called by 'applicationIMAPP.jl'
"

using XLSX # for extracting the matrices out of excel files
using JLD2 # for storing and loading data

include("types.jl")
include("backgroundFunctions.jl")

"GET DATA"
# Open Excel file
file = XLSX.readxlsx("matrices/21Nov25_Arctic.xlsx")
# Get worksheet
sheet = file["PredPreyMatrixArctic"]

# Extract whole foodweb + infos from excel file
# A = Float64.(coalesce.(XLSX.getdata(sheet, "D11:BT79"), 0.0)) #transforms into float and missing entries -> 0
A = Float64.(coalesce.(XLSX.getdata(sheet, "D18:CE97") .|> x -> x === nothing ? missing : x, 0.0))
A=copy(transpose(A)) #collumn eats row
N = size(A,1)
all_names = XLSX.getdata(sheet, "D7:CE7")
all_ids = vec(Int64.(coalesce.(XLSX.getdata(sheet, "D8:CE8"), 0))) # Takes ID, sets missing IDs to 0, converts into Int64
all_alpha = XLSX.getdata(sheet, "D13:CE13")
all_producer = XLSX.getdata(sheet, "D14:CE14")
# all_producer = Float64.(all_alpha .== 0)


#define list of all taxa
taxa_list = [Taxa(name=all_names[n], WoRMS=all_ids[n], producer=all_producer[n]) for n in 1:N]

#TODO if given, set other parameters  (alpha, gamma, phi etc.) for single taxa
[t.params.alpha = a for (t,a) in zip(taxa_list, all_alpha)] # set the parameter alpha for each taxa
# final metacommunity is set up
metacomm = Community(name = "Metacommunity", N= N, A=A, taxa_list=taxa_list)


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

# function run_stability(trials, s, e)
#     J = zeros(N, N)
#     stableSettings = stableJacPara[]   # vector of your struct
#     stab_trials = 0

#     for i in 1:trials
#         generalised_jacobian!(J, s, e[i])
#         println(J)
#         EVals = eigvals(J)

#         if maximum(real(EVals)) < -1e-9
#             stab_trials += 1
#             sJP = stableJacPara(copy(J), e[i])
#             push!(stableSettings, sJP)
#         end
#     end

#     return stableSettings, stab_trials
# end
#---------------------------------------------------------------------------------
"""
Define final communities for usage
"""
#---------------------------------------------------------------------------------
alphas = [t.params.alpha for t in metacomm.taxa_list]

# Print groups of species that have same entrie in the matrix and alpha values
find_redundant(metacomm, metacomm.A, alphas)
# make unique groups
unique_inds = find_unique_indices(metacomm.A, alphas)
N = length(unique_inds)
println("Keep indices: ", unique_inds)
unique_WoRMS = [metacomm.taxa_list[i].WoRMS for i in unique_inds]
unique_producer = [metacomm.taxa_list[i].producer for i in unique_inds]
# # define community without redundancies
# unique_comm = makecommunity("uniquecomm", unique_WoRMS, metacomm.A, metacomm.taxa_list) 
# # println([unique_comm.taxa_list[i].WoRMS for i in 1:N])

#define group names 
group_names = String.(getfield.(metacomm.taxa_list[unique_inds], :name) )
group_names[17] = "Spinocalanidae"
group_names[18] = "AetideidaeI"
group_names[19] = "AetideidaeII"
group_names[20] = "AetideidaeIII"
group_names[24] = "Scolecitrichidae"
group_names[25] = "Tharybidae"
group_names[30] = "Augaptilidae"
group_names[34] = "Oithonidae"
group_names[52] = "Hydrozoa"
# println(length(group_names))

group_list = [Taxa(name=group_names[n], WoRMS=unique_WoRMS[n], producer=unique_producer[n]) for n in 1:N]
group_alphas = [metacomm.taxa_list[i].params.alpha for i in unique_inds]
# println(group_alphas)

# set the parameter alpha for each group and define group community
[t.params.alpha = a for (t,a) in zip(group_list, (group_alphas))] 
group_comm = Community(name = "ArcticGroup", N = N, A = metacomm.A[unique_inds, unique_inds], taxa_list = group_list)

# println([group_comm.taxa_list[i].WoRMS for i in 1:N])
# # unique_inds = find_unique_indices(metacomm.A, alphas)
# # println("Keep indices: ", unique_inds)

# # sort matrix from small to big (high to low turnoverrate)
# degree = zeros(group_comm.N)
# for n in 1:group_comm.N
#     degree[n]=sum(group_comm.A[:,n])+sum(group_comm.A[n,:])
# end
# println(degree)
# println(getfield.(group_comm.taxa_list, :name))
A, sorted_inds = sort_matrix_by_alpha(group_comm.A, group_alphas)
# for n in 1:group_comm.N
#     degree[n]=sum(A[:,n])+sum(A[n,:])
# end
# println(degree)
# # println(sorted_inds)
sort_group_comm = Community(name = "SortedGroupCommunity", N= group_comm.N, A=A, taxa_list=group_comm.taxa_list[sorted_inds])
# println(getfield.(sort_group_comm.taxa_list, :name))

# #attach turnover (alpha) value
# alphas = alphas[sorted_inds]

# WoRMS 3 - Detritus
ex_taxa_id = 1
# ex_ind = findfirst(t -> t.WoRMS == ex_taxa_id, sort_group_comm.taxa_list) #find index of extincting taxon
noDetri_comm = ex_community(ex_taxa_id, sort_group_comm)

#without Toxic Mesoplankton
# ex_ind = findfirst(t -> t.WoRMS == 8, sort_group_comm.taxa_list) #find index of extincting taxon
noToxic_comm = ex_community(8, noDetri_comm)

# # make unique groups
# unique_inds = find_unique_indices(sortmetacomm.A, alphas)
# # println("Keep indices: ", unique_inds)
# unique_WoRMS = [sortmetacomm.taxa_list[i].WoRMS for i in unique_inds]
# unique_comm = makecommunity("uniquecomm", unique_WoRMS, sortmetacomm.A, sortmetacomm.taxa_list) # by alpha sorted community without redundancies




# Save different versions of community
communities = [metacomm, group_comm, sort_group_comm, noDetri_comm, noToxic_comm]

@save "communities_arctic.jld2" communities

