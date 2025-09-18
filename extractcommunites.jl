"
Takes information out of the provided excel files in the matrices folder
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
file = XLSX.readxlsx("matrices/Arctic.xlsx")
# Get worksheet
sheet = file["PredPreyMatrixArctic1"]

# Extract whole foodweb + infos from excel file
A = Float64.(coalesce.(XLSX.getdata(sheet, "D4:BS71"), 0.0)) #transforms into float and missing entries -> 0
N = size(A,1)
all_names = XLSX.getdata(sheet, "A4:A71")
all_ids = vec(Int64.(coalesce.(XLSX.getdata(sheet, "B4:B71"), 0))) # Takes ID, sets missing IDs to 0, converts into Int64
all_alpha = XLSX.getdata(sheet, "BU4:BU71")
all_producer = Float64.(all_alpha .== 0)


#define list of all taxa
taxa_list = [Taxa(name=all_names[n], WoRMS=all_ids[n], producer=all_producer[n]) for n in 1:N]

#TODO if given, set other parameters  (alpha, gamma, phi etc.) for single taxa
all_alpha = ifelse.(all_alpha .== 0, 1.0, 1.0 ./ (2 .* all_alpha)) # define test alpha
[t.params.alpha = a for (t,a) in zip(taxa_list, all_alpha)] # set the parameter alpha for each taxa


"""
Extracting a sub community based on WoRMS Id of single taxa - the order of Ids can be arbitrary
AND saving the communities in jld2 file
"""
comm_1 = Community(name = "Metacommunity", N= N, A=A, taxa_list=taxa_list)
comm_ids = [362051, 139178, 104465, 104852, 104292, 1507271, 104632, 110709, 1337, 137129] #test vector
comm_2 = makecommunity("test1", comm_ids, A, taxa_list)
comm_ids = [104152, 103259, 115484]
comm_3 = makecommunity("test2", comm_ids, A, taxa_list)

communities = [comm_1, comm_2, comm_3]

# Save
@save "communities.jld2" communities

