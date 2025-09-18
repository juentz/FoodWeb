#TODO check parameter alpha
# N = amount of different species
# A= N times N Food web matrix consisting of entries representing who eats who
# BM = body mass
# M = amount of tryies for randaom parameters
# producer = Vector that tells us if species n is a primary producer - 0 if primary producer

# function structural_parameters(N::Int64, A::Matrix{Float64}, M::Vector{Float64}, n::Vector{Float64}, R::Float64, producer::Vector{Bool})
function structural_parameters(N::Int64, A::Matrix{Float64}, BM::Vector{Float64}, producer::Vector{Float64})

    #proportion of predation of j done by i - zero if i does not eat j
    # Julie: denotes the
    # relative weight of the contribution of species m (row) to preda-
    # tive loss rate of species n (collum)  - /beta[m,n] =0 if m does not eat n
    β = A ./ (norm.(eachcol(A), 1)') # ' = transpose

    # proportion of diet of i made up by j; zero if j does not eat i
    # Julie: measures the relative con-
    # tribution of population n to the total amount of food that
    # is available to species m - /xi[m,n] =0 if m does not eat n
    χ = A ./ (norm.(eachrow(A), 1))

    # sets all NaN elements (because of devision by 0) to 0.0
    β[isnan.(β)] .= 0.0
    χ[isnan.(χ)] .= 0.0

      
    #Julie: rate of biomass flow in a steady state - Vector; Why this exponent?; could devide by alpha[1] then other entries are the ration of live expectancies 
    # is expected to increase in warmer water
    α = (BM .^ -0.25)

    #growth from production; one if species is primary producer
    ρ̃ = producer #1 .- ρ

    # growth from predation 
    # Julie: should be fraction not total. describes what fraction of the growth rate
    # (for every species) is gained from predation; entry of ro is zero, if species is primary producer
    # here: only 0 or 1
    ρ = 1 .- ρ̃ 
    #ρ = 1 .* .!producer
    
    
    
    "maybe better if fraction if ende der nahrungskette sigma = 0, sonst etwas zwischen 0.5 und 1"
    # relative loss from predation listed in the foodweb A
    # here: either 0 or 1
    σ = 1 .* (sum(A,dims=1) .!= 0)[:]
    # σ = rand(Uniform(0.5, 1.0), 1, iteration) .* (sum(A,dims=1) .!= 0)[:]

    #relative loss from (other) mortality
    σ̃ = 1 .- σ

    return StructuralParameters(N,A,α,β,χ,ρ,ρ̃,σ,σ̃)
end

# structural_parameters(com::Community) = structural_parameters(com.N, com.A, [x.M[1] for x = com.sp], com.n, com.R, [x.producer[1] for x = com.sp])
# structural_parameters(com::Community) = structural_parameters(com.N, com.A, [x.M[1] for x = com.sp], [x.producer[1] for x = com.sp])

# N = amount of species
# M = amount of trials
function random_parameters(N::Int64, M::Int64)
    #exponent
    γ = rand(Uniform(0.5, 1.5), N, M) #[0.8, 1.5]  nonlinearity of the predation rate on n with respect to prey density; prey abundant then gamma = 0 why not up to 2?
    λ = ones(N,N) # 1 lambda[m,n] nonlinearity of the contribution of population m to the diet of population n - 1 if predator does not distuinguish between prey - 2 if predator is adapted to prey - 0 if n does not eat m
    μ = rand(Uniform(1.0, 2.0), N, M) #[1.0, 2.0] exponent of mortality (linear 1.0, quadratic 2.0) 
    ϕ = rand(Uniform(0.0, 1.0), N, M) #[0.0, 1.0] nutrient availablility - 0 if limited nutrients - 1 if abundant nutrients and no other limiting factors - phi is in the range between 0 and 1
    ψ = rand(Uniform(0.5,1.5), N, M) #[0.5, 1.2] why not [0,1] psi=1 if desity of species n has no effect on predation rate of n itself (linear) - Holling Type

    return [ExponentialParameters(γ[:,i], λ, μ[:,i], ϕ[:,i], ψ[:,i]) for i = 1:M]
end

# random_parameters(com::Community, M::Int64) = random_parameters(com.N, M)

"""
    generalised_parameters(com, K)

Randomly samples a set of generalised parameters from a given interaction matrix and set of niches 
"""
function generalised_parameters(com::Community, K::Int64, fexp::Function)
    sp = structural_parameters(com)
    ep = fexp(com.N , K)

    return GeneralisedParameters.(Ref(sp), ep)
end

"""
    generalised_parameters!(com::Community, K::Int64, sp::StructuralParameters)

TBW
"""
function generalised_parameters!(com::Community, K::Int64, sp::StructuralParameters, fexp::Function)
    ep = fexp(com.N , K)
    return GeneralisedParameters.(Ref(sp), ep)
end

"""
    generalised_jacobian!(J,p::GeneralisedParameters)

Constructs the jaccobian in place as defined by the generalised model. Takes parameter values (including network structure) from the `GeneralisedParameters` object `p`.
"""
function generalised_jacobian!(J,s::StructuralParameters, e::ExponentialParameters)
    #construct generalised_jacobian
    for i = 1:s.N
        for j = 1:s.N
            if i == j #intra Julie: computes diagonal entries of the jacobian 
                J[i,i] = s.ρ̃[i] * e.ϕ[i] + #primary production
                         s.ρ[i] * (e.γ[i] * s.χ[i,i] * e.λ[i,i] + e.ψ[i]) - #production by predation
                         s.σ̃[i] * e.μ[i] - #mortality
                         s.σ[i] * s.β[i,i] * e.ψ[i] #predation loss
                for k = 1:s.N
                    J[i,i] -= s.σ[i] * s.β[k,i] * e.λ[k,i] * ((e.γ[k] - 1) * s.χ[k,i] + 1) #predation loss
                end
            else #inter
                J[i,j] = 0

                if s.χ[i,j] != 0 #means i eats j
                    J[i,j] = s.ρ[i] * e.γ[i] * s.χ[i,j] * e.λ[i,j] 
                end

                if s.β[j,i] != 0
                    J[i,j] -= s.σ[i] * s.β[j,i] * e.ψ[j]
                end

                for k = 1:s.N
                    if (s.β[k,i] != 0) && (s.χ[k,j] != 0)
                        J[i,j] -= s.σ[i] * (s.β[k,i] * e.λ[k,j] * (e.γ[k] - 1) * s.χ[k,j])
                    end
                end
            end

            J[i,j] *= s.α[i]
        end
    end
end

# generalised_jacobian!(J,p::GeneralisedParameters) = generalised_jacobian!(J, p.s, p.e)


"""
    generalised_jacobian(p::GeneralisedParameters)

generate jacobian from parameter set p
"""
function generalised_jacobian(s::StructuralParameters, e::ExponentialParameters)
    J = zeros(s.N, s.N)
    generalised_jacobian!(J,s,e) # mutates J
    return J 
end

# generalised_jacobian(p::GeneralisedParameters) = generalised_jacobian(p.s, p.e)








# #PSW
# function proportion_stable_webs!(J::Matrix{Float64}, com::Community, fexp::Function;  N_trials::Int = 100)
#     #get parameter set
#     s = structural_parameters(com)
#     psw = 0
#     #calculate jacobian
#     p = generalised_parameters!(com, N_trials, s, fexp)
#     for i = 1:N_trials
#         generalised_jacobian!(J,p[i]) # gives a new J in dependence of p[i]
#         λ = real(eigvals!(J)[end]) # is the largest real part of all eigenvalues
#         psw += λ < 0 #adds 1 to psw, if the community with these random parameters had stable states
#     end

#     return psw / N_trials
# end

# function prop_stable_webs(s::StructuralParameters, N_trials::Int)
#     psw = 0
#     reEV = Float64.(zeros(N_trials))
#     #calculate jacobian
#     # p = generalised_parameters!(com, N_trials, s, fexp)
#     e = random_parameters(N, iteration)
#     J = zeros(s.N, s.N)
#     for i = 1:N_trials
#         generalised_jacobian!(J, s, e[i]) # gives a new J in dependence of p[i]
#         reEV[i] = real(eigvals!(J)[end]) # is the largest real part of all eigenvalues
#         psw += reEV[i] < 0 #adds 1 to psw, if the community with these random parameters had stable states
#     end

#     return e, reEV, psw
# end

# """
#     proportion_stable_webs(com::Community, fexp::Function, N_trials::Int = 100)

# Calculates the proportion stable parameter configurations for a Community.   
# """
# function proportion_stable_webs(com::Community, fexp::Function; N_trials::Int = 100)
#     J = zeros(com.N, com.N)
#     return proportion_stable_webs!(J,com, fexp; N_trials = N_trials)
# end
 
# """
#     proportion_stable_webs(mc::MetaCommunity; N_trials::Int = 100)

# Calculates proprtion stable webs across a metacommuntiy, returning an array of proportions. 
# """
# function proportion_stable_webs(mc::MetaCommunity, fexp::Function; N_trials::Int = 100)
#     props = zeros(size(mc.coms))
#     for (i,c) = enumerate(mc.coms)
#         if length(c.sp) > 0
#             J = zeros(size(c.A))
#             props[i] = proportion_stable_webs!(J,c, fexp, N_trials = N_trials)
#         end
#     end
#     return props 
# end


# function stable_parameterisation(com::Community, fexp::Function, N_trials = 100)
#     k = 1
#     J = zeros(size(com.A))
#     while k < N_trials
 
#         p = generalised_parameters(com, 1, fexp)[1]
#         generalised_jacobian!(J,p) 
        
#         if real(eigvals!(J)[end]) < 0
#             return parameterised_community(com, p)
#         end
#         k += 1
#     end
#     return com
# end

# function time2stable(com)
#     t_start = time()

#     if stable_parameterisation(com, Inf) == false
#         return Inf
#     end

#     t_end = time()

#     return t_end - t_start
# end

# function metacommuntiy_stability(mc)
#     [communtiy_stability(c) for c = mc.coms]
# end


"Correlation between single parameters and stability
 returns vector with correlations of that parameter per taxa"
function correlation(e::Vector{ExponentialParameters}, s::Vector{Float64}, N::Int64)
    corr_γ = zeros(N)
    corr_μ = zeros(N)
    corr_ϕ = zeros(N)
    corr_ψ = zeros(N)

    Γ = hcat([ex.γ for ex in e]...)
    Μ = hcat([ex.μ for ex in e]...)
    Φ = hcat([ex.ϕ for ex in e]...)
    Ψ = hcat([ex.ψ for ex in e]...)

    for i = 1:N
        corr_γ[i] = cor(Γ[i, :], s)
        corr_μ[i] = cor(Μ[i, :], s)
        corr_ϕ[i] = cor(Φ[i, :], s)
        corr_ψ[i] = cor(Ψ[i, :], s)
    end

    return corr_γ, corr_μ, corr_ϕ, corr_ψ
end


""
function plot_correlations(corr_γ, corr_μ, corr_ϕ, corr_ψ)
    N = length(corr_γ)
    taxa = 1:N
    
    p1 = bar(taxa, corr_γ, title="γ", xlabel="taxa", ylabel="correlation", label = false)
    p2 = bar(taxa, corr_μ, title="μ", xlabel="taxa", ylabel="correlation", label = false)
    p3 = bar(taxa, corr_ϕ, title="ϕ", xlabel="taxa", ylabel="correlation", label = false)
    p4 = bar(taxa, corr_ψ, title="ψ", xlabel="taxa", ylabel="correlation", label = false)

    plot(p1, p2, p3, p4; layout=(4,1), size=(600,800))
end

function biorthonormalize!(EVecsT::Matrix{ComplexF64}, EVecs::Matrix{ComplexF64}, EVals::Vector{ComplexF64})
    n = size(EVecs, 2)
    for j in 1:n
        # dot product of left and right eigenvector j
        norm_factor = sqrt(dot(EVecsT[:, j], EVecs[:, j]))
        
        # normalize both so that <w_j, v_j> = 1
        EVecsT[:, j] ./= conj(norm_factor)
        EVecs[:, j] ./= norm_factor
        EVals[j] = EVals[j]*norm_factor*norm_factor
    end
    return EVecsT, EVecs, EVals
end


"Calculates Sensitivity and Influence of each taxa"
function  calculateSensandInf(N::Int64, J::Matrix{Float64}, EVals::Vector{ComplexF64}, EVecs::Matrix{ComplexF64},p1::Vector{Int64})
    @assert eigenvalues_are_distinct(N, ComplexF64.(EVals)) == true "trial $i: Sens and Inf Ana not possible. J is not diagonalizable."

    ET = eigen(J')   # transpose, to find left EVal and EVec
    EValsT, EVecsT = ET.values, ET.vectors
    # Sort both spectra consistently (by real part, then imaginary part)
    # p2 = sortperm(EValsT, by = x -> (real(x), imag(x)))
    EValsT, EVecsT = EValsT[p1], EVecsT[:, p1]
    
    # Sanity check: eigenvalues(J) ≈ eigenvalues(Jᵀ)
    if any(abs.(EVals .- EValsT) .> 1e-8)
        @warn "J and transpose(J) differ beyond tolerance = 1e-8. This can be roundoff or ordering."
    end
    # print("$EVals \n")
    # # #Check right norms 
    # biorthonormalize!(EVecsT, EVecs) #ensures the right normalization

    #for the normalization EVecsT i times EVecs j = Kronecker delta
    # EVecsT, EVecs, EVals = biorthonormalize!(EVecsT, EVecs, EVals)
    #check if already given, if about one it is fine.
    # dist = maximum(abs.(J .- EVecs * Diagonal(EVals) * EVecsT'))
    # print("DISTANCE $dist\n")
    #check if eigenvectors are normalized
    # cols_abs = [norm(EVecs[:,k]) for k in 1:N]
    # println(cols_abs)



    # EVNorms = sqrt.([abs(sum(EVecsT[:, n] .* EVecs[:, n])) for n in 1:N]) 
    # print(abs.(EVecsT' * EVecs)) #test for biorthonormalize
    Sens = zeros(N)
    Infl = zeros(N)
    
    for n = 1:N
        Sens[n] = log(sum( abs(EVecs[n,k]) * abs((1/EVals[k]))   for k = 1:N))
        Infl[n] = log(sum( abs(EVecsT[n,k]) * abs((1/EVals[k]))   for k = 1:N))
        # Sens[n] = log(sum( abs(EVecs[n,k]) * abs((1/EVals[k]))   for k = 1:N))
        # Infl[n] = log(sum( abs(EVecsT[n,k]) * abs((1/EVals[k]))   for k = 1:N))
    
        #Aufderheide
        # Sens[n] = log(sum( abs(EVecs[n,k]) * abs(real(1/EVals[k])) / EVNorms[k]  for k = 1:N))
        # Infl[n] = log(sum( abs(EVecsT[n,k]) * abs(real(1/EVals[k])) / EVNorms[k]  for k = 1:N))
    end
    
    return Sens, Infl
end

"""
Takes list of WoRMS/Ids that identify the taxa in the wanted community, overall foddweb A, overall taxa_list
Returns the community
"""
function makecommunity(name::String, comm_ids::Vector{Int64}, A::Matrix{Float64}, taxa_list::Vector{Taxa})
    all_ids = getfield.(taxa_list, :WoRMS) 
    # Check if all comm_ids are contained in all_ids
    @assert all(in(all_ids), comm_ids) == true "Given Ids are not known."
    
    positions = findall(x -> x in comm_ids, all_ids) #sorted from small to big
    
    return  Community(name =name, N=length(comm_ids), A=A[positions,positions], taxa_list = taxa_list[positions])
end

function ex_community(ex_taxa_id::Int64, comm::Community)
    #check extinct species is in community
    @assert any(t -> t.WoRMS == ex_taxa_id, comm.taxa_list) "WoRMS ID $ex_taxa_id not found in the community"
    #find index of the extinct taxon
    ex_ind = findfirst(t -> t.WoRMS == ex_taxa_id, comm.taxa_list) 
    # update community after extinction
    ex_name = "$(comm.name), $(comm.taxa_list[ex_ind].name) extinct"
    ex_N = comm.N -1
    ex_A = comm.A[setdiff(1:comm.N, ex_ind), setdiff(1:comm.N, ex_ind)]
    ex_taxa_list = comm.taxa_list[setdiff(1:comm.N, ex_ind)]
    return Community(name = ex_name, N= ex_N, A=ex_A, taxa_list=ex_taxa_list)
end

function inv_community(inv_taxon_id::Int64, comm::Community, metacomm::Community)
    # check that invader exists in the metacommunity
    @assert any(t -> t.WoRMS == inv_taxon_id, metacomm.taxa_list) "WoRMS ID $inv_taxon_id not found in metacommunity"

    # check that invader is NOT already in the subcommunity
    @assert !any(t -> t.WoRMS == inv_taxon_id, comm.taxa_list) "WoRMS ID $inv_taxon_id already exists in community"

    # find index of the invading taxon in the metacommunity
    inv_ind = findfirst(t -> t.WoRMS == inv_taxon_id, metacomm.taxa_list)

    # update community after invasion
    inv_name = "$(comm.name), $(metacomm.taxa_list[inv_ind].name) invading"
    inv_N = comm.N + 1

    # expand adjacency matrix: keep existing interactions + add row/col for invader
    old_inds = [t.WoRMS for t in comm.taxa_list]                # IDs already in community
    comm_inds = [findfirst(t -> t.WoRMS == id, metacomm.taxa_list) for id in old_inds]
    new_inds = vcat(comm_inds, inv_ind)                    # indices of new subcommunity

    inv_A = metacomm.A[new_inds, new_inds]
    inv_taxa = metacomm.taxa_list[new_inds]

    return Community(name = inv_name, N = inv_N, A = inv_A, taxa_list=inv_taxa)
end


function eigenvalues_are_distinct(N::Int64, EVals::Vector{ComplexF64})
    for i in 1:N-1
        for j in i+1:N
            if abs(EVals[i] - EVals[j]) < 1e-18
                return false
            end
        end
    end
    return true
end

function communityAna(comm::Community, trials::Int64, SensandInfAna::Int64)
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

    #TODO check for given parameters and excange

    # predefined variables for later use
    J = zeros(N, N)
    stability = zeros(trials)
    Sens = zeros(N, trials)
    Infl = zeros(N, trials)

    for i = 1:trials
        #compute Jacobean for every set of parameters
        generalised_jacobian!(J,s,e[i])
        # get left and right eigenvectors and corresponding values
        E  = eigen(J)              # E.values :: Vector{Complex}, E.vectors :: Matrix{Complex}
        EVals,  EVecs  = E.values,  E.vectors

        "test`````````````````````````````````"
        # ET = eigen(J')   # transpose, to find left EVal and EVec
        # EValsT, EVecsT = ET.values, ET.vectors
        # #Check right norms 
        # print(" shall be 1: $(dot(EVecsT[:, 1], EVecs[:, 1]))\n")
        # print(" shall be 0: $(dot(EVecsT[:, 2], EVecs[:, 1]))\n")
        "test`````````````````````````````````"

        # Sort both spectra consistently (by real part, then imaginary part)
        p1 = sortperm(EVals,  by = x -> (real(x), imag(x)))
        EVals,  EVecs  = EVals[p1],  EVecs[:, p1]

        # check asymptotically stable
        if real(EVals[end]) < -10e-6 # not 0 to avoid numerical error
            # count stable conditions
            stability[i] = 1# λ < 0 ? 1 : 0

            # Compute value for sensitivity and influence for every taxa and trial
            if SensandInfAna == 1
                Sens[:,i], Infl[:,i] = calculateSensandInf(N, J, ComplexF64.(EVals), ComplexF64.(EVecs), p1)    
            end
        end
    end
    if SensandInfAna == 1
        # Extract only entries where stability == 1
        Sens = Sens[:, stability .== 1]
        Infl = Infl[:, stability .== 1]
    end

    return sum(stability), Sens, Infl
end

function communityAnaPert(comm::Community, trials::Int64, SensandInfAna::Int64, Perturbation::Int64, K::Vector{Float64})
    N = comm.N
    A = comm.A
    # print(A)
    producer = getfield.(comm.taxa_list, :producer)
    #turnover rates per taxa the higher the trophie the smaller
    alpha = [t.params.alpha for t in comm.taxa_list] #takes alpha from each taxa in the community
    print("\nAmount of species = $N.\n")

    # define (random) parameters for each of the N taxa
    s = structural_parameters(N, A, alpha, producer)
    e = random_parameters(N, trials)

    #TODO check for given parameters and excange

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
        if real(EVals[end]) < -10e-6 # not 0 to avoid numerical error
            # count stable conditions
            stability[i] = 1# λ < 0 ? 1 : 0

            # Compute value for sensitivity and influence for every taxa and trial
            if SensandInfAna == 1
                Sens[:,i], Infl[:,i] = calculateSensandInf(N, J, ComplexF64.(EVals), ComplexF64.(EVecs), p1)
            end

            #pertubation analysis
            if Perturbation == 1
                    @assert size(J,1) == size(J,2) "J must be square"
                    @assert size(J,1) == length(K) "J and K dimensions are incompatible"
                    Impact[:, i] = -J \ K # solves Jx=K for x
            end
        end
    end
    if SensandInfAna == 1
        # Extract only entries where stability == 1
        Sens = Sens[:, stability .== 1]
        Infl = Infl[:, stability .== 1]
    end
    if Perturbation == 1
         Impact = Impact[:, stability .== 1]
    end

    return sum(stability), Sens, Infl, Impact
end

"Makes the mean for the Impact which can be negative as well."
function signed_geomean(x::AbstractVector)
    if isempty(x)
        return NaN
    end

    # 1. Identify positive and negative indices
    pos_idx = findall(>(0), x)
    neg_idx = findall(<(0), x)

    # 2. Take logs
    log_x = similar(x, Float64)
    log_x[pos_idx] .= log.(x[pos_idx])
    log_x[neg_idx] .= log.(abs.(x[neg_idx]))

    # 3. Restore signs for negative entries
    log_x[neg_idx] .= log_x[neg_idx] .* -1

    # 4. Mean over all entries
    mean_log = mean(log_x)

    # 5. Exponentiate - dont to that want to keep negative impacts negative.
    return mean_log
end