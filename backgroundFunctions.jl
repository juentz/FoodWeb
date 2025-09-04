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








#PSW
function proportion_stable_webs!(J::Matrix{Float64}, com::Community, fexp::Function;  N_trials::Int = 100)
    #get parameter set
    s = structural_parameters(com)
    psw = 0
    #calculate jacobian
    p = generalised_parameters!(com, N_trials, s, fexp)
    for i = 1:N_trials
        generalised_jacobian!(J,p[i]) # gives a new J in dependence of p[i]
        λ = real(eigvals!(J)[end]) # is the largest real part of all eigenvalues
        psw += λ < 0 #adds 1 to psw, if the community with these random parameters had stable states
    end

    return psw / N_trials
end

function prop_stable_webs(s::StructuralParameters, N_trials::Int)
    psw = 0
    reEV = Float64.(zeros(N_trials))
    #calculate jacobian
    # p = generalised_parameters!(com, N_trials, s, fexp)
    e = random_parameters(N, iteration)
    J = zeros(s.N, s.N)
    for i = 1:N_trials
        generalised_jacobian!(J, s, e[i]) # gives a new J in dependence of p[i]
        reEV[i] = real(eigvals!(J)[end]) # is the largest real part of all eigenvalues
        psw += reEV[i] < 0 #adds 1 to psw, if the community with these random parameters had stable states
    end

    return e, reEV, psw
end

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


function stable_parameterisation(com::Community, fexp::Function, N_trials = 100)
    k = 1
    J = zeros(size(com.A))
    while k < N_trials
 
        p = generalised_parameters(com, 1, fexp)[1]
        generalised_jacobian!(J,p) 
        
        if real(eigvals!(J)[end]) < 0
            return parameterised_community(com, p)
        end
        k += 1
    end
    return com
end

# function time2stable(com)
#     t_start = time()

#     if stable_parameterisation(com, Inf) == false
#         return Inf
#     end

#     t_end = time()

#     return t_end - t_start
# end

function metacommuntiy_stability(mc)
    [communtiy_stability(c) for c = mc.coms]
end


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


"Calculates Sensitivity and Influence of each taxa"
function  calculateSensandInf(N::Int64, J::Matrix{Float64}, EVals::Vector{ComplexF64}, EVecs::Matrix{ComplexF64},p1::Vector{Int64})
    ET = eigen(transpose(J))   # transpose, to find left EVal and EVec
    EValsT, EVecsT = ET.values, ET.vectors
    # Sort both spectra consistently (by real part, then imaginary part)
    # p2 = sortperm(EValsT, by = x -> (real(x), imag(x)))
    EValsT, EVecsT = EValsT[p1], EVecsT[:, p1]
    
    # Sanity check: eigenvalues(J) ≈ eigenvalues(Jᵀ)
    if any(abs.(EVals .- EValsT) .> 1e-8)
        @warn "J and transpose(J) differ beyond tolerance = 1e-8. This can be roundoff or ordering."
    end

    # EVNorms = sqrt.([abs(sum(EVecsT[:, n] .* EVecs[:, n])) for n in 1:N])
    # EVNorms = sqrt.(abs.(sum(EVecsT .* EVecs, dims=1)))[:]
    # print("EVNorms: $(size(EVNorms)) should be $N")
    Sens = zeros(N)
    Infl = zeros(N)
    for n = 1:N
        Sens[n] = log(sum(abs(EVecs[n,k]) * abs(real(1/EVals[k]))  for k = 1:N))
        Infl[n] = log(sum(abs(EVecsT[n,k]) * abs(real(1/EVals[k]))  for k = 1:N))

        # Aufderheide
        # Sens[n, i] = log.(sum(abs(EVecs[n,k]) * abs(real(1/EVals[k])) / EVNorms[k]  for k = 1:N)) # Aufderheide
        # Infl[n, i] = log.(sum(abs(EVecsT[n,k]) * abs(real(1/EVals[k])) /EVNorms[k] for k = 1:N))
    end

    return Sens, Infl
end

"""
Takes list of WoRMS/Ids that identify the taxa in the wanted community, overall foddweb A, overall taxa_list
Returns the community
"""
function makecommunity(comm_ids::Vector{Int64}, A::Matrix{Float64}, taxa_list::Vector{Taxa})
    all_ids = getfield.(taxa_list, :WoRMS) 
    # Check if all comm_ids are contained in all_ids
    @assert all(in(all_ids), comm_ids) == true "Given Ids are not known."
    
    positions = findall(x -> x in comm_ids, all_ids) #sorted from small to big
    
    return  Community(N=length(comm_ids), A=A[positions,positions], taxa_list = taxa_list[positions])
end