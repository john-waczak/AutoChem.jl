abstract type RxnJacobian end


struct BimolecularJacobianTerm <: RxnJacobian
    # du[i]/du[j]
    i::Int
    j::Int
    idxs_in::Vector{Int}  # indices of reactants
    idx_k::Int            # reaction index
    prefac::Int           # ±1 (could use bool for storage)
end

struct TrimolecularJacobianTerm <: RxnJacobian
    # du[i]/du[j]
    i::Int
    j::Int
    idxs_in::Vector{Int}  # indices of reactants
    idx_k::Int            # reaction index
    prefac::Int           # ±1 (could use bool for storage)
end

struct PhotolysisJacobianTerm <: RxnJacobian
    # du[i]/du[j]
    i::Int
    j::Int
    idxs_in::Vector{Int}  # indices of reactants
    idx_k::Int            # reaction index
    prefac::Int           # ±1 (could use bool for storage)
end




function get_jacobian_terms(drxn::BimolecularDerivativeTerm)
    jac_terms = []

    i = drxn.idx_du

    for j ∈ drxn.idxs_in
        push!(
            jac_terms,
            BimolecularJacobianTerm(
                i,
                j,
                [idx for idx ∈ drxn.idxs_in if idx != j],  # since du[i]/du[i] = 1
                drxn.idx_k,
                drxn.prefac
            )
        )
    end

    return jac_terms
end


function get_jacobian_terms(drxn::TrimolecularDerivativeTerm)
    jac_terms = []

    i = drxn.idx_du

    for j ∈ drxn.idxs_in
        push!(
            jac_terms,
            TrimolecularJacobianTerm(
                i,
                j,
                [idx for idx ∈ drxn.idxs_in if idx != j],  # since du[i]/du[i] = 1
                drxn.idx_k,
                drxn.prefac
            )
        )
    end

    return jac_terms
end


function get_jacobian_terms(drxn::PhotolysisDerivativeTerm)
    jac_terms = []

    i = drxn.idx_du

    for j ∈ drxn.idxs_in
        push!(
            jac_terms,
            PhotolysisJacobianTerm(
                i,
                j,
                [idx for idx ∈ drxn.idxs_in if idx != j],  # since du[i]/du[i] = 1
                drxn.idx_k,
                drxn.prefac
            )
        )
    end

    return jac_terms
end


function get_bimolecular_jacobian_terms(drxns)
    jac_terms = BimolecularJacobianTerm[]

    for drxn ∈ drxns
        terms = get_jacobian_terms(drxn)
        for term ∈ terms
            push!(jac_terms, term)
        end
    end

    return jac_terms
end


function get_trimolecular_jacobian_terms(drxns)
    jac_terms = TrimolecularJacobianTerm[]

    for drxn ∈ drxns
        terms = get_jacobian_terms(drxn)
        for term ∈ terms
            push!(jac_terms, term)
        end
    end

    return jac_terms
end

function get_photolysis_jacobian_terms(drxns)
    jac_terms = PhotolysisJacobianTerm[]

    for drxn ∈ drxns
        terms = get_jacobian_terms(drxn)
        for term ∈ terms
            push!(jac_terms, term)
        end
    end

    return jac_terms
end


function update_jacobian!(
    idx_t::Int,
    Jac,
    u,
    jac_term::BimolecularJacobianTerm,
    K_matrix,
    prod_temp,
    U_noint,
    n_integrated
    )

    # multiply concentrations
    for idx ∈ jac_term.idxs_in
        prod_temp *= get_concentration(idx, idx_t, u, U_noint, n_integrated)
    end

    Jac[jac_term.i, jac_term.j] += jac_term.prefac * K_matrix[jac_term.idx_k, idx_t] * prod_temp

    nothing
end


function update_jacobian!(
    idx_t::Int,
    Jac,
    u,
    jac_term::TrimolecularJacobianTerm,
    K_matrix,
    prod_temp,
    U_noint,
    n_integrated
    )

    # multiply concentrations
    for idx ∈ jac_term.idxs_in
        prod_temp *= get_concentration(idx, idx_t, u, U_noint, n_integrated)
    end

    Jac[jac_term.i, jac_term.j] += jac_term.prefac * K_matrix[jac_term.idx_k, idx_t] * prod_temp

    nothing
end


function update_jacobian!(
    idx_t::Int,
    Jac,
    u,
    jac_term::PhotolysisJacobianTerm,
    K_matrix,
    prod_temp,
    U_noint,
    n_integrated
    )

    # multiply concentrations
    for idx ∈ jac_term.idxs_in
        prod_temp *= get_concentration(idx, idx_t, u, U_noint, n_integrated)
    end

    Jac[jac_term.i, jac_term.j] += jac_term.prefac * K_matrix[jac_term.idx_k, idx_t] * prod_temp

    nothing
end




jac_func = """

function jac!(Jac, u, p, t)
    # get time value and index
    idx_t = get_time_index(t, Δt_step, ts[1])


    # set derivatives to zero
    Jac .= 0.0

    # loop over bimol terms
    prod_temp = 1.0
    @inbounds for i ∈ eachindex(jacobian_terms_bimol)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_jacobian!(
            idx_t,
            Jac,
            u,
            jacobian_terms_bimol[i],
            K_bimol,
            prod_temp,
            U_noint,
            n_integrated
        )
    end

    # loop over trimol derivatives
    prod_temp = 1.0
    @inbounds for i ∈ eachindex(jacobian_terms_trimol)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_jacobian!(
            idx_t,
            Jac,
            u,
            jacobian_terms_trimol[i],
            K_trimol,
            prod_temp,
            U_noint,
            n_integrated
        )
    end


    # loop over photolysis derivatives
    prod_temp = 1.0
    @inbounds for i ∈ eachindex(jacobian_terms_photo)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_jacobian!(
            idx_t,
            Jac,
            u,
            jacobian_terms_photo[i],
            K_photo,
            prod_temp,
            U_noint,
            n_integrated
        )
    end

    nothing
end


"""



function write_jac_func(outdir)
    outpath = joinpath(outdir, "mechanism", "jacobian.jl")

    # if it already exists, remove it so we can recreate it
    if isfile(outpath)
        rm(outpath)
    end

    open(outpath, "w") do f
        println(f, jac_func)
    end

end




function generate_jac_prototype(bimol_terms, trimol_terms, photo_terms, n_species)
    idx_pairs = []

    for term ∈ bimol_terms
        push!(idx_pairs, (term.i, term.j))
    end

    for term ∈ trimol_terms
        push!(idx_pairs, (term.i, term.j))
    end

    for term ∈ photo_terms
        push!(idx_pairs, (term.i, term.j))
    end

    # make the pairs unique
    unique!(idx_pairs)

    I = [idx_pair[1] for idx_pair ∈ idx_pairs]
    J = [idx_pair[2] for idx_pair ∈ idx_pairs]
    V = zeros(size(I))

    Jac = sparse(I,J,V,n_species, n_species)

    println("Sparsity Percentage: ", (1-(length(nonzeros(Jac))/(.*(size(Jac)...))))*100, " %")


    return Jac
end
