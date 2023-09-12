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
        # prod_temp *= u[idx]
        prod_temp *= get_concentration(idx, idx_t, u, U_noint, n_integrated)
    end

    Jac[jac_term.i, jac_term.j] += jac_term.prefac * K_matrix[jac_term.idx_k, idx_t] * prod_temp
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
        # prod_temp *= u[idx]
        prod_temp *= get_concentration(idx, idx_t, u, U_noint, n_integrated)
    end

    Jac[jac_term.i, jac_term.j] += jac_term.prefac * K_matrix[jac_term.idx_k, idx_t] * prod_temp
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
        # prod_temp *= u[idx]

        prod_temp *= get_concentration(idx, idx_t, u, U_noint, n_integrated)
    end

    Jac[jac_term.i, jac_term.j] += jac_term.prefac * K_matrix[jac_term.idx_k, idx_t] * prod_temp
end
