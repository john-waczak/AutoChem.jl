abstract type RxnDerivative end

struct BimolecularDerivativeTerm <: RxnDerivative
    idx_du::Int           # index of species derivative term, i.e. du[i]
    idxs_in::Vector{Int}  # indices of reactants
    idx_k::Int            # reaction index
    prefac::Int           # ±1 (could use bool for storage)
end

struct TrimolecularDerivativeTerm <: RxnDerivative
    idx_du::Int           # index of species derivative term, i.e. du[i]
    idxs_in::Vector{Int}  # indices of reactants
    idx_k::Int            # reaction index
    prefac::Int           # ±1 (could use bool for storage)
end

struct PhotolysisDerivativeTerm <: RxnDerivative
    idx_du::Int           # index of species derivative term, i.e. du[i]
    idxs_in::Vector{Int}  # indices of reactants
    idx_k::Int            # reaction index
    prefac::Int           # ±1 (could use bool for storage)
end




function get_derivative_terms(rxn::BimolecularReaction, idx_k::Int)
    dts = []

    # loop over reactants
    for idx_du ∈ rxn.reactants
        push!(
            dts,
            BimolecularDerivativeTerm(
                idx_du,
                rxn.reactants,
                idx_k,
                -1
            )
        )
    end

    # loop over prodcuts
    for idx_du ∈ rxn.products
        push!(
            dts,
            BimolecularDerivativeTerm(
                idx_du,
                rxn.reactants,
                idx_k,
                1
            )
        )
    end

    return dts
end


function get_derivative_terms(rxn::TrimolecularReaction, idx_k::Int)
    dts = []

    # loop over reactants
    for idx_du ∈ rxn.reactants
        push!(
            dts,
            TrimolecularDerivativeTerm(
                idx_du,
                rxn.reactants,
                idx_k,
                -1
            )
        )
    end

    # loop over prodcuts
    for idx_du ∈ rxn.products
        push!(
            dts,
            TrimolecularDerivativeTerm(
                idx_du,
                rxn.reactants,
                idx_k,
                1
            )
        )
    end

    return dts
end



function get_derivative_terms(rxn::FittedPhotolysisReaction, idx_k::Int)
    dts = []

    # loop over reactants
    for idx_du ∈ rxn.reactants
        push!(
            dts,
            PhotolysisDerivativeTerm(
                idx_du,
                rxn.reactants,
                idx_k,
                -1
            )
        )
    end

    # loop over prodcuts
    for idx_du ∈ rxn.products
        push!(
            dts,
            PhotolysisDerivativeTerm(
                idx_du,
                rxn.reactants,
                idx_k,
                1
            )
        )
    end

    return dts
end



function get_bimolecular_derivatives(db)
    dts = BimolecularDerivativeTerm[]

    for i ∈ 1:length(db)
        terms = get_derivative_terms(db[i], i)
        for term ∈ terms
            push!(dts, term)
        end
    end

    return dts
end


function get_trimolecular_derivatives(db)
    dts = TrimolecularDerivativeTerm[]

    for i ∈ 1:length(db)
        terms = get_derivative_terms(db[i], i)
        for term ∈ terms
            push!(dts, term)
        end
    end

    return dts
end

function get_photolysis_derivatives(db)
    dts = PhotolysisDerivativeTerm[]

    for i ∈ 1:length(db)
        terms = get_derivative_terms(db[i], i)
        for term ∈ terms
            push!(dts, term)
        end
    end

    return dts
end

