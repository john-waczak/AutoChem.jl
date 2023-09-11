abstract type RxnDerivative end

struct DerivativeTerm <: RxnDerivative
    idx_du::Int           # index of derivative term, i.e. du[i]
    idxs_in::Vector{Int}  # indices of reactants
    idx_k::Int            # reaction index
    prefac::Int           # ±1 (could use bool for storage)
end
