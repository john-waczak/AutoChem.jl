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




function get_derivative_terms(rxn::BimolecularReaction, idx_k::Int, idx_integrated)
    dts = []

    # loop over reactants
    for idx_du ∈ rxn.reactants
        if idx_du ∈ idx_integrated
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
    end

    # loop over prodcuts
    for idx_du ∈ rxn.products
        if idx_du ∈ idx_integrated
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
    end

    return dts
end


function get_derivative_terms(rxn::TrimolecularReaction, idx_k::Int, idx_integrated)
    dts = []

    # loop over reactants
    for idx_du ∈ rxn.reactants
        if idx_du ∈ idx_integrated
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
    end

    # loop over prodcuts
    for idx_du ∈ rxn.products
        if idx_du ∈ idx_integrated
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
    end

    return dts
end



function get_derivative_terms(rxn::FittedPhotolysisReaction, idx_k::Int, idx_integrated)
    dts = []

    # loop over reactants
    for idx_du ∈ rxn.reactants
        if idx_du ∈ idx_integrated
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
    end

    # loop over prodcuts
    for idx_du ∈ rxn.products
        if idx_du ∈ idx_integrated
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
    end

    return dts
end



function get_bimolecular_derivatives(db, df_species)
    dts = BimolecularDerivativeTerm[]

    idx_integrated = findall(df_species.is_integrated .== 1)

    for i ∈ 1:length(db)
        terms = get_derivative_terms(db[i], i, idx_integrated)
        for term ∈ terms
            push!(dts, term)
        end
    end

    return dts
end


function get_trimolecular_derivatives(db, df_species)
    dts = TrimolecularDerivativeTerm[]

    idx_integrated = findall(df_species.is_integrated .== 1)

    for i ∈ 1:length(db)
        terms = get_derivative_terms(db[i], i, idx_integrated)
        for term ∈ terms
            push!(dts, term)
        end
    end

    return dts
end

function get_photolysis_derivatives(db, df_species)
    dts = PhotolysisDerivativeTerm[]

    idx_integrated = findall(df_species.is_integrated .== 1)

    for i ∈ 1:length(db)
        terms = get_derivative_terms(db[i], i, idx_integrated)
        for term ∈ terms
            push!(dts, term)
        end
    end

    return dts
end


# 21 μs
function update_k_bimol!(k_bimol, db_bimol, T, P, Mnow)
    Threads.@threads for i ∈ eachindex(k_bimol)
        @inbounds k_bimol[i] = db_bimol[i](T, P, Mnow)
    end
end

# 14 μs
function update_k_trimol!(k_trimol, db_trimol, T, P, Mnow)
    Threads.@threads for i ∈ eachindex(k_trimol)
        @inbounds k_trimol[i] = db_trimol[i](T, P, Mnow)
    end
end

# 182 μs
function update_k_photo!(k_photo, db_, T, P, Is)
    Threads.@threads for i ∈ eachindex(k_photo)
        @inbounds k_photo[i] = db_photo[i](T, P, Is)
    end
end





# NOTE: we use round() here to match calculation of pseudo observations

function get_time_index(t::Float64, Δt_step::Float64, tmin::Float64)
    return round(Int, (t-tmin + Δt_step)/Δt_step)
end


function get_concentration(idx, idx_t, u, U_noint, n_integrated)
    return (idx > n_integrated) ? U_noint[idx-n_integrated, idx_t] : u[idx]
end



function update_derivative!(idx_t::Int,
                            du::Vector{Float64},
                            u::Vector{Float64},
                            deriv_term::BimolecularDerivativeTerm,
                            K_matrix::Matrix{Float64},
                            prod_temp::Float64,
                            U_noint::Matrix{Float64},
                            n_integrated::Int
                            )

    for idx ∈ deriv_term.idxs_in
        prod_temp *= get_concentration(idx, idx_t, u, U_noint, n_integrated)
    end

    du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[deriv_term.idx_k, idx_t] * prod_temp
end



function update_derivative!(idx_t::Int,
                            du::Vector{Float64},
                            u::Vector{Float64},
                            deriv_term::TrimolecularDerivativeTerm,
                            K_matrix::Matrix{Float64},
                            prod_temp::Float64,
                            U_noint::Matrix{Float64},
                            n_integrated::Int
                            )

    for idx ∈ deriv_term.idxs_in
        prod_temp *= get_concentration(idx, idx_t, u, U_noint, n_integrated)
    end

    du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[deriv_term.idx_k, idx_t] * prod_temp
end


function update_derivative!(idx_t::Int,
                            du::Vector{Float64},
                            u::Vector{Float64},
                            deriv_term::PhotolysisDerivativeTerm,
                            K_matrix::Matrix{Float64},
                            prod_temp::Float64,
                            U_noint::Matrix{Float64},
                            n_integrated::Int
                            )

    for idx ∈ deriv_term.idxs_in
        prod_temp *= get_concentration(idx, idx_t, u, U_noint, n_integrated)
    end

    du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[deriv_term.idx_k, idx_t] * prod_temp
end



rhs_func = """

function rhs!(du, u, p, t)
    # get time value and index
    idx_t = get_time_index(t, Δt_step, ts[1])

    # set derivatives to zero
    du .= 0.0

    # loop over bimol derivatives
    prod_temp = 1.0
    @inbounds for i ∈ 1:length(derivatives_bimol)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives_bimol[i],
            K_bimol,
            prod_temp,
            U_noint,
            n_integrated
        )
    end

    # loop over trimol derivatives
    prod_temp = 1.0
    @inbounds for i ∈ 1:length(derivatives_trimol)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives_trimol[i],
            K_trimol,
            prod_temp,
            U_noint,
            n_integrated
        )
    end


    # loop over photolysis derivatives
    prod_temp = 1.0
    @inbounds for i ∈ 1:length(derivatives_photo)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives_photo[i],
            K_photo,
            prod_temp,
            U_noint,
            n_integrated
        )
    end
end


"""


function write_rhs_func(;model_name::String="autochem-w-ions")
    outpath = "./models/$(model_name)/mechanism/rhs.jl"

    # if it already exists, remove it so we can recreate it
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./models/$(model_name)")
        @warn "WARNING: Model directory, ./models/$(model_name), does not exists!"
        mkdir("./models/$(model_name)")
    end

    open(outpath, "w") do f
        println(f, rhs_func)
    end

end

