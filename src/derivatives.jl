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



function update_derivative!(idx_t::Int,
                            du::Vector{Float64},
                            u::Vector{Float64},
                            deriv_term::BimolecularDerivativeTerm,
                            K_matrix::Matrix{Float64},
                            prod_temp::Float64
                            )

    for i ∈ axes(deriv_term.idxs_in,1)
        prod_temp *= u[deriv_term.idxs_in[i]]
    end

    du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[deriv_term.idx_k, idx_t] * prod_temp
end

function update_derivative!(idx_t::Int,
                            du::Vector{Float64},
                            u::Vector{Float64},
                            deriv_term::TrimolecularDerivativeTerm,
                            K_matrix::Matrix{Float64},
                            prod_temp::Float64
                            )

    for i ∈ axes(deriv_term.idxs_in,1)
        prod_temp *= u[deriv_term.idxs_in[i]]
    end

    du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[deriv_term.idx_k, idx_t] * prod_temp
end


function update_derivative!(idx_t::Int,
                            du::Vector{Float64},
                            u::Vector{Float64},
                            deriv_term::PhotolysisDerivativeTerm,
                            K_matrix::Matrix{Float64},
                            prod_temp::Float64
                            )

    for i ∈ axes(deriv_term.idxs_in,1)
        prod_temp *= u[deriv_term.idxs_in[i]]
    end

    du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[deriv_term.idx_k, idx_t] * prod_temp
end




# function update_derivative!(idx_t::Int,
#                             du::Vector{Float64},
#                             u::Vector{Float64},
#                             deriv_term::DerivativeTermRO2,
#                             ro2_ratio::Float64,
#                             K_matrix::Matrix{Float64},
#                             Δt_step::Float64,
#                             prod_temp::Float64
#                             )

#     # prod_temp = u[deriv_term.idxs_in[1]]
#     # for i ∈ 2:size(deriv_term.idxs_in,1)
#     for i ∈ axes(deriv_term.idxs_in,1)
#         prod_temp *= u[deriv_term.idxs_in[i]]
#     end

#     du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[idx_t, deriv_term.idx_k] * prod_temp * ro2_ratio

#     # du[deriv_term.idx_du] += deriv_term.prefac * K_matrix[idx_t, deriv_term.idx_k] * prod(u[deriv_term.idxs_in]) * ro2_ratio
# end




# rhs_func = """
# function rhs!(du, u, p, t)
#     # set everything to sero
#     du .= 0.0

#     # get the time index
#     tval,idx_t = findmin(x -> abs.(x.- t), ts)

#     # get the current ro2_ratio
#     ro2_ratio = sum(u[idx_ro2])
#     ro2_ratio = ro2_ratio/RO2ᵢ
#     ro2_ratio = ro2_ratio > 0.0 ? ro2_ratio : 1.0  # make sure we don't have issues with ratio of 0

#     # set up product temporary value:
#     prod_temp = 1.0

#     # update derivatives
#     @inbounds for i ∈ 1:length(derivatives)
#         prod_temp = 1.0  # <-- need to start fresh for each derivative
#         update_derivative!(
#             idx_t,
#             du,
#             u,
#             derivatives[i],
#             ro2_ratio,
#             K_matrix,
#             Δt_step,
#             prod_temp
#         )
#     end

#     prod_temp = 1.0

#     @inbounds for i ∈ 1:length(derivatives_ro2)
#         prod_temp = 1.0  # <-- need to start fresh for each derivative
#         update_derivative!(
#             idx_t,
#             du,
#             u,
#             derivatives_ro2[i],
#             ro2_ratio,
#             K_matrix,
#             Δt_step,
#             prod_temp
#         )
#     end
# end
# """

# function write_rhs_func(;model_name::String="mcm")
#     outpath = "./models/$(model_name)/rhs.jl"

#     # if it already exists, remove it so we can recreate it
#     if isfile(outpath)
#         rm(outpath)
#     end

#     if !isdir("./models/$(model_name)")
#         mkdir("./models/$(model_name)")
#     end

#     open(outpath, "w") do f
#         println(f, rhs_func)
#     end

# end

