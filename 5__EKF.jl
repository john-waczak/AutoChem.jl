using AutoChem
using DelimitedFiles, CSV, DataFrames
using JSON
using LinearAlgebra
using DifferentialEquations
using Sundials, OrdinaryDiffEq
using Zygote, ForwardDiff, DiffResults
using SciMLSensitivity
using ParameterHandling, EarlyStopping
using StableRNGs
using BenchmarkTools
using ProgressMeter
using Measurements




# set up model output directory
collection_id = "empty"
unc_ext = "_std"
# model_name = "autochem-w-ions"
model_name = "methane"
model_path= "models"


@info "Setting up file paths..."
outpath = joinpath(model_path, model_name, "runs", collection_id)
@assert ispath(outpath)

docs_path = joinpath(model_path, model_name, "docs")
@assert ispath(docs_path)

if !ispath(joinpath(outpath, "figures"))
    joinpath(model_path, model_name, "figures")
end


if ! ispath(joinpath(outpath, "EKF"))
    mkpath(joinpath(outpath, "EKF"))
end



# fetch species list
@info "Loading data into DataFrames"

df_species = CSV.read(joinpath(outpath, "mechanism", "species.csv"), DataFrame);
# fetch data
df_params= CSV.read(joinpath(outpath, "mechanism", "state_parameters.csv"), DataFrame);
df_params_ϵ = CSV.read(joinpath(outpath, "mechanism", "state_parameters_ϵ.csv"), DataFrame);
df_nd = CSV.read(joinpath(outpath, "mechanism", "number_densities.csv"), DataFrame);
df_nd_ϵ = CSV.read(joinpath(outpath, "mechanism", "number_densities_ϵ.csv"), DataFrame);
df_ions = CSV.read(joinpath(outpath, "mechanism", "ions.csv"), DataFrame)
df_ions_ϵ = CSV.read(joinpath(outpath, "mechanism", "ions_ϵ.csv"), DataFrame)

df_u₀ = CSV.read(joinpath(outpath, "4d-var", "u0_final.csv"), DataFrame);
df_u0 = CSV.read(joinpath(outpath, "4d-var", "u0.csv"), DataFrame);
df_u0_orig = CSV.read(joinpath(outpath, "mechanism", "u0.csv"), DataFrame);


# set initial conditions
@info "Getting initial condition vector"
const u₀ = df_u₀.u0[:]
@assert !(any(u₀ .< 0.0))



# generate integrated vs non-integrated indices
@info "Generating species indices"
const idx_noint = findall(df_species.is_integrated .== 2)
const n_integrated = sum(df_species.is_integrated .== 1)
const idx_pos = get_positive_indices(df_species)
const idx_neg = get_negative_indices(df_species)
const Δt_step = 15.0


# read reaction databases
@info "Loading reaction databases"
const db_bimol = read_bimol(joinpath(outpath, "mechanism", "bimol.json"));
const db_trimol = read_trimol(joinpath(outpath, "mechanism", "trimol.json"));
const db_photo = read_fitted_photolysis(joinpath(outpath, "mechanism", "fitted_photolysis.json"));


# generate derivatives list
@info "Generating derivatives terms"
const derivatives_bimol = get_bimolecular_derivatives(db_bimol, df_species)
const derivatives_trimol = get_trimolecular_derivatives(db_trimol, df_species)
const derivatives_photo = get_photolysis_derivatives(db_photo, df_species)


# generate jacobian list
@info "Generating jacobian terms"
const jacobian_terms_bimol = get_bimolecular_jacobian_terms(derivatives_bimol)
const jacobian_terms_trimol = get_trimolecular_jacobian_terms(derivatives_trimol)
const jacobian_terms_photo = get_photolysis_jacobian_terms(derivatives_photo)


# read reaction rate coefficients
@info "Reading reaction rate coefficients"
const K_bimol = readdlm(joinpath(outpath, "mechanism", "K_bimol.csv"), ',')
const K_trimol = readdlm(joinpath(outpath, "mechanism", "K_trimol.csv"), ',')
const K_photo = readdlm(joinpath(outpath, "mechanism", "K_photo.csv"), ',')

# read temperature and pressure data
const temperatures = df_params.temperature
const pressures = df_params.pressure

# generate lookup table for non-integrated species
@info "generating non-integrated species lookup tables"
const U_noint = zeros(length(idx_noint), nrow(df_params))

noint_dict = Dict(
    "Ar" => Ar.(temperatures, pressures),
    "O2" => O2.(temperatures, pressures),
    "N2" => N2.(temperatures, pressures),
    "Photon" => ones(Float64, nrow(df_params)),
    "m" => M.(temperatures, pressures),
)

for i ∈ 1:length(idx_noint)
    U_noint[i,:] .= noint_dict[df_species.varname[idx_noint[i]]]
end


# generate global constants
@info "Generating measurement matrices"
#measurements_to_ignore = ["C2H6", "SO2", "t", "w_ap"]
measurements_to_ignore = ["t", "w_ap"]

df_nd_to_use = df_nd[:, Not(measurements_to_ignore)]
df_nd_to_use_ϵ = df_nd_ϵ[:, Not(measurements_to_ignore)]

const W = zeros(ncol(df_nd_to_use)+2, nrow(df_nd_to_use))
const meas_ϵ = zeros(ncol(df_nd_to_use_ϵ)+2, nrow(df_nd_to_use_ϵ))

W[1:end-2, :] .= Matrix(df_nd_to_use)'
W[end-1:end, :] .= Matrix(df_ions)'

meas_ϵ[1:end-2, :] .= Matrix(df_nd_to_use_ϵ)'
meas_ϵ[end-1:end, :] .= Matrix(df_ions_ϵ)'

# reduce uncertainty in ion counts to force assimilation to adhere to it more
# meas_ϵ[end-1:end,:] .= 0.1 .* meas_ϵ[end-1:end,:]


@info "generating measurement indices"
idx_measurements = []
for meas ∈ names(df_nd_to_use)
    println(meas)
    idx = findfirst(df_species.varname .== meas)
    push!(idx_measurements, idx)
end

const idx_meas = idx_measurements


# generating other global constants
@info "Generating other global constants"

const ts = df_params.t

const fudge_fac::Float64 = 0.5
# const fudge_fac::Float64 = 1.0

const tmin::Float64 = minimum(ts)
const tmax::Float64 = maximum(ts)
const abstol::Float64 = 1e-3
const reltol::Float64 = 1e-3
const tspan = (tmin, tmax)

const ϵ::Float64 = 1
const ϵ_min::Float64 = 1e-12


# load rhs and jacobian functions
@info "Loading rhs and jac functions"
include(joinpath(outpath, "mechanism", "rhs.jl"))
include(joinpath(outpath, "mechanism", "jacobian.jl"))



# Set up sensitivity Algorithm
@info "Choosing sensitivity algorithm and solver"
sensealg_dict = Dict(
    :QuadratureAdjoint => QuadratureAdjoint(),
    :BacksolveAdjoint => BacksolveAdjoint(),
    :InterpolatingAdjoint => InterpolatingAdjoint(),
    :ZygoteAdjoint => ZygoteAdjoint(),
    :ForwardDiffSensitivity => ForwardDiffSensitivity()
)

sensealg = sensealg_dict[:ForwardDiffSensitivity]



# idx_0 = findfirst(ts .== 0) + 1
idx_0 = 1



# define the ODE function
@info "Defining ODE function"
fun = ODEFunction(rhs!; jac=jac!)

# we need to figure out if starting at a different time will cause problems...
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀ , (ts[idx_0], ts[end]))
# ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun,df_u0_orig.u0, tspan)

@info "Trying a solve with default u₀"
@benchmark solve(ode_prob; alg_hints=[:stiff], saveat=Δt_step, reltol=reltol, abstol=abstol)
sol = solve(ode_prob; alg_hints=[:stiff], saveat=Δt_step, reltol=reltol, abstol=abstol)

# set up observation observation operator and it's jacobian
@info "Testing observation operator"
Obs(u₀, idx_meas, idx_pos, idx_neg)
JObs(u₀, idx_meas, idx_pos, idx_neg)

@info "Testing Observation Covariance Matrix and its Inverse"
Rmat(1, meas_ϵ; fudge_fac=fudge_fac)
Rinv(1, meas_ϵ; fudge_fac=fudge_fac)


# Setting up matrices for EKF
@info "Pre-allocating matrices for EKF"
const P::Matrix{Float64} = zeros(nrow(df_species), nrow(df_species))
const P_diag::Matrix{Float64} = zeros(nrow(df_species), length(ts)) # i.e. P_diag[i] == P[i,i]
const Q::Matrix{Float64} = zeros(size(P))
const uₐ::Matrix{Float64} = zeros(length(u₀), length(ts))



# Initialize Covariance Matrices
@info "Initializing covariance matrices"

for i ∈ 1:length(u₀)
    if i ∈ idx_pos
        P[i,i] = (ϵ * u₀[i])^2 + ϵ_min^2
        #P[i,i] = (0.1 * ϵ * u₀[i])^2 + ϵ_min^2
        #P[i,i] = (10 * ϵ * u₀[i])^2 + ϵ_min^2
    elseif i ∈ idx_neg
        P[i,i] = (ϵ * u₀[i])^2 + ϵ_min^2
        #P[i,i] = (0.1 * ϵ * u₀[i])^2 + ϵ_min^2
        #P[i,i] = (10 * ϵ * u₀[i])^2 + ϵ_min^2
    else
        P[i,i] = (ϵ * u₀[i])^2 + ϵ_min^2
    end

    P_diag[i,1] = P[i,i]
end


# initially, set Q to match P
Q .= P


# set the first value to the background estimate
#uₐ[:,1] .= u₀
uₐ[:,idx_0] .= u₀

# Establish forward model function
@info "Testing model propagator"

function model_forward(u_now, t_now)
    _prob = remake(ode_prob, u0=u_now, tspan=(t_now, t_now+Δt_step))
    solve(_prob; alg_hints=[:stiff], reltol=reltol, abstol=abstol, dense=false,save_everystep=false,save_start=false, sensealg=sensealg)[:,end]
end


@info "Testing out Jacobian Determination"

@benchmark model_forward(u₀, ts[1])


res = DiffResults.JacobianResult(u₀);
ForwardDiff.jacobian!(res, u->model_forward(u,ts[1]), u₀);  # ~ 10 s


# @benchmark Zygote.withjacobian(model_forward, u₀, ts[1])
# u_test, DM_test = Zygote.withjacobian(model_forward, u₀, ts[1])

res.value
res.derivs[1]


u_test = res.value
DM_test = res.derivs[1]

# Perform Assimilation Loop

const u_now::Vector{Float64} = zeros(size(uₐ[:,1]))
const DM::Matrix{Float64} = DM_test

@info "Starting assimilation..."

any(isnan.(W))



#@showprogress for k ∈ 1:length(ts)-1  # because we always update the *next* value
@showprogress for k ∈ idx_0:length(ts)-1  # because we always update the *next* value
    # k=1

    # collect current model estimate
    u_now .= uₐ[:,k]  # should preallocate this

    # --------------------------
    # Forecast Step
    # --------------------------

    # run model forward one Δt

    # u_next, DM_tup = Zygote.withjacobian(model_forward, u_now, ts[k])  # <-- can't make this mutating
    # DM .= DM_tup[1]  # result is a 1-element tuple, so we index it
    ForwardDiff.jacobian!(res, u->model_forward(u, ts[k]), u_now);

    u_next = res.value
    DM .= res.derivs[1]

    # collect observations
    u_h = Obs(u_next, idx_meas, idx_pos, idx_neg)

    # update loop for the model covariance matrix, Q
    if k > 1
        Q .= 0.0
        for j ∈ axes(Q,1)
            # define error growth rate
            eg = abs(u_next[j]-u_now[j])/max(1.0, u_next[j])  # note: we should actually do the sens. analysis here

            # clip the growth rate to 1%-50% range
            eg = max(min(eg, 0.5), 0.01)


            # if we have a measurement, further manipulate the growth rate
            if j ∈ idx_meas #_nonan
                #u_h_j = u_h[idx_meas_nonan .== j][1]
                u_h_j = u_h[findfirst(idx_meas .== j)]
                diff = abs(u_next[j] - u_h_j)
                ratio = diff/u_h_j

                if ratio > 1
                    eg = max(eg, 0.5)
                elseif ratio > 0.5
                    eg = max(eg, 0.2)
                elseif ratio > 0.2
                    eg = max(eg, 0.1)
                else
                    continue
                end
            end

            # finally, update Q
            Q[j,j] = eg * u_next[j]^2 + ϵ_min^2
        end
    end


    # update the background covariance matrix
    P .= DM*P*DM' + Q

    # --------------------------
    # Analysis Step
    # --------------------------

    DH = JObs(u_next, idx_meas, idx_pos, idx_neg)

    # R = Rmat_nonan(k+1, is_meas_not_nan, meas_ϵ; fudge_fac=fudge_fac)
    R = Rmat(k+1, meas_ϵ; fudge_fac=fudge_fac)

    denom = DH*P*DH' + R

    #Kalman = P*DH'/denom
    Kalman = P*DH'* inv(denom)  # <-- this seems slightly faster for now

    u_next .= u_next .+ Kalman*(W[:, k+1] - u_h)

    P .= (I(length(u_next)) - Kalman*DH)*P


    # do some clipping to make sure we stay reasonable
    for j ∈ axes(P,1)
        P[j,j] = min(max(0, u_next[j]^2), P[j,j])
        P[j,j] = max(0.05*u_next[j]^2, P[j,j])
    end

    # filter negative values to zero
    u_next[u_next .≤ 0.0] .= 0.0

    # update the analysis vector
    #update_uₐ!(uₐ, u_next, k+1)
    uₐ[:,k+1] .= u_next
    P_diag[:,k+1] .= [P[i,i] for i ∈ axes(P,1)]
end



# --------------------------------------------------------------------------------------------------------------------------
# 12. Post Process
# --------------------------------------------------------------------------------------------------------------------------



# combine concentration with uncertainty
uₐ_nd = uₐ[:,idx_0:end] .± sqrt.(P_diag[:,idx_0:end])

# convert final output into mixing ratios
d = df_params.M[idx_0:end] .± (fudge_fac .* df_params_ϵ.M[idx_0:end])

uₐ_mr = to_mixing_ratio(uₐ_nd, d)

# chop off values and uncertainties for easier plotting
ua_mr_vals = Measurements.value.(uₐ_mr)
ua_mr_ϵ = Measurements.uncertainty.(uₐ_mr)


# save to output file
df_ekf = DataFrame()
df_ekf_ϵ = DataFrame()



@showprogress for i ∈ axes(ua_mr_vals, 1)
    df_ekf[!, df_species[i, "printname"]] = ua_mr_vals[i,:]
    df_ekf_ϵ[!, df_species[i, "printname"]] = ua_mr_ϵ[i,:]
end

df_ekf[!, :times] = ts[idx_0:end]
df_ekf_ϵ[!, :times] = ts[idx_0:end]

CSV.write(joinpath(outpath, "EKF", "ekf_output.csv"), df_ekf)
CSV.write(joinpath(outpath, "EKF", "ekf_ϵ_output.csv"), df_ekf_ϵ)


# combine measurements with uncertainties
W_mr = W[:,idx_0:end] .± (fudge_fac .* meas_ϵ[:,idx_0:end])
W_mr = to_mixing_ratio(W_mr, d)


W_mr_val = Measurements.value.(W_mr)
W_mr_ϵ = Measurements.uncertainty.(W_mr)

# save measurements to csv files for final output
size(W_mr_val)
idx_meas

df_w = DataFrame()
df_w_ϵ = DataFrame()


@showprogress for i ∈ 1:length(idx_meas)
    df_w[!, df_species[idx_meas[i], "printname"]] = W_mr_val[i,:]
    df_w_ϵ[!, df_species[idx_meas[i], "printname"]] = W_mr_ϵ[i,:]
end


df_w[!,:times] = ts[idx_0:end]
df_w_ϵ[!,:times] = ts[idx_0:end]

CSV.write(joinpath(outpath, "EKF", "ekf_measurements.csv"), df_w)
CSV.write(joinpath(outpath, "EKF", "ekf_measurements_ϵ.csv"), df_w_ϵ)



# add ion data separately since we measure these as :

df_ions = DataFrame()
df_ions_ϵ = DataFrame()

pos_tot = sum(uₐ_nd[idx_pos, :], dims=1)[:]
pos_tot_val = Measurements.value.(pos_tot)
pos_tot_ϵ= Measurements.uncertainty.(pos_tot)

neg_tot = sum(uₐ_nd[idx_neg, :], dims=1)[:]
neg_tot_val = Measurements.value.(neg_tot)
neg_tot_ϵ= Measurements.uncertainty.(neg_tot)

df_ions[!, "Total Positive Ions (modeled)"] = pos_tot_val
df_ions[!, "Total Positive Ions (measured)"] = W[end-1,idx_0:end]
df_ions[!, "Total Negative Ions (modeled)"] = neg_tot_val
df_ions[!, "Total Negative Ions (measured)"] = W[end,idx_0:end]

df_ions_ϵ[!, "Total Positive Ions (modeled)"] = pos_tot_ϵ
df_ions_ϵ[!, "Total Positive Ions (measured)"] = meas_ϵ[end-1,idx_0:end]
df_ions_ϵ[!, "Total Negative Ions (modeled)"] = neg_tot_ϵ
df_ions_ϵ[!, "Total Negative Ions (measured)"] = meas_ϵ[end,idx_0:end]

CSV.write(joinpath(outpath, "EKF", "ions.csv"), df_ions)
CSV.write(joinpath(outpath, "EKF", "ions_ϵ.csv"), df_ions_ϵ)




# --------------------------------------------------------------------------------------------------------------------------
# 14. Compute Lifetime
# --------------------------------------------------------------------------------------------------------------------------

# we write the lifetime for species i as

# τᵢ = ∑ⱼτᵢⱼ where j are all reactions for which i is a reactant

# the form for each τᵢⱼ depends on the specific reaction type

# Photodissociation reaction:
# X + hν ⟶ products
# Ẋ = -jX  [molecules/cm³/s]
# τ = 1\j  [s]

# Collision reaction:
# X + ∑ⱼYⱼ ⟶ products
# Ẋ = -kX⋅ΠⱼYⱼ   [molecules/cm³/s]
# τ = 1\(kΠⱼYⱼ)
#     or
# τ = X\(kXΠⱼYⱼ)

# now we need to combine
τs = ones(size(uₐ))  # preallocate matrix to hold values
ℓ_mat = zeros(size(uₐ))  # loss rate
#ℓ = 1.0

size(uₐ)
size(uₐ_nd)

K_bimol_view = @view K_bimol[:, idx_0:end]
K_trimol_view = @view K_trimol[:, idx_0:end]
K_photo_view = @view K_photo[:, idx_0:end]

@showprogress for d ∈ 1:length(derivatives_bimol)
    derivative = derivatives_bimol[d]
    if derivative.prefac < 0.0
        for idx_t ∈ axes(K_bimol_view, 2)
            ℓ = K_bimol_view[derivative.idx_k, idx_t]
            for i ∈ derivative.idxs_in
                if i != derivative.idx_du
                    if i ∈ idx_noint
                        ℓ *= U_noint[i-n_integrated,idx_0+idx_t-1]
                    else
                        ℓ *= uₐ[i, idx_t]
                    end
                end
            end
            ℓ_mat[derivative.idx_du, idx_t] += ℓ
        end
    end
end


@showprogress for d ∈ 1:length(derivatives_trimol)
    derivative = derivatives_trimol[d]
    if derivative.prefac < 0.0
        for idx_t ∈ axes(K_trimol_view, 2)
            ℓ = K_trimol_view[derivative.idx_k, idx_t]
            for i ∈ derivative.idxs_in
                if i != derivative.idx_du
                    if i ∈ idx_noint
                        ℓ *= U_noint[i-n_integrated,idx_0+idx_t-1]
                    else
                        ℓ *= uₐ[i, idx_t]
                    end
                end
            end

            ℓ_mat[derivative.idx_du, idx_t] += ℓ
        end
    end
end


@showprogress for d ∈ 1:length(derivatives_photo)
    derivative = derivatives_photo[d]
    if derivative.prefac < 0.0
        for idx_t ∈ axes(K_trimol_view, 2)
            ℓ = K_photo_view[derivative.idx_k, idx_t]
            ℓ_mat[derivative.idx_du, idx_t] += ℓ
        end
    end
end


τs
for j ∈ axes(τs, 2), i ∈ axes(τs,1)
    # if isinf(τs[i,j]/ℓ_mat[i,j] ) || isnan(τs[i,j]/ℓ_mat[i,j] )
    #     println("idx: ", (i,j), "\tu:\t", τs[i,j], "\tℓ:\t",ℓ_mat[i,j], "\tτ:\t", τs[i,j]/ℓ_mat[i,j] )
    # end
    τs[i,j] = τs[i,j] / ℓ_mat[i,j]
end

# heatmap(τs)


# generate lifetime dataframes for output
df_τs = DataFrame()
@showprogress for i ∈ axes(τs, 1)
    df_τs[!, df_species[i, "printname"]] = τs[i,:]
end

df_τs[!, :times] = ts[idx_0:end]

CSV.write(joinpath(outpath, "EKF", "lifetimes.csv"), df_τs)

df_τs_means = DataFrame()
@showprogress for i ∈ axes(τs, 1)
    df_τs_means[!, df_species[i, "printname"]] = [mean(τs[i,:])]
end

df_τs_means

τs_means = Matrix(df_τs_means)
idx_sort = sortperm(τs_means, dims=2, rev=true)

spec_name = []
mean_lifetime = []
for idx ∈ idx_sort
    push!(spec_name, names(df_τs_means)[idx])
    push!(mean_lifetime,df_τs_means[1, idx] )
end

df_τs_means_sorted = DataFrame(:species => spec_name, :τ_seconds => mean_lifetime)

df_τs_means_sorted.τ_minutes = (df_τs_means_sorted.τ_seconds ./ 60)
df_τs_means_sorted.τ_hours = (df_τs_means_sorted.τ_minutes ./ 60)
df_τs_means_sorted.τ_days = (df_τs_means_sorted.τ_hours ./ 24)
df_τs_means_sorted.τ_weeks = (df_τs_means_sorted.τ_days ./7)
df_τs_means_sorted.τ_years = (df_τs_means_sorted.τ_days ./365)


idx_good = []
idx_bad = []
for i ∈ 1:nrow(df_τs_means_sorted)
    row = df_τs_means_sorted[i, :]
    if isnan(row.τ_seconds) || isinf(row.τ_seconds)
        push!(idx_bad, i)
    else
        push!(idx_good, i)
    end
end


df_τs_means_sorted = df_τs_means_sorted[idx_good, :]

CSV.write(joinpath(outpath, "EKF", "mean_lifetimes.csv"), df_τs_means_sorted)


df_τs_means_sorted


df_τs_means_sorted[35, :]

# --------------------------------------------------------------------------------------------------------------------------
# 15. Compute Production and Destruction Fractions
# --------------------------------------------------------------------------------------------------------------------------
du = zeros(size(uₐ,1))

du_prod = zeros(size(uₐ))
du_dest = zeros(size(uₐ))
du_bi = zeros(length(derivatives_bimol), size(uₐ,2))
du_tri = zeros(length(derivatives_trimol), size(uₐ,2))
du_photo = zeros(length(derivatives_photo), size(uₐ,2))

# we now need to loop through all the times and recompute these values
# we probably need new versions of "update derivative" that will also update the matrices

ts[idx_0:end]



function rhs2!(du, u, p, t, du_bi, du_tri, du_photo)
    # get time value and index
    idx_t = get_time_index(t, Δt_step, ts[1])

    # set derivatives to zero
    du .= 0.0

    # loop over bimol derivatives
    prod_temp = 1.0
    @inbounds for i ∈ 1:length(derivatives_bimol)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_derivative!(
            i,
            idx_t,
            du,
            u,
            derivatives_bimol[i],
            K_bimol,
            prod_temp,
            U_noint,
            n_integrated,
            du_bi,
            du_tri,
            du_photo
        )
    end

    # loop over trimol derivatives
    prod_temp = 1.0
    @inbounds for i ∈ 1:length(derivatives_trimol)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_derivative!(
            i,
            idx_t,
            du,
            u,
            derivatives_trimol[i],
            K_trimol,
            prod_temp,
            U_noint,
            n_integrated,
            du_bi,
            du_tri,
            du_photo
        )
    end


    # loop over photolysis derivatives
    prod_temp = 1.0
    @inbounds for i ∈ 1:length(derivatives_photo)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_derivative!(
            i,
            idx_t,
            du,
            u,
            derivatives_photo[i],
            K_photo,
            prod_temp,
            U_noint,
            n_integrated,
            du_bi,
            du_tri,
            du_photo
        )
    end

    nothing
end


@benchmark rhs2!(du, uₐ[:,idx_0], nothing, ts[idx_0], du_bi, du_tri, du_photo)

for k ∈ idx_0:length(ts)
    rhs2!(du, uₐ[:,k], nothing, ts[k], du_bi, du_tri, du_photo)
end

K_photo

du_bi

# convert these to production and loss rates in units of Hz
# f_bi = copy(du_bi)
# f_tri = copy(du_tri)
# f_photo = copy(du_photo)

# # loop through reactions and divide by concentration
# size(derivatives_bimol)
# for i ∈ 1:length(derivatives_bimol)
#     idx_du = derivatives_bimol[i].idx_du
#     # f_bi[i,:] .= f_bi[i,:] ./ uₐ[idx_du,:]
#     f_bi[i,:] .= f_bi[i,:]
# end

# for i ∈ 1:length(derivatives_trimol)
#     idx_du = derivatives_trimol[i].idx_du
#     # f_tri[i,:] .= f_tri[i,:] ./ uₐ[idx_du,:]
#     f_tri[i,:] .= f_tri[i,:]
# end


# for i ∈ 1:length(derivatives_photo)
#     idx_du = derivatives_photo[i].idx_du
#     # f_photo[i,:] .= f_photo[i,:] ./ uₐ[idx_du,:]
#     f_photo[i,:] .= f_photo[i,:]
# end

# maximum(du)
# minimum(du)

# writedlm(joinpath(outpath, "EKF", "freq_bi.csv"), f_bi[:, idx_0:end], ",")
# writedlm(joinpath(outpath, "EKF", "freq_tri.csv"), f_tri[:, idx_0:end], ",")
# writedlm(joinpath(outpath, "EKF", "freq_photo.csv"), f_photo[:, idx_0:end], ",")

writedlm(joinpath(outpath, "EKF", "du_bi.csv"), du_bi[:, idx_0:end], ",")
writedlm(joinpath(outpath, "EKF", "du_tri.csv"), du_tri[:, idx_0:end], ",")
writedlm(joinpath(outpath, "EKF", "du_photo.csv"), du_photo[:, idx_0:end], ",")



