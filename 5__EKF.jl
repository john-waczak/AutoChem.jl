using AutoChem

using DelimitedFiles, CSV, DataFrames
using JSON

using LinearAlgebra

using DifferentialEquations
using Sundials, OrdinaryDiffEq
using Zygote, ForwardDiff
using SciMLSensitivity


using ParameterHandling, EarlyStopping

using StableRNGs

using BenchmarkTools

using ProgressMeter

using CairoMakie
using MintsMakieRecipes

set_theme!(mints_theme)
update_theme!(
    figure_padding=30,
    Axis = (
        xticklabelsize=20,
        yticklabelsize=20,
        xlabelsize=22,
        ylabelsize=22,
        titlesize=25,
    ),
    Colorbar = (
        ticklabelsize=20,
        labelsize=22
    )
)



# set up model output directory
collection_id = "empty"
unc_ext = "_std"
model_name = "autochem-w-ions"
# model_name = "qroc-methane-intel"
model_path= "models"

if ! ispath(joinpath(model_path, model_name, "EKF"))
    mkpath(joinpath(model_path, model_name, "EKF"))
end



# fetch species list
@info "Loading data into DataFrames"

df_species = CSV.read(joinpath(model_path, model_name, "mechanism", "species.csv"), DataFrame);
# fetch data
df_params= CSV.read(joinpath(model_path, model_name, "mechanism", "state_parameters.csv"), DataFrame);
df_params_ϵ = CSV.read(joinpath(model_path, model_name, "mechanism", "state_parameters_ϵ.csv"), DataFrame);
df_nd = CSV.read(joinpath(model_path, model_name, "mechanism", "number_densities.csv"), DataFrame);
df_nd_ϵ = CSV.read(joinpath(model_path, model_name, "mechanism", "number_densities_ϵ.csv"), DataFrame);
df_ions = CSV.read(joinpath(model_path, model_name, "mechanism", "ions.csv"), DataFrame)
df_ions_ϵ = CSV.read(joinpath(model_path, model_name, "mechanism", "ions_ϵ.csv"), DataFrame)

df_u₀ = CSV.read(joinpath(model_path, model_name, "4d-var", "u0_final.csv"), DataFrame);
df_u0 = CSV.read(joinpath(model_path, model_name, "4d-var", "u0.csv"), DataFrame);
df_u0_orig = CSV.read(joinpath(model_path, model_name, "mechanism", "u0.csv"), DataFrame);

# set initial conditions
@info "Getting initial condition vector"
const u₀ = df_u₀.u0[:]



# generate integrated vs non-integrated indices
@info "Generating species indices"
const idx_noint = findall(df_species.is_integrated .== 2)
const n_integrated = sum(df_species.is_integrated .== 1)
const idx_pos = get_positive_indices(df_species)
const idx_neg = get_negative_indices(df_species)
const Δt_step = 15.0


# read reaction databases
@info "Loading reaction databases"
const db_bimol = read_bimol(joinpath(model_path, model_name, "mechanism", "bimol.json"));
const db_trimol = read_trimol(joinpath(model_path, model_name, "mechanism", "trimol.json"));
const db_photo = read_fitted_photolysis(joinpath(model_path, model_name, "mechanism", "fitted_photolysis.json"));


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
const K_bimol = readdlm(joinpath(model_path, model_name, "mechanism", "K_bimol.csv"), ',')
const K_trimol = readdlm(joinpath(model_path, model_name, "mechanism", "K_trimol.csv"), ',')
const K_photo = readdlm(joinpath(model_path, model_name, "mechanism", "K_photo.csv"), ',')

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
measurements_to_ignore = ["C2H6", "SO2", "t", "w_ap"]

df_nd_to_use = df_nd[:, Not(measurements_to_ignore)]
df_nd_to_use_ϵ = df_nd_ϵ[:, Not(measurements_to_ignore)]

const W = zeros(ncol(df_nd_to_use)+2, nrow(df_nd_to_use))
const meas_ϵ = zeros(ncol(df_nd_to_use_ϵ)+2, nrow(df_nd_to_use_ϵ))

W[1:end-2, :] .= Matrix(df_nd_to_use)'
W[end-1:end, :] .= Matrix(df_ions)'

meas_ϵ[1:end-2, :] .= Matrix(df_nd_to_use_ϵ)'
meas_ϵ[end-1:end, :] .= Matrix(df_ions_ϵ)'


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

ts = df_params.t

#const fudge_fac::Float64 = 0.1
#const fudge_fac::Float64 = 0.5
const fudge_fac::Float64 = 1.0

const tmin::Float64 = minimum(ts)
const tmax::Float64 = maximum(ts)
const abstol::Float64 = 1e-3
const reltol::Float64 = 1e-3
const tspan = (tmin, tmax)

const ϵ::Float64 = 0.5
const ϵ_min::Float64 = 1e-12


# load rhs and jacobian functions
@info "Loading rhs and jac functions"
include("models/$model_name/mechanism/rhs.jl")
include("models/$model_name/mechanism/jacobian.jl")



# Set up sensitivity Algorithm
@info "Choosing sensitivity algorithm and solver"
sensealg_dict = Dict(
    :QuadratureAdjoint => QuadratureAdjoint(),
    :BacksolveAdjoint => BacksolveAdjoint(),
    :InterpolatingAdjoint => InterpolatingAdjoint(),
    :ZygoteAdjoint => ZygoteAdjoint()
)

sensealg = sensealg_dict[:QuadratureAdjoint]
solver = TRBDF2()

# define the ODE function
@info "Defining ODE function"
fun = ODEFunction(rhs!; jac=jac!)

# we need to figure out if starting at a different time will cause problems...
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀ , tspan)
# ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun,df_u0_orig.u0, tspan)

@info "Trying a solve with default u₀"
# sol = solve(ode_prob, QNDF(); saveat=Δt_step, reltol=reltol, abstol=abstol)
# sol = solve(ode_prob; saveat=Δt_step, reltol=reltol, abstol=abstol)
# sol = solve(ode_prob, CVODE_BDF(); saveat=Δt_step, reltol=reltol, abstol=abstol)
# sol = solve(ode_prob, QNDF(); saveat=Δt_step, reltol=reltol, abstol=abstol)

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



# # Set up loss function for 4d-var
# @info "Setting up loss function..."
# const u0b::Vector{Float64} = copy(u0a) # i.e. "background guess"
# const B::Matrix{Float64} = diagm((ϵ .* (u₀)) .^2  .+ ϵ_min^2)
# const Binv::Matrix{Float64} = inv(B)

# const use_background_cov::Bool = false



# Initialize Covariance Matrices
@info "Initializing covariance matrices"

for i ∈ 1:length(u₀)
    P[i,i] = (ϵ * u₀[i])^2 + ϵ_min^2
    P_diag[i,1] = P[i,i]
end

# initially, set Q to match P
Q .= P


# set the first value to the background estimate
uₐ[:,1] .= u₀


# Establish forward model function
@info "Testing model propagator"
function model_forward(u_now, t_now)
    _prob = remake(ode_prob, u0=u_now, tspan=(t_now, t_now+Δt_step))
    solve(_prob, TRBDF2(), reltol=reltol, abstol=abstol, dense=false,save_everystep=false,save_start=false, sensealg=sensealg)[:,end]
end

@info "Testing out Jacobian Determination"
u_test, DM_test = Zygote.withjacobian(model_forward, u₀, ts[1])
u_test, DM_test = Zygote.jacobian(model_forward, u₀, ts[1])




# Perform Assimilation Loop

const u_now::Vector{Float64} = zeros(size(uₐ[:,1]))
const DM::Matrix{Float64} = DM_test[1]

@info "Starting assimilation..."



@showprogress for k ∈ 1:length(ts)-1  # because we always update the *next* value
    # k=2

    # collect current model estimate
    u_now .= uₐ[:,k]  # should preallocate this

    # --------------------------
    # Forecast Step
    # --------------------------

    # run model forward one Δt
    u_next, DM_tup = Zygote.withjacobian(model_forward, u_now, ts[k])  # <-- can't make this mutating

    DM .= DM_tup[1]  # result is a 1-element tuple, so we index it

    # collect observations
    local is_meas_not_nan = get_idxs_not_nans(W[:,k+1])
    local idx_meas_nonan = idx_meas[is_meas_not_nan]
    local u_h = Obs(u_next, idx_meas_nonan)


    # update loop for the mode covariance matrix, Q
    if k > 1
        Q .= 0.0
        for j ∈ axes(Q,1)
            # define error growth rate
            eg = abs(u_next[j]-u_now[j])/max(1.0, u_next[j])  # note: we should actually do the sens. analysis here

            # clip the growth rate to 1%-50% range
            eg = max(min(eg, 0.5), 0.01)


            # if we have a measurement, further manipulate the growth rate
            if j ∈ idx_meas_nonan
                u_h_j = u_h[idx_meas_nonan .== j][1]
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

    DH = JObs(u_h, u_next, idx_meas_nonan)
    R = Rmat_nonan(k+1, is_meas_not_nan, meas_ϵ; fudge_fac=fudge_fac)
    denom = DH*P*DH' + R

    #Kalman = P*DH'/denom
    Kalman = P*DH'* inv(denom)  # <-- this seems slightly faster for now

    u_next .= u_next .+ Kalman*(W[is_meas_not_nan, k+1] - u_h)

    P .= (I(length(u_next)) - Kalman*DH)*P


    # do some clipping to make sure we stay reasonable
    for j ∈ axes(P,1)
        P[j,j] = min(max(0, u_next[j]^2), P[j,j])
        P[j,j] = max(0.05*u_next[j]^2, P[j,j])
    end


    # filter negative values to zero
    u_next[u_next .≤ 0.0] .= 0.0

    # update the analysis vector
    update_uₐ!(uₐ, u_next, k+1)
    P_diag[:,k+1] .= [P[i,i] for i ∈ axes(P,1)]
end



# --------------------------------------------------------------------------------------------------------------------------
# 12. Post Process
# --------------------------------------------------------------------------------------------------------------------------


# combine output w/ uncertainty from diagonal
uₐ_nd = uₐ .± sqrt.(P_diag)

# convert final output into mixing ratios
M = df_params.M .± (fudge_fac .* df_params_ϵ.M)

uₐ_mr = to_mixing_ratio(uₐ_nd, M)

# chop off values and uncertainties for easier plotting
ua_mr_vals = Measurements.value.(uₐ_mr)
ua_mr_ϵ = Measurements.uncertainty.(uₐ_mr)


# save to output file
df_ekf = DataFrame()
df_ekf_ϵ = DataFrame()

@showprogress for i ∈ axes(ua_mr_vals, 1)
    df_ekf[!, df_species[i, "MCM Name"]] = ua_mr_vals[i,:]
    df_ekf_ϵ[!, df_species[i, "MCM Name"]] = ua_mr_ϵ[i,:]
end

df_ekf[!, :times] = ts
df_ekf_ϵ[!, :times] = ts

CSV.write("models/$model_name/EKF/ekf_output.csv", df_ekf)
CSV.write("models/$model_name/EKF/ekf_ϵ_output.csv", df_ekf_ϵ)


# combine measurements with uncertainties
W_mr = W .± (fudge_fac .* meas_ϵ)
W_mr = to_mixing_ratio(W_mr, M)

W_mr_val = Measurements.value.(W_mr)
W_mr_ϵ = Measurements.uncertainty.(W_mr)

# save measurements to csv files for final output
size(W_mr_val)
idx_meas

df_w = DataFrame()
df_w_ϵ = DataFrame()
@showprogress for i ∈ axes(W_mr_val, 1)
    df_w[!, df_species[idx_meas[i], "MCM Name"]] = W_mr_val[i,:]
    df_w_ϵ[!, df_species[idx_meas[i], "MCM Name"]] = W_mr_ϵ[i,:]
end

CSV.write("models/$model_name/EKF/ekf_measurements.csv", df_w)
CSV.write("models/$model_name/EKF/ekf_measurements_ϵ.csv", df_w_ϵ)



# --------------------------------------------------------------------------------------------------------------------------
# 13. Plots
# --------------------------------------------------------------------------------------------------------------------------


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

# Collision reaction w/ RO2 (i.e. all termolecular reactions) will look the same as above since M is included already inside of our computed k.

derivatives
derivatives_ro2

# now we need to combine

ua_nd =  Measurements.value.(uₐ_nd)
size(ua_nd)

τs = copy(ua_nd)  # preallocate matrix to hold values
ℓ_mat = zeros(size(ua_nd))  # loss rate
#ℓ = 1.0

@showprogress for d ∈ 1:length(derivatives)
    derivative = derivatives[d]
    if derivative.prefac < 0.0 # i.e. if it's negative so that we have a reactant not product
        for idx_t ∈ axes(K_matrix,1)
            ℓ = K_matrix[idx_t,derivative.idx_k]
            for i ∈ derivative.idxs_in
                ℓ  *= ua_nd[i, idx_t]
            end

            ℓ_mat[derivative.idx_du, idx_t] += ℓ
        end
    end
end

@showprogress for d ∈ 1:length(derivatives_ro2)
    derivative = derivatives_ro2[d]
    if derivative.prefac < 0.0 # i.e. if it's negative so that we have a reactant not product
        for idx_t ∈ axes(K_matrix,1)
            ℓ = K_matrix[idx_t,derivative.idx_k]
            for i ∈ derivative.idxs_in
                ℓ  *= ua_nd[i, idx_t]
            end

            ℓ_mat[derivative.idx_du, idx_t] += ℓ
        end
    end
end


for j ∈ axes(τs, 2), i ∈ axes(τs,1)
    if isinf(τs[i,j]/ℓ_mat[i,j] ) || isnan(τs[i,j]/ℓ_mat[i,j] )
        println("idx: ", (i,j), "\tu:\t", τs[i,j], "\tℓ:\t",ℓ_mat[i,j], "\tτ:\t", τs[i,j]/ℓ_mat[i,j] )
    end

    τs[i,j] = τs[i,j] / ℓ_mat[i,j]
end



# generate lifetime dataframes for output
df_τs = DataFrame()
@showprogress for i ∈ axes(τs, 1)
    df_τs[!, df_species[i, "MCM Name"]] = τs[i,:]
end

df_τs.t = df_params.t

CSV.write("models/$model_name/EKF/lifetimes.csv", df_τs)


df_τs_means = DataFrame()
@showprogress for i ∈ axes(τs, 1)
    df_τs_means[!, df_species[i, "MCM Name"]] = [mean(τs[i,:])]
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


CSV.write("models/$model_name/EKF/mean_lifetimes.csv", df_τs_means_sorted[7:end,:])

