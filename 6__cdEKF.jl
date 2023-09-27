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

const ts = df_params.t

#const fudge_fac::Float64 = 0.1
#const fudge_fac::Float64 = 0.5
const fudge_fac::Float64 = 0.5

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

# set up observation observation operator and it's jacobian
@info "Testing observation operator"
Obs(u₀, idx_meas, idx_pos, idx_neg)
JObs(u₀, idx_meas, idx_pos, idx_neg)

@info "Testing Observation Covariance Matrix and its Inverse"
Rmat(1, meas_ϵ; fudge_fac=fudge_fac)
Rinv(1, meas_ϵ; fudge_fac=fudge_fac)


@info "Pre-allocating matrices for EKF"
const P::Matrix{Float64} = zeros(nrow(df_species), nrow(df_species))
const P_diag::Matrix{Float64} = zeros(nrow(df_species), length(ts)) # i.e. P_diag[i] == P[i,i]
const Q::Matrix{Float64} = zeros(size(P))
const uₐ::Matrix{Float64} = zeros(length(u₀), length(ts))


# Initialize Covariance Matrices
@info "Initializing covariance matrices"

for i ∈ 1:length(u₀)
    P[i,i] = (ϵ * u₀[i])^2 + ϵ_min^2
    P_diag[i,1] = P[i,i]
end

# initially, set Q to match P
Q .= P


# Set up sensitivity Algorithm
@info "Choosing sensitivity algorithm and solver"
sensealg_dict = Dict(
    :QuadratureAdjoint => QuadratureAdjoint(),
    :BacksolveAdjoint => BacksolveAdjoint(),
    :InterpolatingAdjoint => InterpolatingAdjoint(),
    :ZygoteAdjoint => ZygoteAdjoint(),
    :ForwardDiffSensitivity => ForwardDiffSensitivity()
)

sensealg = sensealg_dict[:BacksolveAdjoint]
solver = TRBDF2()






# define the ODE function
@info "Defining ODE function"
fun = ODEFunction(rhs!; jac=jac!)

# we need to figure out if starting at a different time will cause problems...
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀ , (ts[1], ts[2]))
# ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun,df_u0_orig.u0, tspan)

@info "Trying a solve with default u₀"
sol = solve(ode_prob, TRBDF2(); reltol=reltol, abstol=abstol)
# sol = solve(ode_prob; saveat=Δt_step, reltol=reltol, abstol=abstol)
# sol = solve(ode_prob, CVODE_BDF(); saveat=Δt_step, reltol=reltol, abstol=abstol)
# sol = solve(ode_prob, QNDF(); saveat=Δt_step, reltol=reltol, abstol=abstol)


const F = zeros(length(u₀), length(u₀))
jac!(F, sol(ts[1]), nothing, ts[1])


const Mult_temp = zeros(length(u₀), length(u₀))

function P_rhs!(dP, P, p, t)
    sol = p
    # update ODE Jacobian
    jac!(F, sol(t), nothing, t)

    # initialize to zero
    dP .= 0.0

    # perform matrix multiplications
    mul!(Mult_temp, P, F')
    mul!(dP, Jac, Mult_temp)

    dP .+ Q
end

dP = zeros(length(u₀), length(u₀))
P_rhs!(dP, P, sol, ts[1])



cov_prob = ODEProblem{true, SciMLBase.FullSpecialize}(P_rhs!, P, (ts[1], ts[2]), sol)
sol = solve(cov_prob; reltol=reltol, abstol=abstol)



function updateQ!(Q, u_now, u_next, u_h, idx_meas)
    Q .= 0.0

    # only mess with diagonals
    for j ∈ axes(Q,1)
        # define error growth rate
        eg = abs(u_next[j]-u_now[j])/max(1.0, u_next[j])

        # clip to 1%-50% range
        eg = clamp(eg, 0.01, 0.5)

        # if we have a measurement for the species, further manipulate the growth rate
        if j ∈ idx_meas
            u_h_j = u_h[findfirst(idx_meas .== j)]
            diff = abs(u_next[j] - u_h_j)
            ratio = diff/u_h_j

            if ratio > 1
                eg = max(eg, 0.5)
            elseif ratio > 0.5
                eg = max(eg. 0.2)
            elseif ratio > 0.2
                eg = max(eg, 0.1)
            else
                continue
            end
        end
        Q[j,j] = eg * u_next[j]^2 + ϵ_min^2
    end
end




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
model_forward(u₀, ts[1])
@benchmark ForwardDiff.jacobian!(res, u->model_forward(u,ts[1]), u₀)

# u_test, DM_test = Zygote.withjacobian(model_forward, u₀)  # way too long

u_test = res.value
DM_test = res.derivs[1]

# Perform Assimilation Loop

const u_now::Vector{Float64} = zeros(size(uₐ[:,1]))
const DM::Matrix{Float64} = DM_test

@info "Starting assimilation..."

any(isnan.(W))

@showprogress for k ∈ 1:length(ts)-1  # because we always update the *next* value
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
    # local is_meas_not_nan = get_idxs_not_nans(W[:,k+1])
    # local idx_meas_nonan = idx_meas[is_meas_not_nan]
    # local u_h = Obs(u_next, idx_meas, idx_pos, idx_neg)
    u_h = Obs(u_next, idx_meas, idx_pos, idx_neg)


    # update loop for the mode covariance matrix, Q
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
    P .= DM*P*DM' # + Q

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


df_species
lines(ts, uₐ[82,:])


# --------------------------------------------------------------------------------------------------------------------------
# 12. Post Process
# --------------------------------------------------------------------------------------------------------------------------


using Measurements

# combine output w/ uncertainty from diagonal
uₐ_nd = uₐ .± sqrt.(P_diag)

fig, ax, b = band(ts, uₐ[82,:] .- sqrt.(P_diag[82,:]), uₐ[82,:] .+ sqrt.(P_diag[82,:]); color=(mints_colors[1],0.2))
l = lines!(ax, ts, uₐ[82,:], lw=3)
fig

fig, ax, b = band(ts, uₐ[idx_meas[4],:] .- sqrt.(P_diag[idx_meas[4],:]), uₐ[idx_meas[4],:] .+ sqrt.(P_diag[idx_meas[4],:]); color=(mints_colors[1],0.2))
l = lines!(ax, ts, uₐ[idx_meas[4],:], lw=3)
errorbars!(ax, ts, W[4,:], meas_ϵ[4,:])
scatter!(ax, ts, W[4,:])
fig

df_species[idx_meas,:]


