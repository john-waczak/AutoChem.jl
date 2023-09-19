using AutoChem

using DelimitedFiles, CSV, DataFrames
using JSON

using LinearAlgebra

using DifferentialEquations
using Sundials, OrdinaryDiffEq
using Zygote, ForwardDiff
using SciMLSensitivity

using Optimization
using OptimizationOptimJL
using OptimizationFlux

using ParameterHandling, EarlyStopping

using StableRNGs

using BenchmarkTools

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

df_u₀ = CSV.read(joinpath(model_path, model_name, "mechanism", "u0.csv"), DataFrame);

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

const fudge_fac::Float64 = 0.5

const tmin::Float64 = minimum(ts)
const tmax::Float64 = 0.0 # maximum(ts)
const abstol::Float64 = 1e-3
const reltol::Float64 = 1e-3
const tspan = (tmin, tmax)

const ϵ::Float64 = 0.5
const ϵ_min::Float64 = 1e-12


# load rhs and jacobian functions
@info "Loading rhs and jac functions"
include("models/$model_name/mechanism/rhs.jl")
include("models/$model_name/mechanism/jacobian.jl")

# define the ODE function
@info "Defining ODE function"
fun = ODEFunction(rhs!; jac=jac!)

#ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀ , tspan)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀ , (tmin, tmin+Δt_step)

@info "Trying a solve with default u₀"
# sol = solve(ode_prob, QNDF(); saveat=Δt_step, reltol=reltol, abstol=abstol)

sol = solve(ode_prob; saveat=Δt_step)


# set up observation observation operator and it's jacobian
@info "Testing observation operator"
Obs(u₀, idx_meas, idx_pos, idx_neg)
JObs(u₀, idx_meas, idx_pos, idx_neg)

@info "Testing Observation Covariance Matrix and its Inverse"
Rmat(1, meas_ϵ; fudge_fac=fudge_fac)
Rinv(1, meas_ϵ; fudge_fac=fudge_fac)


# Establish analysis vector with small initial offset
@info "initializing analysis vector"
current_u0a = (;u0 = [positive(u+100.0) for u ∈ u₀ ])
u0a, unflatten = ParameterHandling.value_flatten(current_u0a)
# NOTE: `unflatten` both reconstructs the NamedTuple with our paramters and applies inverse transform to Positive

# Set up loss function for 4d-var
@info "Setting up loss function..."
const u0b::Vector{Float64} = copy(u0a) # i.e. "background guess"
const B::Matrix{Float64} = diagm((ϵ .* (u₀)) .^2  .+ ϵ_min^2)
const Binv::Matrix{Float64} = inv(B)

const use_background_cov::Bool = false


sensealg_dict = Dict(
    :QuadratureAdjoint => QuadratureAdjoint(),
    :BacksolveAdjoint => BacksolveAdjoint(),
    :InterpolatingAdjoint => InterpolatingAdjoint(),
    :ZygoteAdjoint => ZygoteAdjoint()
)

sensealg = sensealg_dict[:QuadratureAdjoint]


#function loss(log_u0a)
function loss(u0a)
    u0a_now = unflatten(u0a)


    # remake problem using current value
    _prob = remake(ode_prob; u0=u0a_now.u0)

    # integrate the model forward
    sol = solve(
        _prob,
        TRBDF2();
        saveat=Δt_step,
        reltol=reltol,
        abstol=abstol,
        sensealg=sensealg,
        verbose=false
    )

    # l = sum((sol))
    # compute loss
    l = 0.0

    for j ∈ 1:length(sol.t)
        #idx_t = get_time_index(sol.t[j], Δt_step, tmin)
        # l += 0.5 * ((W[:,idx_t] .- Obs(sol[j], idx_meas, idx_pos, idx_neg))' * Rinv(idx_t, meas_ϵ; fudge_fac=fudge_fac) * (W[:,idx_t] .- Obs(sol[j], idx_meas, idx_pos, idx_neg)))[1]
        l += 0.5 * ((W[:,j] .- Obs(sol[j], idx_meas, idx_pos, idx_neg))' * Rinv(j, meas_ϵ; fudge_fac=fudge_fac) * (W[:,j] .- Obs(sol[j], idx_meas, idx_pos, idx_neg)))[1]
    end

    # optionally, add additional loss term quantifying our belief in the inital condition vector
    if use_background_cov
        l += 0.5*(u0a-u0b)'*Binv*(u0a-u0b)
    end

    return l
end

@info "Testing loss function..."
loss(u0a)


Zygote.gradient(loss, u0a)


# Set up the optimization problem(s)
@info "Setting up optimization problem(s)"
# set up stoppers for each callback function
stopper1 = EarlyStopper(Disjunction(Warmup(Patience(n=5);n=10), Warmup(NumberSinceBest(n=10);n=10), TimeLimit(t=24.0)))

losses1::Vector{Float64} = []
losses2::Vector{Float64} = []

df_out = DataFrame(unflatten(u0a))
if ! ispath(joinpath(model_path, model_name, "4d-var"))
    mkpath(joinpath(model_path, model_name, "4d-var"))
end

CSV.write(joinpath(model_path, model_name, "4d-var", "u0.csv"), df_out)

function callback(u0a, lossval)
    println("current loss: ", lossval)
    df_out.u0 = unflatten(u0a).u0
    CSV.write(joinpath(model_path, model_name, "4d-var", "u0.csv"), df_out)
    push!(losses1, lossval)
    done!(stopper1, lossval)  # return false unless we've matched a stopping criterion
end

Zygote.gradient(loss, u0a)

@info "Trying out gradient of loss function"
res = Zygote.gradient(loss, u0a)

sum(res[1])
    # b = @benchmark Zygote.gradient(loss, u0a)

