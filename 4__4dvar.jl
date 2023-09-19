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
# model_name = "autochem-w-ions"
model_name = "qroc-methane-intel"
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
U_noint[1,:] .= Ar.(temperatures, pressures)
U_noint[2,:] .= O2.(temperatures, pressures)
U_noint[3,:] .= N2.(temperatures, pressures)
U_noint[4,:] .= ones(Float64, nrow(df_params))
U_noint[5,:] .= M.(temperatures, pressures)



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

# extrema(K_bimol)
# extrema(K_trimol)  # <-- is this the issue ???
# extrema(K_photo)


# heatmap(K_bimol)
# argmax(K_trimol)
# heatmap(K_trimol)  # <-- troubles seem to be in 50 - 65 rnage


# K_bi_mean = mean(K_bimol, dims=2)[:,1]
# K_tri_mean = mean(K_trimol, dims=2)[:,1]
# K_photo_mean = mean(K_photo, dims=2)[:,1]



# barplot(1:length(K_tri_mean), K_tri_mean)
# idx_rxns_bad = findall(K_tri_mean .> 0.5)

# db_trimol[idx_rxns_bad]


# idx_41_bi = []
# for i ∈ 1:length(db_bimol)
#     rxn = db_bimol[i]
#     if 41 ∈ rxn.reactants
#         push!(idx_41_bi, i)
#     end
# end

# idx_41_tri = []
# for i ∈ 1:length(db_trimol)
#     rxn = db_trimol[i]
#     if 41 ∈ rxn.reactants
#         push!(idx_41_tri, i)
#     end
# end

# idx_41_photo = []
# for i ∈ 1:length(db_photo)
#     rxn = db_photo[i]
#     if 41 ∈ rxn.reactants
#         push!(idx_41_photo, i)
#     end
# end


# db_bimol[idx_41_bi]
# db_trimol[idx_41_tri]
# db_photo[idx_41_photo]

# K_tri_mean[idx_41_tri]*get_concentration(94,1,u₀, U_noint, n_integrated)

# K_tri_mean[idx_41_tri]


# K_tri_mean[idx_41_tri]  #
# K_photo_mean[idx_41_photo]



# get_concentration(12, 1, u₀, U_noint, n_integrated)
# get_concentration(33, 1, u₀, U_noint, n_integrated)
# get_concentration(94, 1, u₀, U_noint, n_integrated)
# get_concentration(93, 1, u₀, U_noint, n_integrated)


# K_bi_mean[idx_41_bi]

# db_trimol[51] # <-- is this the offending reaction ?

# idx_bad = [41, 81, 82, 85, 87, 88]

# idx_good = [i for i ∈ 1:length(u₀) if !(i∈idx_bad)]

# df_species[idx_bad,:]


# db_
# test_u0 = copy(u₀)
# test_u0[idx_good] .+= 1
# #test_u0[41] += 1

# du = zeros(length(u₀))
# rhs!(du, test_u0, nothing, tmin)

# minimum(du)

# # du[41]

# # maximum(du)
# # minimum(du)


#ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀ , tspan)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, test_u0 , tspan)

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
