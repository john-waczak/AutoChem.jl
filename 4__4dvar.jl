ENV["GKSwstype"] = 100

@info "Setting up Julia Environment..."
using Pkg
Pkg.activate(".")
Pkg.instantiate()
@info "\t...finished"

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
using ArgParse

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




# set up function with ArgParse macros
# to parse the command line flags
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--data_basepath"
            help = "Path to data files to be used for testing"
            arg_type = String
            default = "data/intertek-emergency-testing"
        "--collection_id"
            help = "Name of collection to analyze"
            arg_type = String
            default = "empty"
        "--unc_ext"
            help = "Extension for uncertainty files."
            arg_type = String
            default = "_std"
        "--model_name"
            help = "Name for the resulting model used in output paths"
            arg_type = String
            #default = "autochem-w-ions"
            default = "methane"
        "--time_step"
            help = "The time step used during integration of mechanism (in minutes)."
            arg_type = Float64
            default = 15.0
        "--use_background_cov"
            help = "Whether or not to use background covariance matrix in loss"
            action = :store_true
        "--fudge_fac"
            help = "A fudge factor for manipulating scale of measurement uncertainties"
            arg_type = Float64
            # default = 0.5
            default = 1.0
        "--epsilon"
            help = "Estimated background uncertainty for diagonal of B matrix, i.e. uncertainty in initial condition"
            arg_type = Float64
            default = 0.5
        "--sensealg"
            help = "Method for computing sensitivities of loss function w.r.t. initial condition vector"
            arg_type = Symbol
            default = :ForwardDiffSensitivity
    end


    parsed_args = parse_args(ARGS, s; as_symbols=true)

    return parsed_args
end




@info "Parsing command line flags..."

parsed_args = parse_commandline()

data_basepath = parsed_args[:data_basepath]
model_name = parsed_args[:model_name]
collection_id = parsed_args[:collection_id]
unc_ext = parsed_args[:unc_ext]
model_path= "models/"

@info "Setting up file paths..."
outpath = joinpath(model_path, model_name, "runs", collection_id)
@assert ispath(outpath)

docs_path = joinpath(model_path, model_name, "docs")
@assert ispath(docs_path)

if !ispath(joinpath(outpath, "figures"))
    joinpath(model_path, model_name, "figures")
end




# fetch species list
@info "Fetching data"

df_species = CSV.read(joinpath(outpath, "mechanism", "species.csv"), DataFrame);
df_params= CSV.read(joinpath(outpath, "mechanism", "state_parameters.csv"), DataFrame);
df_params_ϵ = CSV.read(joinpath(outpath, "mechanism", "state_parameters_ϵ.csv"), DataFrame);
df_nd = CSV.read(joinpath(outpath, "mechanism", "number_densities.csv"), DataFrame);
df_nd_ϵ = CSV.read(joinpath(outpath, "mechanism", "number_densities_ϵ.csv"), DataFrame);
df_ions = CSV.read(joinpath(outpath, "mechanism", "ions.csv"), DataFrame)
df_ions_ϵ = CSV.read(joinpath(outpath, "mechanism", "ions_ϵ.csv"), DataFrame)


df_u₀ = CSV.read(joinpath(outpath, "mechanism", "u0.csv"), DataFrame);

# set initial conditions
@info "Getting initial condition vector"
u₀ = df_u₀.u0[:]


# generate integrated vs non-integrated indices
@info "Generating species indices"
const idx_noint = findall(df_species.is_integrated .== 2)
const n_integrated = sum(df_species.is_integrated .== 1)
const idx_pos = get_positive_indices(df_species)
const idx_neg = get_negative_indices(df_species)
const Δt_step = parsed_args[:time_step]


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

# decrease uncertainty of ion counts
meas_ϵ[end-1:end,:] .= 0.1 .* meas_ϵ[end-1:end,:]

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
const fudge_fac::Float64 = parsed_args[:fudge_fac]
#const fudge_fac::Float64 = 1.0

const tmin::Float64 = minimum(ts)
const tmax::Float64 = 0.0 # maximum(ts)
const abstol::Float64 = 1e-3
const reltol::Float64 = 1e-3
const tspan = (tmin, tmax)

const ϵ::Float64 = parsed_args[:epsilon]
const ϵ_min::Float64 = 1e-12


# load rhs and jacobian functions
@info "Loading rhs and jac functions"
include(joinpath(outpath, "mechanism", "rhs.jl"))
include(joinpath(outpath, "mechanism", "jacobian.jl"))

# define the ODE function
@info "Defining ODE function"
fun = ODEFunction(rhs!; jac=jac!)

# ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀ , tspan)

#idx_0 = findfirst(ts .== 0) + 1
idx_0 = 1
#n_steps_to_use = 4  # e.g. 1 hour worth
n_steps_to_use = length(ts)  # e.g. 1 hour worth


# set inital values for ions to be equal
idx_pos
idx_neg
W[end-1, idx_0]
W[end, idx_0]

u₀[idx_pos] .= (W[end-1, idx_0] / length(idx_pos))
u₀[idx_neg] .= (W[end, idx_0] / length(idx_neg))

# add all negative ion counts to O2-
# u₀[82] = W[end, idx_0]



#ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀ .+ 100 , (ts[idx_0], ts[idx_0+n_steps_to_use]))
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀ .+ 100 , (ts[idx_0], ts[end]))

@info "Trying a solve with default u₀"
sol = solve(ode_prob; alg_hints=[:stiff], saveat=Δt_step, reltol=reltol, abstol=abstol)



# set up observation observation operator and it's jacobian
@info "Testing observation operator"
Obs(u₀, idx_meas, idx_pos, idx_neg)
JObs(u₀, idx_meas, idx_pos, idx_neg)

@info "Testing Observation Covariance Matrix and its Inverse"
Rmat(1, meas_ϵ; fudge_fac=fudge_fac)
Rinv(1, meas_ϵ; fudge_fac=fudge_fac)


# Establish analysis vector with small initial offset
@info "initializing analysis vector"
current_u0a = (;u0 = [positive(u+1.0) for u ∈ u₀ ])
u0a, unflatten = ParameterHandling.value_flatten(current_u0a)
# NOTE: `unflatten` both reconstructs the NamedTuple with our paramters and applies inverse transform to Positive

# Set up loss function for 4d-var
@info "Setting up loss function..."
const u0b::Vector{Float64} = copy(u0a) # i.e. "background guess"
const B::Matrix{Float64} = diagm((ϵ .* (u₀)) .^2  .+ ϵ_min^2)
const Binv::Matrix{Float64} = inv(B)

const use_background_cov::Bool = parsed_args[:use_background_cov]


sensealg_dict = Dict(
    :QuadratureAdjoint => QuadratureAdjoint(),
    :BacksolveAdjoint => BacksolveAdjoint(),
    :InterpolatingAdjoint => InterpolatingAdjoint(),
    :ZygoteAdjoint => ZygoteAdjoint(),
    :ForwardDiffSensitivity => ForwardDiffSensitivity(),
    :ForawrdSensitivity => ForwardSensitivity(),
)

sensealg = sensealg_dict[parsed_args[:sensealg]]

# Both of these worked:
# sensealg = sensealg_dict[:ForwardDiffSensitivity]
# sensealg = sensealg_dict[:BacksolveAdjoint]


#function loss(log_u0a)
function loss(u0a)
    u0a_now = unflatten(u0a)


    # remake problem using current value
    _prob = remake(ode_prob; u0=u0a_now.u0)

    # integrate the model forward
    sol = solve(
        _prob;
        alg_hints=[:stiff],
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
        idx_t = get_time_index(sol.t[j], Δt_step, tmin)
        l += 0.5 * ((W[:,idx_t] .- Obs(sol[j], idx_meas, idx_pos, idx_neg))' * Rinv(idx_t, meas_ϵ; fudge_fac=fudge_fac) * (W[:,idx_t] .- Obs(sol[j], idx_meas, idx_pos, idx_neg)))[1]
        #l += 0.5 * ((W[:,j] .- Obs(sol[j], idx_meas, idx_pos, idx_neg))' * Rinv(j, meas_ϵ; fudge_fac=fudge_fac) * (W[:,j] .- Obs(sol[j], idx_meas, idx_pos, idx_neg)))[1]
    end

    # optionally, add additional loss term quantifying our belief in the inital condition vector
    if use_background_cov
        l += 0.5*(u0a-u0b)'*Binv*(u0a-u0b)
    end

    return l
end

@info "Testing loss function..."
loss(u0a)



@info "Trying out gradient of loss function"
Zygote.gradient(loss, u0a)

@benchmark Zygote.gradient(loss, u0a)  # :BacksolveAdjoint ~ 108 s

# @benchmark Zygote.gradient(loss, u0a)  # :ForwardDiffSensitivity ~6.826 seconds



# Set up the optimization problem(s)
@info "Setting up optimization problem(s)"
# set up stoppers for each callback function
#stopper1 = EarlyStopper(Disjunction(Warmup(Patience(n=5);n=10), Warmup(NumberSinceBest(n=10);n=10), TimeLimit(t=24.0)))
stopper1 = EarlyStopper(Disjunction(Warmup(Patience(n=5);n=10), Warmup(NumberSinceBest(n=10);n=10), TimeLimit(t=1.0)))

losses1::Vector{Float64} = []
losses2::Vector{Float64} = []

df_out = DataFrame(unflatten(u0a))
if ! ispath(joinpath(outpath, "4d-var"))
    mkpath(joinpath(outpath, "4d-var"))
end

CSV.write(joinpath(outpath, "4d-var", "u0.csv"), df_out)

function callback(u0a, lossval)
    println("current loss: ", lossval)
    df_out.u0 = unflatten(u0a).u0
    CSV.write(joinpath(outpath, "4d-var", "u0.csv"), df_out)
    push!(losses1, lossval)
    done!(stopper1, lossval)  # return false unless we've matched a stopping criterion
end


# define the optimization function and declare we are using Zygote for AutoGrad
optf = OptimizationFunction((x,p)->loss(x), Optimization.AutoZygote())
opt_prob = Optimization.OptimizationProblem(optf, u0a)

@info "Starting first round of optimization:"
# solve first with ADAM which is fast but can get stuck in local minimum

method1 = ADAM()
method2 = BFGS(initial_stepnorm=0.01)

# #method2=LBFGS()  # <-- took bad steps
# if want_restart
#     method1 = ADAM(0.005)
# end

opt_sol = solve(
    opt_prob,
    method1,
    callback=callback,
    maxiters=2e5,
)




@info "Starting second round of optimization"
u0a = opt_sol.u
prob2 = remake(opt_prob, u0=u0a)

stopper2 = EarlyStopper(Disjunction(Warmup(Patience(n=5);n=3), Warmup(NumberSinceBest(n=30);n=3), TimeLimit(t=24.0)))

function callback2(u0a, lossval)
    println("current loss: ", lossval)
    df_out.u0 = unflatten(u0a).u0
    CSV.write(joinpath(outpath, "4d-var", "u0.csv"), df_out)
    push!(losses2, lossval)
    done!(stopper2, lossval)  # return false unless we've matched a stopping criterion
end


opt_sol = solve(prob2,
                method2,
                callback=callback2,
                allow_f_increases=false,
                )

# see example here: https://docs.sciml.ai/DiffEqFlux/stable/examples/neural_ode/

u0a_final = unflatten(opt_sol.u).u0


fig = Figure();
ax = Axis(fig[1,1], xlabel="Iteration", ylabel="Loss", title="4D-Var Training");

iterations1 = 1:length(losses1)
iterations2 = (iterations1[end]):(iterations1[end]+length(losses2)-1)

l1 = lines!(ax, iterations1, losses1, linewidth=3)
l2 = lines!(ax, iterations2, losses2, linewidth=3)

labels = ["ADAM", "BFGS"]

axislegend(ax, [l1, l2], labels, position=:rt)

fig

save(joinpath(outpath, "4d-var", "losses.png"), fig)
save(joinpath(outpath, "4d-var", "losses.svg"), fig)
save(joinpath(outpath, "4d-var", "losses.eps"), fig)
save(joinpath(outpath, "4d-var", "losses.pdf"), fig)


# save final values
length(u0a_final)

fig = Figure();
ax = Axis(fig[1,1], ylabel="Percent change from u₀", xticks = (1:length(idx_meas), String.(df_species.varname[idx_meas])), xticklabelrotation=π/2, xticklabelsize=20);
barplot!(
    ax,
    1:length(idx_meas),
    (u0a_final[idx_meas] .- u₀[idx_meas])./(u₀[idx_meas]) .* 100
)
fig

save(joinpath(outpath, "4d-var", "u0-change.png"), fig)
save(joinpath(outpath, "4d-var", "u0-change.svg"), fig)
save(joinpath(outpath, "4d-var", "u0-change.eps"), fig)
save(joinpath(outpath, "4d-var", "u0-change.pdf"), fig)


df_out.u0 = u0a_final
CSV.write(joinpath(outpath, "4d-var", "u0_final.csv"), df_out)

df_species[end-37:end,:]
u0a_final[end-37:end]


# remake problem using current value
# _prob = remake(ode_prob; u0=u0a_final, tspan=(ts[idx_0], ts[end]))

# # integrate the model forward
# sol = solve(
#     _prob;
#     alg_hints=[:stiff],
#     saveat=Δt_step,
#     reltol=reltol,
#     abstol=abstol,
# )

# println(minimum(sol[:,:]))


# u0a_final

# u0a_final[end-40:end]
# df_species[end-40:end,:]

# df_species[end-15:end,:]

# u0a_final[end-15]

# df_species
# lines(ts[idx_0:end], sol[74, :])

