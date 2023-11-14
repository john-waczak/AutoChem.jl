ENV["GKSwstype"] = 100

@info "Setting up Julia Environment..."
using Pkg
Pkg.activate(".")
Pkg.instantiate()
@info "\t...finished"


using AutoChem
using DelimitedFiles, CSV, DataFrames
using JSON
using DifferentialEquations
using Sundials, OrdinaryDiffEq
using BenchmarkTools, ArgParse

using CairoMakie
using MintsMakieRecipes
set_theme!(mints_theme)


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
            #default = "high_primed"
            default = "empty"
        "--unc_ext"
            help = "Extension for uncertainty files."
            arg_type = String
            default = "_std"
        "--model_name"
            help = "Name for the resulting model used in output paths"
            arg_type = String
            default = "autochem-w-ions"
        "--time_step"
            help = "The time step used during integration of mechanism (in minutes)."
            arg_type = Float64
            default = 15.0
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
const Δt_step = parsed_args[:time_step]  # time step in minutes
model_path= "models/"

@info "Setting up file paths..."

outpath = joinpath(model_path, model_name, "runs", collection_id)
@assert ispath(outpath)

docs_path = joinpath(model_path, model_name, "docs")
@assert ispath(docs_path)

if !ispath(joinpath(outpath, "figures"))
    joinpath(model_path, model_name, "figures")
end




# -----
# 1. fetch species and create indices
# -----
@info "Fetching species list"
df_species = CSV.read(joinpath(outpath, "mechanism", "species.csv"), DataFrame)

@info "Generating ion indices"
const idxs_positive = get_positive_indices(df_species)
const idxs_negative = get_negative_indices(df_species)



# -----
# 2. generate lookup tables for measurements
# -----

# print the non-integrated species
df_species[df_species.is_integrated .== 2, :]
# Ar is 0.94% of air: https://www.rsc.org/periodic-table/element/18/argon#:~:text=Argon%20makes%20up%200.94%25%20of,third%20most%20abundant%20atmospheric%20gas.


replacement_dict = Dict(
    "CH3OH" => "MeOH"
)

generate_densities(
    joinpath(data_basepath, "number_densities", collection_id, "number_densities.csv"),
    joinpath(data_basepath, "number_densities", collection_id, "number_densities"*unc_ext*".csv"),
    outpath;
    replacement_dict=replacement_dict
)


# now let's copy over the ion concentrations
@info "Create Ion Concentration tables"
df_ion1 = CSV.read(joinpath(data_basepath, "pseudo", collection_id, "ion1.csv"), DataFrame)
df_ion2 = CSV.read(joinpath(data_basepath, "pseudo", collection_id, "ion2.csv"), DataFrame)

# now we need to decide which is positive and which is negative
sign_1 = "positive_ions"
if df_ion1.channel[1] < 0
    sign_1 = "negative_ions"
end

sign_2 = "positive_ions"
if df_ion2.channel[1] < 0
    sign_2 = "negative_ions"
end

df_ions = DataFrame(sign_1=>df_ion1.count, sign_2=>df_ion2.count)
df_ions_ϵ = DataFrame(sign_1=>df_ion1[:, "count"*unc_ext], sign_2=>df_ion2[:, "count"*unc_ext])

# re-sort so positive is always first regardless of which was in pos/neg mode
df_ions = df_ions[:, ["positive_ions", "negative_ions"]]
df_ions_ϵ = df_ions_ϵ[:, ["positive_ions", "negative_ions"]]

# save the results
CSV.write(joinpath(outpath, "mechanism", "ions.csv"), df_ions)
CSV.write(joinpath(outpath, "mechanism", "ions_ϵ.csv"), df_ions_ϵ)


df_params = CSV.read(joinpath(outpath, "mechanism", "state_parameters.csv"), DataFrame)
df_number_densities = CSV.read(joinpath(outpath, "mechanism", "number_densities.csv"), DataFrame)


# -----
# 5. generate sane initial conditions
# -----
@info "Generate best guess initial condition"
init_path = "./assets/initial_concentrations/full.txt"
@assert isfile(init_path) "Cant find ./assets/initial_concentrations/full.txt"

init_dict = generate_init_dict(init_path, df_params.M[1])


for (key,val) ∈ init_dict
    println(key,":\t",val)
end

# this should be put in a config file
measurements_to_ignore = [:C2H6, :SO2]  # skip any with nans or negative values



const idx_noint = findall(df_species.is_integrated .== 2)
const n_integrated = sum(df_species.is_integrated .== 1)

n_integrated

u₀ = zeros(Float64, nrow(df_species))

for (key, val) ∈ init_dict
    try
        println("$(key): $(val)")
        idx = df_species[df_species[!, "varname"] .== key, :idx_species][1]
        if idx ≤ n_integrated
            #idx = findfirst(x -> x == key, species)
            println("idx: ", idx)
            u₀[idx] = val
        end
    catch e
        println(e)
    end
end


# next update with initial values of our measured species
df_nd_init = df_number_densities[1, Not([measurements_to_ignore..., :t, :w_ap])]
for name ∈ names(df_nd_init)
    if name ∈ df_species.varname
        println("$(name)")
        idx = df_species[df_species.varname .== name, :].idx_species[1]
        println("\tOld val: ", u₀[idx])
        if idx ≤ n_integrated
            u₀[idx] = df_nd_init[name]
            println("\tNew val: ", u₀[idx])
        end
    else
        println("$(name) not in species list")
    end
end

sum(u₀ .!= 0.0)

idx_nonzero = findall(u₀ .!= 0.0)
println(df_species.varname[idx_nonzero])
println(names(df_nd_init))


df_u0 = DataFrame(:u0 => u₀)
CSV.write(joinpath(outpath, "mechanism", "u0.csv"), df_u0)




# -----
# 6. generate reaction dbs, derivative dbs, and jacobian dbs
# -----
@info "Generate derivative and jacobian lists"
# read in reaction databases
const db_bimol = read_bimol(joinpath(outpath, "mechanism", "bimol.json"));
const db_trimol = read_trimol(joinpath(outpath, "mechanism", "trimol.json"));
const db_photo = read_fitted_photolysis(joinpath(outpath, "mechanism", "fitted_photolysis.json"));


# need to fix this
N = generate_stoich_mat(df_species, db_bimol, db_trimol, db_photo, outpath);



const derivatives_bimol = get_bimolecular_derivatives(db_bimol, df_species)
const derivatives_trimol = get_trimolecular_derivatives(db_trimol, df_species)
const derivatives_photo = get_photolysis_derivatives(db_photo, df_species)

const jacobian_terms_bimol = get_bimolecular_jacobian_terms(derivatives_bimol)
const jacobian_terms_trimol = get_trimolecular_jacobian_terms(derivatives_trimol)
const jacobian_terms_photo = get_photolysis_jacobian_terms(derivatives_photo)



# -----
# 7. pre-alloc reaction rate vectors
# -----
@info "Loading reaction rate tables..."
const K_bimol = readdlm(joinpath(outpath, "mechanism", "K_bimol.csv"), ',')
const K_trimol = readdlm(joinpath(outpath, "mechanism", "K_trimol.csv"), ',')
const K_photo = readdlm(joinpath(outpath, "mechanism", "K_photo.csv"), ',')


# -----
# 8. compute non-integrated species concentrations
# -----
@info "Compute non-integrated species concentrations"

df_params = CSV.read(joinpath(data_basepath, "number_densities", collection_id, "number_densities.csv"), DataFrame)

const ts = df_params.t
const temperatures = df_params.temperature
const pressures = df_params.pressure

# generate constant indices for N2, O2, Ar, and other non-integrated species and use those to do the update.
nrow(df_species)
println(length(idx_noint))
df_species[idx_noint, :]


length(u₀)


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



@info "Testing functions for getting concentration"
@benchmark get_concentration(1,1,u₀, U_noint, n_integrated) # 23 ns
@benchmark get_concentration(n_integrated+1, 1, u₀, U_noint, n_integrated)

@info "Testing derivative update"

du = copy(u₀)
prod_temp = 1.0

@benchmark update_derivative!(1, du, u₀, derivatives_bimol[1], K_bimol, prod_temp, U_noint, n_integrated)
@benchmark update_derivative!(1, du, u₀, derivatives_trimol[1], K_trimol, prod_temp, U_noint, n_integrated)
@benchmark update_derivative!(1, du, u₀, derivatives_photo[1], K_photo, prod_temp, U_noint, n_integrated)



# -----
# 9. generate rhs func
# -----
@info "Generating right-hand-side function"
write_rhs_func(outpath)
include(joinpath(outpath, "mechanism", "rhs.jl"))

# -----
# 10. generate jacobian func
# -----
@info "Generating jacobian function"
write_jac_func(outpath)
include(joinpath(outpath, "mechanism", "jacobian.jl"))



# -----
# 11. generate jacobian prototype
# -----
@info "Generating sparse jacobian prototype"
n_species = nrow(df_species)

jac_prototype = generate_jac_prototype(
    jacobian_terms_bimol,
    jacobian_terms_bimol,
    jacobian_terms_bimol,
    n_species
)



# -----
# 12. Test out integration
# -----
test_u₀ = copy(u₀)
test_jac = zeros(length(u₀), length(u₀))

# test jacobian update
@info "Testing jacobian update"
@benchmark update_jacobian!(1, test_jac, u₀, jacobian_terms_bimol[1], K_bimol, prod_temp, U_noint, n_integrated)


du = zeros(length(u₀))
u₀_test = copy(u₀) .+ 1

rhs!(du, u₀_test, nothing, ts[1])

@assert !all(du .== 0.0)  # make sure we are actually updating the du


@info "Testing rhs and jacobian functions"
@benchmark rhs!(du, u₀, nothing, ts[1])  # 20 μs
@benchmark jac!(test_jac, u₀, nothing, ts[1])  # 48 μs

# 81, 82, 85, 87, 88
# idx_bad = [81, 82, 85, 87, 88]
# df_species[idx_bad,:]

@info "Test integration of model"

const tspan = (ts[1], ts[end])
test_u₀ = copy(u₀) .+ 1000.0
# test_u₀[idx_bad] .= 0.0
# test_u₀ .+= 1.0
# test_u₀[[41,85,86,87,88]] .= 0.0

# define ODE function
fun = ODEFunction(rhs!; jac=jac!)
fun2 = ODEFunction(rhs!; jac=jac!, jac_prototype=jac_prototype)

# fun = ODEFunction(rhs!) # ; jac=jac!) #, jac_prototype=jac_prototype)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, test_u₀, tspan)
ode_prob2 = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun2, test_u₀, tspan)

# @benchmark sol = solve(ode_prob, CVODE_BDF(); saveat=15.0)  # 3.702 s
# @benchmark sol = solve(ode_prob, QNDF(); saveat=15.0)
# @benchmark sol = solve(ode_prob; saveat=15.0)  # 83.530 ms
# @benchmark sol = solve(ode_prob2, TRBDF2(); saveat=Δt_step)  # 405 ms, i.e. not sparse enough for performance gain

#  @benchmark sol = solve(ode_prob, TRBDF2(); saveat=Δt_step)  # 161 ms

# @benchmark solve(ode_prob; alg_hints=[:stiff], saveat=Δt_step, reltol=1e-3, abstol=1e-3) # 111 ms

sol = solve(ode_prob; alg_hints=[:stiff], saveat=Δt_step) #, reltol=1e-3, abstol=1e-3)# , reltol=1e-5, abstol=1e-5)  # 161 ms


println("min val: ", minimum(sol[:,:]))

@assert minimum(sol[:,:]) > -1.0 "Minimum value > -1.0"


# test out jacobian calculation
@info "Testing model jacobian via ForwardDiff.jl"
using ForwardDiff, DiffResults
using SciMLSensitivity

function test_f(u_now)
    _prob = remake(ode_prob, u0=u_now, tspan=(ts[1], ts[2]))
    solve(_prob; alg_hints=[:stiff], reltol=1e-3, abstol=1e-3, dense=false,save_everystep=false,save_start=false, sensealg=:QuadratureAdjoint)[:,end]
end

result = DiffResults.JacobianResult(u₀);
result = ForwardDiff.jacobian!(result, test_f, u₀);


@benchmark ForwardDiff.jacobian!(result, test_f, u₀) # 1.5 s



result.value  # the value
result.derivs[1]  # the jacobian



@info "Testing observation operator and jacobian"
idx_meas = Int[1, 2, 3]
idxs_positive
idxs_negative

h = zeros(length(idx_meas)+2)
dhdu = zeros(length(h), length(u₀))

@benchmark Obs!(h, u₀, idx_meas, idxs_positive, idxs_negative)  # 172 ns
@benchmark JObs!(dhdu, u₀, idx_meas, idxs_positive, idxs_negative)  # 86 ns

dhdu
