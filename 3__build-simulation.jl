using AutoChem
using DelimitedFiles, CSV, DataFrames
using JSON
using CairoMakie
using MintsMakieRecipes
using BenchmarkTools

using DifferentialEquations
using Sundials, OrdinaryDiffEq


set_theme!(mints_theme)


# set up model output directory
data_basepath = "data/intertek-emergency-testing"
collection_id = "empty"
unc_ext = "_std"
model_name = "autochem-w-ions"
model_path= "models"
Δt_step = 15.0


if !ispath(joinpath(model_path, model_name, "figures"))
    joinpath(model_path, model_name, "figures")
end



# -----
# 1. fetch species and create indices
# -----
df_species = CSV.read(joinpath(model_path, model_name, "mechanism", "species.csv"), DataFrame)


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
    model_name=model_name,
    replacement_dict=replacement_dict
)



df_params = CSV.read(joinpath(model_path, model_name, "mechanism", "state_parameters.csv"), DataFrame)
df_number_densities = CSV.read(joinpath(model_path, model_name, "mechanism", "number_densities.csv"), DataFrame)



# -----
# 3. Generate lookup table for Intensities
# -----
Is = readdlm(joinpath(data_basepath, "rates", collection_id, "Intensities.csv"), ',', Float64)


# visualize the intensities
db_photo = read_fitted_photolysis(joinpath(model_path, model_name, "mechanism", "fitted_photolysis.json"));

fig = Figure()
ax = Axis(fig[1,1], xlabel="time (s)", ylabel="λ (nm)")
hm = heatmap!(ax, df_params.t, db_photo[1].λs, Is')
cb = Colorbar(fig[1,2], hm; label="Irradiance (photons  s⁻¹ cm⁻² nm⁻¹)")
fig

save(joinpath(model_path, model_name, "figures", "Intensities.png"), fig)
#save(joinpath(model_path, model_name, "figures", "Intensities.svg"), fig)
save(joinpath(model_path, model_name, "figures", "Intensities.eps"), fig)
save(joinpath(model_path, model_name, "figures", "Intensities.pdf"), fig)

fig = Figure()
ax = Axis(fig[1,1], xlabel="time (s)", ylabel="λ (nm)")
hm = heatmap!(ax, df_params.t, db_photo[1].λs, log10.(Is'))
cb = Colorbar(fig[1,2], hm; label="log10(Irradiance)")
fig

save(joinpath(model_path, model_name, "figures", "log10-Intensities.png"), fig)
#save(joinpath(model_path, model_name, "figures", "log10-Intensities.svg"), fig)
save(joinpath(model_path, model_name, "figures", "log10-Intensities.eps"), fig)
save(joinpath(model_path, model_name, "figures", "log10-Intensities.pdf"), fig)

fig

# test out computation of photolysis rate
@benchmark db_photo[1](df_params.temperature[1], df_params.pressure[1], Is[:,1])  # 13.174 μs

# -----
# 4. generate indices for ro2 sum
# -----

# skipping for now as it's not relevant except for MCM



# -----
# 5. generate sane initial conditions
# -----
init_path = "./assets/initial_concentrations/full.txt"
@assert isfile(init_path) "Cant find ./assets/initial_concentrations/full.txt"

init_dict = generate_init_dict(init_path, df_params.M[1])


for (key,val) ∈ init_dict
    println(key,":\t",val)
end

# this should be put in a config file
measurements_to_ignore = [:C2H6, :SO2]  # skip any with nans or negative values

u₀ = zeros(Float64, nrow(df_species))

for (key, val) ∈ init_dict
    try
        println("$(key): $(val)")
        idx = df_species[df_species[!, "varname"] .== key, :idx_species][1]
        #idx = findfirst(x -> x == key, species)
        println("idx: ", idx)
        u₀[idx] = val

    catch e

    end
end

# next update with initial values of our measured species
df_nd_init = df_number_densities[1, Not([measurements_to_ignore..., :t, :w_ap])]
for name ∈ names(df_nd_init)
    if name ∈ df_species.varname
        println("$(name)")
        idx = df_species[df_species.varname .== name, :].idx_species[1]
        println("\tOld val: ", u₀[idx])
        u₀[idx] = df_nd_init[name]
        println("\tNew val: ", u₀[idx])
    else
        println("$(name) not in species list")
    end
end

df_species.varname
sum(u₀ .!= 0.0)

df_species

for i ∈ 1:nrow(df_species[idx_integrated,:])
    if u₀[i] != 0.0
        println(df_species.varname[i], "\t", u₀[i])
    end
end

df_u0 = DataFrame(:u0 => u₀)
CSV.write(joinpath(model_path, model_name, "mechanism", "u0.csv"), df_u0)

# -----
# 6. generate stoichiometry matrix for later visualization
# -----

# read in reaction databases
db_bimol = read_bimol(joinpath(model_path, model_name, "mechanism", "bimol.json"));
db_trimol = read_trimol(joinpath(model_path, model_name, "mechanism", "trimol.json"));
db_photo = read_fitted_photolysis(joinpath(model_path, model_name, "mechanism", "fitted_photolysis.json"));


# need to fix this
N = generate_stoich_mat(df_species, db_bimol, db_trimol, db_photo);


# -----
# 7. generate rhs func
# -----

const derivatives_bimol = get_bimolecular_derivatives(db_bimol, df_species)
const derivatives_trimol = get_trimolecular_derivatives(db_trimol, df_species)
const derivatives_photo = get_photolysis_derivatives(db_photo, df_species)


# generate constant indices for N2, O2, Ar, and other non-integrated species and use those to do the update.
df_species[df_species.is_integrated .== 2, :]

const idx_Ar = findfirst(df_species.varname .== "Ar")
const idx_O2 = findfirst(df_species.varname .== "O2")
const idx_N2 = findfirst(df_species.varname .== "N2")
const idx_M = findfirst(df_species.varname .== "m")


# pre-alloc raction rate vectors
const K_bimol = readdlm(joinpath(model_path, model_name, "mechanism", "K_bimol.csv"), ',')
const K_trimol = readdlm(joinpath(model_path, model_name, "mechanism", "K_trimol.csv"), ',')
const K_photo = readdlm(joinpath(model_path, model_name, "mechanism", "K_photo.csv"), ',')

minimum(K_bimol)
maximum(K_bimol)
mean(K_bimol)
median(K_bimol)

df_params = CSV.read(joinpath(data_basepath, "number_densities", collection_id, "number_densities.csv"), DataFrame)

const ts = df_params.t
const temperatures = df_params.temperature
const pressures = df_params.pressure

du = zeros(Float64, size(u₀))

nrow(df_params)
size(K_bimol)

prod_temp = 1.0

@benchmark update_derivative!(1, du, u₀, derivatives_bimol[1], K_bimol, prod_temp)
@benchmark update_derivative!(1, du, u₀, derivatives_trimol[1], K_trimol, prod_temp)
@benchmark update_derivative!(1, du, u₀, derivatives_photo[1], K_photo, prod_temp)



function rhs!(du, u, p, t)
    # get time value and index
    idx_t = get_time_index(t, Δt_step, ts[1])

    # update values for number density, O2, N2, Ar
    u[idx_M] = M(temperatures[idx_t], pressures[idx_t])
    u[idx_O2] = O2(temperatures[idx_t], pressures[idx_t])
    u[idx_N2] = N2(temperatures[idx_t], pressures[idx_t])
    u[idx_Ar] = Ar(temperatures[idx_t], pressures[idx_t])

    # set derivatives to zero
    du .= 0.0

    # loop over bimol derivatives
    prod_temp = 1.0
    @inbounds for i ∈ eachindex(derivatives_bimol)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives_bimol[i],
            K_bimol,
            prod_temp
        )
    end

    # loop over trimol derivatives
    prod_temp = 1.0
    @inbounds for i ∈ eachindex(derivatives_trimol)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives_trimol[i],
            K_trimol,
            prod_temp
        )
    end


    # loop over photolysis derivatives
    prod_temp = 1.0
    @inbounds for i ∈ eachindex(derivatives_photo)
        prod_temp = 1.0 # <-- start fresh for each derivative
        update_derivative!(
            idx_t,
            du,
            u,
            derivatives_photo[i],
            K_photo,
            prod_temp
        )
    end
end



test_u₀ = copy(u₀)
test_u₀ .+ 1e4

@benchmark rhs!(du, u₀, nothing, ts[1])


const tspan = (ts[1], ts[end])

# function get_species_idx()

# end


write_rhs_func(model_name=model_name)
include("models/$model_name/rhs.jl")

# create derivative struct
# create functions to parse each reaction type to derivative struct
# create callable rhs function


# -----
# 8. generate jacobian func
# -----
write_jac_func(model_name=model_name)
include("models/$model_name/jacobian.jl")



fun = ODEFunction(rhs!; jac=jac!, jac_prototype=jac_prototype)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan)
solve(ode_prob, QNDF(); saveat=15.0, reltol=1e-3, abstol=1e-3)

