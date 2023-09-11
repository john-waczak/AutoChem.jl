using AutoChem
using DelimitedFiles, CSV, DataFrames
using JSON
using CairoMakie
using MintsMakieRecipes
using BenchmarkTools

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

idx_integrated = df_species.is_integrated .== 1
u₀ = zeros(Float64, nrow(df_species[idx_integrated, :]))

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
N



# -----
# 7. generate rhs func
# -----


n_bimol = length(db_bimol)
n_trimol = length(db_trimol)
n_photo = length(db_photo)


res = get_derivative_terms(db_bimol[1], 1)
res = get_derivative_terms(db_trimol[1], 1)
res = get_derivative_terms(db_photo[1], 1)


list = eltype(db_bimol)[]

get_bimolecular_derivatives(db_bimol)
get_trimolecular_derivatives(db_trimol)
get_photolysis_derivatives(db_photo)

# the trick now will be to update the values of u for the non-integrated species. We can do this at the beginning of the function when we fetch the value of temperatures, pressures, etc...

# generate constant indices for N2, O2, Ar, and other non-integrated species and use those to do the update.

df_species[df_species.is_integrated .== 2, :]


idx_Ar = findfirst(df_species.varname .== "Ar")
# and so on for O2, N2, etc...


# function get_species_idx()

# end


# function DerivativeTerms(rxn::BimolecularReaction)
#     dts = DerivativeTerm[]

# end

#     struct BimolecularReaction{T0<:Integer, T1<:AbstractString, T2<:AbstractString, T3<:Real, T4<:Real}<: Reaction
#         idx::T0
#         source::String
#         reactants::AbstractVector{T1}
#         products::AbstractVector{T2}
#         prod_stoich::AbstractVector{T3}
#         a1::T4
#         a2::T4
#         a3::T4
#         a4::T4
#         a5::T4
#         contains_OH::Bool
#         contains_HONO2::Bool
#         contains_CO::Bool
#         all_reactants_HO2::Bool
#     end




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



