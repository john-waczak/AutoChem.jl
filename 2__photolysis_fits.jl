using AutoChem
using DelimitedFiles, CSV, DataFrames
using JSON
using CairoMakie
using MintsMakieRecipes
using DataInterpolations
using BenchmarkTools

set_theme!(mints_theme)


# set up model output directory
model_name = "autochem-w-ions"
#model_name = "qroc-methane-intel"
model_path= "models/"

data_basepath = "data/intertek-emergency-testing"
collection_id = "empty"

λs_spec = CSV.read("./assets/hr4000_wavelengths.txt", DataFrame).λ
λ_min = 225.0
λ_max = 1000.0

λs_spec = λs_spec[λs_spec .≥ λ_min .&& λs_spec .≤ λ_max]


# for now, let's build up the fits using our current mechanism files
df_species = CSV.read(joinpath(model_path, model_name, "mechanism", "species.csv"), DataFrame)
db_photolysis = read_photolysis(joinpath(model_path, model_name, "mechanism", "photolysis.json"))

autochem_basepath = joinpath(AutoChem.assets_path, "autochem", "data", "CrossSections")
readdir(autochem_basepath)

σ_basepath = joinpath(AutoChem.assets_path, "mpi-mainz-uvviz", "joined", "cross-sections")
Φ_basepath = joinpath(AutoChem.assets_path, "mpi-mainz-uvviz", "joined", "quantum-yields")



db = FittedPhotolysisReaction[]


for idx ∈ 1:length(db_photolysis)
    # load
    data_autochem = readdlm(joinpath(autochem_basepath, db_photolysis[idx].autochem_files[1]))
    λ_ac = Float64.(data_autochem[2:end,2])
    σ_ac = Float64.(data_autochem[2:end,3])
    Φ_ac = Float64.(data_autochem[2:end,4])

    # fit
    itp = LinearInterpolation(σ_ac, λ_ac)
    σ_out = itp.(λs_spec)

    itp = LinearInterpolation(Φ_ac, λ_ac)
    Φ_out = itp.(λs_spec)

    # save
    push!(db, FittedPhotolysisReaction(
        idx,
        db_photolysis[idx].reactants,
        db_photolysis[idx].products,
        db_photolysis[idx].prod_stoich,
        λs_spec,
        σ_out,
        Φ_out
    ))
end


# save the resulting database
open(joinpath(model_path, model_name, "mechanism", "fitted_photolysis.json"), "w") do f
    JSON.print(f, db, 2)
end



# precompute photolysis rates for each reaction for each time

db_bimol = read_bimol(joinpath(model_path, model_name, "mechanism", "bimol.json"))
db_trimol = read_trimol(joinpath(model_path, model_name, "mechanism", "trimol.json"))
db_photo = read_fitted_photolysis(joinpath(model_path, model_name, "mechanism", "fitted_photolysis.json"))

Is = readdlm(joinpath(data_basepath, "rates", collection_id, "Intensities.csv"), ',')

df_params = CSV.read(joinpath(data_basepath, "number_densities", collection_id, "number_densities.csv"), DataFrame)

Ts = df_params.temperature
Ps = df_params.pressure
Ms = df_params.M

@assert all(Ms .== M.(Ts, Ps))

K_bimol = zeros(Float64, length(db_bimol), nrow(df_params))
K_trimol = zeros(Float64, length(db_trimol), nrow(df_params))
K_photo = zeros(Float64, length(db_photo), nrow(df_params))

for j ∈ axes(K_bimol,2), i ∈ axes(K_bimol,1)
    K_bimol[i,j] = db_bimol[i](Ts[j], Ps[j], Ms[j])
end

for j ∈ axes(K_trimol,2), i ∈ axes(K_trimol,1)
    K_trimol[i,j] = db_trimol[i](Ts[j], Ps[j], Ms[j])
end

for j ∈ axes(K_photo,2), i ∈ axes(K_photo,1)
    K_photo[i,j] = db_photo[i](Ts[j], Ps[j], Is[:,j])
end


writedlm(joinpath(model_path, model_name, "mechanism", "K_bimol.csv"), K_bimol, ',')
writedlm(joinpath(model_path, model_name, "mechanism", "K_trimol.csv"), K_trimol, ',')
writedlm(joinpath(model_path, model_name, "mechanism", "K_photo.csv"), K_photo, ',')


# # ----------------------------
# # Br2 + Photon ⟶ Br + Br
# # ----------------------------

# idx = 1

# println("source: ", db_photolysis[idx].source)
# println("AutoChem files: ", db_photolysis[idx].autochem_files[1])
# println("σ files: ", db_photolysis[idx].crosssection_files[1])
# println("Φ files: ", db_photolysis[idx].quantumyield_files[1])


# data_autochem = readdlm(joinpath(autochem_basepath, db_photolysis[idx].autochem_files[1]))
# λ_ac = Float64.(data_autochem[2:end,2])
# σ_ac = Float64.(data_autochem[2:end,3])
# Φ_ac = Float64.(data_autochem[2:end,4])

# df_σ = CSV.read(joinpath(σ_basepath, db_photolysis[idx].crosssection_files[1]), DataFrame)
# df_Φ = DataFrame()

# # visualize σ

# fig = Figure();
# ax = Axis(fig[1,1], xlabel="λ (nm)", ylabel="σ (cm²)", limits=(nothing, nothing, 0, nothing), ytickformat="{:.3e}")
# ax2 = Axis(fig[1,2], xlabel="λ (nm)", ylabel="log10(σ)")

# sc_ac = scatter!(ax, λ_ac, σ_ac)
# ls_ac = lines!(ax, λ_ac, σ_ac; linewidth=3)
# sc_db = scatter!(ax, df_σ.λs, df_σ.σs, alpha=0.1)

# σ_log = log10.(σ_ac)
# idx_inf = isinf.(σ_log)
# σ_log[idx_inf] .= -22

# σ_log2 = log10.(df_σ.σs)
# idx_inf = isinf.(σ_log2)
# σ_log2[idx_inf] .= -22

# sc_ac_lg = scatter!(ax2, λ_ac, σ_log)
# ls_ac_lg = lines!(ax2, λ_ac, σ_log; linewidth=3)
# sc_db_lg = scatter!(ax2, df_σ.λs, σ_log2; alpha=0.1)

# # fit
# itp = LinearInterpolation(σ_ac, λ_ac)
# σ_out = itp.(λs_spec)

# itp = LinearInterpolation(Φ_ac, λ_ac)
# Φ_out = itp.(λs_spec)


# lines!(ax, λs_spec, σ_out; color=:purple)

# fig


# # save
# push!(db, FittedPhotolysisReaction(
#     idx,
#     db_photolysis[idx].reactants,
#     db_photolysis[idx].products,
#     db_photolysis[idx].prod_stoich,
#     λs_spec,
#     σ_out,
#     Φ_out
# ))


