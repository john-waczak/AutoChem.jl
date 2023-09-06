using AutoChem
using DelimitedFiles, CSV, DataFrames
using JSON
using CairoMakie
using MintsMakieRecipes
set_theme!(mints_theme)


# set up model output directory
model_name = "autochem-w-ions"
model_path= "models/"


# for now, let's build up the fits using our current mechanism files
df_species = CSV.read(joinpath(model_path, model_name, "mechanism", "species.csv"), DataFrame)
db_photolysis = read_photolysis(joinpath(model_path, model_name, "mechanism", "photolysis.json"))

autochem_basepath = joinpath(AutoChem.assets_path, "autochem", "data", "CrossSections")
readdir(autochem_basepath)

σ_basepath = joinpath(AutoChem.assets_path, "mpi-mainz-uvviz", "joined", "cross-sections")
Φ_basepath = joinpath(AutoChem.assets_path, "mpi-mainz-uvviz", "joined", "quantum-yields")



db = FittedPhotolysisReaction[]


# ----------------------------
# Br2 + Photon ⟶ Br + Br
# ----------------------------

idx = 1

println("source: ", db_photolysis[idx].source)
println("AutoChem files: ", db_photolysis[idx].autochem_files[1])
println("σ files: ", db_photolysis[idx].crosssection_files[1])
println("Φ files: ", db_photolysis[idx].quantumyield_files[1])


data_autochem = readdlm(joinpath(autochem_basepath, db_photolysis[idx].autochem_files[1]))
λ_ac = Float64.(data_autochem[2:end,2])
σ_ac = Float64.(data_autochem[2:end,3])
Φ_ac = Float64.(data_autochem[2:end,4])

df_σ = CSV.read(joinpath(σ_basepath, db_photolysis[idx].crosssection_files[1]), DataFrame)
df_Φ = DataFrame()

# fit σ

fig = Figure();
ax = Axis(fig[1,1], xlabel="λ (nm)", ylabel="σ (cm²)", limits=(nothing, nothing, 0, nothing), ytickformat="{:.3e}")
ax2 = Axis(fig[1,2], xlabel="λ (nm)", ylabel="log10(σ)")

sc_ac = scatter!(ax, λ_ac, σ_ac)
ls_ac = lines!(ax, λ_ac, σ_ac; linewidth=3)
sc_db = scatter!(ax, df_σ.λs, df_σ.σs, alpha=0.1)

σ_log = log10.(σ_ac)
idx_inf = isinf.(σ_log)
σ_log[idx_inf] .= -22

σ_log2 = log10.(df_σ.σs)
idx_inf = isinf.(σ_log2)
σ_log2[idx_inf] .= -22

sc_ac_lg = scatter!(ax2, λ_ac, σ_log)
ls_ac_lg = lines!(ax2, λ_ac, σ_log; linewidth=3)
sc_db_lg = scatter!(ax2, df_σ.λs, σ_log2; alpha=0.1)

fig


# db[idx] = FittedPhotolysisReaction(
#     idx,
#     db_photolysis[idx].reactants,
#     db_photolysis[idx].products
# )
