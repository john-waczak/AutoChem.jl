ENV["GKSwstype"] = 100

@info "Setting up Julia Environment..."
using Pkg
Pkg.activate(".")
Pkg.instantiate()
@info "\t...finished"

using AutoChem
using DelimitedFiles, CSV, DataFrames
using JSON
using DataInterpolations
using BenchmarkTools
using ArgParse

# using CairoMakie
# using MintsMakie
# set_theme!(mints_them)

# ^ we should create plots for σ,Φ, and reaction rates and link to plots

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
            default = "methane"
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
Δt_step = parsed_args[:time_step]  # time step in minutes


model_path= "models/"

@info "Setting up file paths..."

outpath = joinpath(model_path, model_name, "runs", collection_id)
@assert ispath(outpath)

docs_path = joinpath(model_path, model_name, "docs")
@assert ispath(docs_path)



@info "Loading in spectrometer wavelengths"
λs_spec = CSV.read("./assets/hr4000_wavelengths.txt", DataFrame).λ
λ_min = 225.0
λ_max = 1000.0

λs_spec = λs_spec[λs_spec .≥ λ_min .&& λs_spec .≤ λ_max]


# for now, let's build up the fits using our current mechanism files
@info "Loading species file and photolysis reaction database"
df_species = CSV.read(joinpath(outpath, "mechanism", "species.csv"), DataFrame)
db_photolysis = read_photolysis(joinpath(outpath, "mechanism", "photolysis.json"))

@info "loading σ and Φ data"
autochem_basepath = joinpath(AutoChem.assets_path, "autochem", "data", "CrossSections")
σ_basepath = joinpath(AutoChem.assets_path, "mpi-mainz-uvviz", "joined", "cross-sections")
Φ_basepath = joinpath(AutoChem.assets_path, "mpi-mainz-uvviz", "joined", "quantum-yields")




@info "Generating new photolysis db with fitted σ,Φ"
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
open(joinpath(outpath, "mechanism", "fitted_photolysis.json"), "w") do f
    JSON.print(f, db, 2)
end




# precompute photolysis rates for each reaction for each time
@info "Precomputing reaction rate coefficient tables..."
db_bimol = read_bimol(joinpath(outpath, "mechanism", "bimol.json"))
db_trimol = read_trimol(joinpath(outpath, "mechanism", "trimol.json"))
db_photo = read_fitted_photolysis(joinpath(outpath, "mechanism", "fitted_photolysis.json"))

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

@info "Saving outputs"
writedlm(joinpath(outpath, "mechanism", "K_bimol.csv"), K_bimol, ',')
writedlm(joinpath(outpath, "mechanism", "K_trimol.csv"), K_trimol, ',')
writedlm(joinpath(outpath, "mechanism", "K_photo.csv"), K_photo, ',')





























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


