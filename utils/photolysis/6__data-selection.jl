using CSV, DataFrames
using HDF5
using CairoMakie
using MintsMakieRecipes
set_theme!(mints_theme)

include("data-utils.jl")


basepath = "../../assets/mpi-mainz-uvviz/raw"
outpath = "../../assets/mpi-mainz-uvviz/joined"
readdir(basepath)

# 1. Cross-section data

σfs = []
for (root, dirs, files) ∈ walkdir(joinpath(basepath, "cross-sections"))
    for file ∈ files
        if endswith(file, ".csv")
            fpath = joinpath(split(joinpath(root, file), "/")[1:end-1])
            push!(σfs, fpath)
        end
    end
end

unique!(σfs)


# 2. Quantum-Yield data

Φfs = []
for (root, dirs, files) ∈ walkdir(joinpath(basepath, "quantum-yields"))
    for file ∈ files
        if endswith(file, ".csv")
            fpath = joinpath(split(joinpath(root, file), "/")[1:end-1])
            push!(Φfs, fpath)
        end
    end
end

unique!(Φfs)



# 3. Process Cross-sections
rootpath = σfs[1]
@assert ispath(rootpath)

df = collate_summary_σ(rootpath, basepath)
names(df)

df.download_url[1]
df.download_url[end]

collate_summary_ϕ(Φfs[1], basepath)
