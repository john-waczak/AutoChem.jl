using CSV, DataFrames
using ProgressMeter
using CairoMakie
using MintsMakieRecipes
set_theme!(mints_theme)

include("data-utils.jl")


basepath = "../../assets/mpi-mainz-uvviz/raw"
outpath = "../../assets/mpi-mainz-uvviz/joined"
if !isdir(outpath)
    mkpath(outpath)
    mkpath(joinpath(outpath, "cross-sections"))
    mkpath(joinpath(outpath, "quantum-yields"))
end

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





# 3. Process Cross-sections
@showprogress for σf ∈ σfs
    try
        df = collate_summary_σ(σf, basepath)
        CSV.write(joinpath(outpath, "cross-sections", df.formula[1] * ".csv"), df)
    catch e
        println(σf)
        println(e)
    end

end


# 4. Process Quantum-yields
@showprogress for Φf ∈ Φfs
    try
        df = collate_summary_Φ(Φf, basepath)
        CSV.write(joinpath(outpath, "quantum-yields", df.reaction[1] * ".csv"), df)
    catch e
        println(Φf)
        println(e)
    end

end

