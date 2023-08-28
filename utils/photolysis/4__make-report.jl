using CSV, DataFrames
using ProgressMeter
using CairoMakie
using MintsMakieRecipes
set_theme!(mints_theme)

basepath   = "../../assets/mpi-mainz-uvviz/raw"
basepath_σ = "../../assets/mpi-mainz-uvviz/raw/cross-sections"
basepath_Φ = "../../assets/mpi-mainz-uvviz/raw/quantum-yields"
@assert ispath(basepath_σ)
@assert ispath(basepath_Φ)

outpath_figures = joinpath(basepath, "figures")
if !ispath(outpath_figures)
    mkpath(outpath_figures)
end

outpath_σ = "../../assets/mpi-mainz-uvviz/raw/cross-section_summary.qmd"
outpath_Φ = "../../assets/mpi-mainz-uvviz/raw/quantum-yield_summary.qmd"



# load a sample folder for cross-sections:
σ_paths = []
for (root, dirs, files) ∈ walkdir(basepath_σ)
    for file ∈ files
        if endswith("summary.csv", file)
            push!(σ_paths, joinpath(root, file))
        end
    end
end


open(outpath_σ, "w") do f
    @showprogress for σ_path ∈ σ_paths
        try
            df = CSV.read(σ_path, DataFrame)

            species = df.formula[1]
            speciesname = df.name[1]

            println(f, "# $(speciesname) ($(species))\n")

            fig = Figure(resolution=(1200, 500));
            ax = Axis(fig[1,1], xlabel="λ (nm)", ylabel="σ (cm²)", title="$(df.formula[1])");
            ax2 = Axis(fig[1,2], xlabel="λ (nm)", ylabel="log10(σ)", title="$(df.formula[1])");
            scatters = []
            authors = []
            for row ∈ eachrow(df)
                println(f, "author:\t$(row.author_year)\n")
                println(f, "- source:\t$(row.download_url)")
                println(f, "- comments:\t$(row.comments)")
                println(f, "- doi:\t$(row.doi)")
                println(f, "- fname:\t$(row.fname)")
                println(f, "- temperature:\t$(row.T1)")
                println(f, "\n")

                df_σ = CSV.read(joinpath(basepath, row.fname), DataFrame)

                df_σ.σ[df_σ.σ .≤ 0.0] .= 0.0;

                s = scatter!(ax, df_σ.λ, df_σ.σ)
                scatter!(ax2, df_σ.λ, log10.(df_σ.σ))
                push!(authors, row.author_year * " T=$(row.T1)")
                push!(scatters, s)
            end
            leg = Legend(
                fig[1, 3],
                scatters,
                authors,
            )

            save(joinpath(outpath_figures, species*"_σ.png"), fig)
            figpath = joinpath("./figures", species*"_σ.png")

            println(f, "![]($(figpath))\n")
        catch e
            println(e)
        end
    end
end



# load a sample folder for cross-sections:
Φ_paths = []
for (root, dirs, files) ∈ walkdir(basepath_Φ)
    for file ∈ files
        if endswith("summary.csv", file)
            push!(Φ_paths, joinpath(root, file))
        end
    end
end


open(outpath_Φ, "w") do f
    @showprogress for Φ_path ∈ Φ_paths
        try
            df = CSV.read(Φ_path, DataFrame)

            rxn = df.reaction[1]
            rxnname = df.name[1]

            println(f, "# $(rxnname) \n")

            fig = Figure();
            ax = Axis(fig[1,1], xlabel="λ (nm)", ylabel="Φ ∈ [0,1]", title="$(df.reaction[1])");
            scatters = []
            authors = []
            for row ∈ eachrow(df)
                println(f, "author:\t$(row.author_year)\n")
                println(f, "- source:\t$(row.download_url)")
                println(f, "- comments:\t$(row.comments)")
                println(f, "- doi:\t$(row.doi)")
                println(f, "- fname:\t$(row.fname)")
                println(f, "- temperature:\t$(row.T1)")
                println(f, "\n")

                df_Φ = CSV.read(joinpath(basepath, row.fname), DataFrame)
                s = scatter!(ax, df_Φ.λ, df_Φ.Φ)
                push!(authors, row.author_year * " T=$(row.T1)")
                push!(scatters, s)
            end
            leg = Legend(
                fig[1, 2],
                scatters,
                authors,
            )

            save(joinpath(outpath_figures, rxnname*"_Φ.png"), fig)
            figpath = joinpath("./figures", rxnname*"_σ.png")

            println(f, "![]($(figpath))\n")
        catch e
            println(e)
        end
    end
end
