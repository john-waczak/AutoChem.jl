# now we need a function to generate a DataFrame with all data
function collate_summary_σ(rootpath, basepath)
    summary_df = CSV.read(joinpath(rootpath, "summary.csv"), DataFrame)

    T₁ = []
    T₂ = []
    author_year = []
    name = []
    doi = []
    formula = []
    download_url = []

    # these will be updated with data from each individual file
    λ = []
    σ = []
    Δσ = []
    source_idx = []

    for i ∈ 1:nrow(summary_df)
        try
            df = CSV.read(joinpath(basepath, summary_df.fname[i]), DataFrame)
            push!(λ, df.λ)
            push!(σ, df.σ)
            push!(Δσ, df.Δσ)
            push!(source_idx, [i for j ∈ 1:nrow(df)])

            push!(T₁, [summary_df.T1[i] for _ ∈ 1:nrow(df)])
            push!(T₂, [summary_df.T2[i] for _ ∈ 1:nrow(df)])
            push!(author_year, [summary_df.author_year[i] for _ ∈ 1:nrow(df)])
            push!(name, [summary_df.name[i] for _ ∈ 1:nrow(df)])
            push!(doi, [summary_df.doi[i] for _ ∈ 1:nrow(df)])
            push!(formula, [summary_df.formula[i] for _ ∈ 1:nrow(df)])
            push!(download_url, [summary_df.download_url[i] for _ ∈ 1:nrow(df)])

        catch e
            println("Couldn't load $(summary_df.fname[i])")
            println(e)
        end
    end


    return DataFrame((
        T1 = vcat(T₁...),
        T2 = vcat(T₂...),
        author_year = vcat(author_year...),
        name = vcat(name...),
        doi = vcat(doi...),
        formula = vcat(formula...),
        download_url = vcat(download_url...),
        λs = vcat(λ...),
        σs = vcat(σ...),
        Δσs = vcat(Δσ...)
    ))

    # # note that name and formula are the same for each row
    # return Float64.(T₁), Float64.(T₂), author_year, comments, name[1], doi, formula[1], fname, download_url, vcat(λ...), vcat(σ...), vcat(Δσ...), vcat(source_idx...)


end



function collate_summary_ϕ(rootpath, basepath)
    summary_df = CSV.read(joinpath(rootpath, "summary.csv"), DataFrame)

    T₁ = []
    T₂ = []
    author_year = []
    name = []
    doi = []
    reaction = []
    download_url = []

    # these will be updated with data from each individual file
    λ = []
    Φ = []
    ΔΦ = []
    species = []
    source_idx = []

    for i ∈ 1:nrow(summary_df)
        try
            df = CSV.read(joinpath(basepath, summary_df.fname[i]), DataFrame)
            push!(λ, df.λ)
            push!(Φ, df.Φ)
            push!(ΔΦ, df.ΔΦ)
            push!(species, df.species)
            push!(source_idx, [i for j ∈ 1:nrow(df)])

            push!(T₁, [summary_df.T1[i] for _ ∈ 1:nrow(df)])
            push!(T₂, [summary_df.T2[i] for _ ∈ 1:nrow(df)])
            push!(author_year, [summary_df.author_year[i] for _ ∈ 1:nrow(df)])
            push!(name, [summary_df.name[i] for _ ∈ 1:nrow(df)])
            push!(doi, [summary_df.doi[i] for _ ∈ 1:nrow(df)])
            push!(reaction, [summary_df.reaction[i] for _ ∈ 1:nrow(df)])
            push!(download_url, [summary_df.download_url[i] for _ ∈ 1:nrow(df)])

        catch e
            println("Couldn't load $(summary_df.fname[i])")
            println(e)
        end
    end

    # note that name and formula are the same for each row
    # return Float64.(T₁), Float64.(T₂), author_year, comments, name[1], doi, reaction[1], fname, download_url, vcat(λ...), vcat(Φ...), vcat(ΔΦ...), vcat(species...), vcat(source_idx...)
    return DataFrame((
        T1 = vcat(T₁...),
        T2 = vcat(T₂...),
        author_year = vcat(author_year...),
        name = vcat(name...),
        doi = vcat(doi...),
        reaction = vcat(reaction...),
        download_url = vcat(download_url...),
        λs = vcat(λ...),
        Φs = vcat(Φ...),
        ΔΦs = vcat(ΔΦ...)
    ))



end


