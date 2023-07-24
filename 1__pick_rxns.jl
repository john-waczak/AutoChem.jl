using AutoChem
using CSV, DataFrames
using JSON



# set up model output directory
model_name = "autochem-w-ions"
model_path= "models/"

if !ispath(joinpath(model_path, model_name))
    mkpath(joinpath(model_path, model_name))
    mkpath(joinpath(model_path, model_name, "mechanism"))
    mkpath(joinpath(model_path, model_name, "figures"))
    mkpath(joinpath(model_path, model_name, "4d-var"))
    mkpath(joinpath(model_path, model_name, "EKF"))
end

outpath = joinpath(model_path, model_name)

# 1. Create Species List from list of qroc files


qroc_list = [
    "qroc-red-NHMC-nagfor",
    "qroc-methane-ion-nagfor"
]


function join_species(qrocs...)
    dfs = DataFrame[]

    for qroc ∈ qrocs
        try
            species = CSV.File(joinpath("src", "species", "species", qroc * ".csv")) |> DataFrame
            push!(dfs, species)
        catch e
            println(e)
        end
    end

    df_species = vcat(dfs...)

    idx_unique = [findfirst(x->x==name, df_species.varname) for name ∈ unique(df_species.varname)]
    df_species = df_species[idx_unique, :]

    sort!(df_species, :is_integrated)

    return df_species
end



function join_totals(qrocs...)
    dfs = DataFrame[]

    for qroc ∈ qrocs
        try
            species = CSV.File(joinpath("src", "species", "totals", qroc * ".csv")) |> DataFrame
            push!(dfs, species)
        catch e
            println(e)
        end
    end

    df_totals = vcat(dfs...)

    idx_unique = [findfirst(x->x==name, df_totals.name) for name ∈ unique(df_totals.name)]
    df_totals = df_totals[idx_unique, :]

    return df_totals
end



function join_ratios(qrocs...)
    dfs = DataFrame[]

    for qroc ∈ qrocs
        try
            species = CSV.File(joinpath("src", "species", "ratios", qroc * ".csv")) |> DataFrame
            push!(dfs, species)
        catch e
            println(e)
        end
    end

    df_ratios = vcat(dfs...)

    idx_unique = [findfirst(x->x==name, df_ratios.species_a) for name ∈ unique(df_ratios.species_a)]
    df_ratios = df_ratios[idx_unique, :]

    return df_ratios
end



# generate the lists
df_species = join_species(qroc_list...)
df_totals = join_totals(qroc_list...)
df_ratios = join_ratios(qroc_list...)


# save them to the mechanism directory
CSV.write(joinpath(outpath, "mechanism", "species.csv"), df_species)
CSV.write(joinpath(outpath, "mechanism", "totals.csv"), df_totals)
CSV.write(joinpath(outpath, "mechanism", "ratios.csv"), df_ratios)


# now we want to fetch the relevant reactions
function get_rxn_databases(path::String)
    qroc_lines = readlines(path)

    idx_bi = findall([occursin("bimolecular", line) && occursin(".d", line) && !occursin("bulk", line) for line ∈ qroc_lines])


    idx_tri = findall([occursin("trimolecular", line) && occursin(".d", line) && !occursin("bulk", line) for line ∈ qroc_lines])

    idx_phot = findall([occursin("photolysis", line) && occursin(".d", line) && !occursin("bulk", line) for line ∈ qroc_lines])


    bimol_dbs = [split(split(l, "'")[2], ".")[1]*".json" for l ∈ qroc_lines[idx_bi]]
    trimol_dbs = [split(split(l, "'")[2], ".")[1]*".json" for l ∈ qroc_lines[idx_tri]]
    photo_dbs = [split(split(l, "'")[2], ".")[1]*".json" for l ∈ qroc_lines[idx_phot]]


    return bimol_dbs, trimol_dbs, photo_dbs
end




function get_rxn_databases(paths::Vector{String})
    bimol_dbs = String[]
    trimol_dbs = String[]
    photo_dbs = String[]

    for path ∈ paths
        b, t, p = get_rxn_databases(path)
        append!(bimol_dbs, b)
        append!(trimol_dbs, t)
        append!(photo_dbs, p)
    end

    base = "src/json-databases/"

    return base .* "bimol/" .* bimol_dbs, base .* "trimol/" .* trimol_dbs, base .* "photolysis/" .* photo_dbs
end




paths = joinpath.("src/autochem-databases/qroc-scripts", qroc_list)
bimol_dbs, trimol_dbs, photo_dbs = get_rxn_databases(paths);

bimol_dbs
isfile.(bimol_dbs)

trimol_dbs
isfile.(trimol_dbs)

photo_dbs
isfile.(photo_dbs)




