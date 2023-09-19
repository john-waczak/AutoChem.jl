using AutoChem
using CSV, DataFrames
using JSON



# set up model output directory
model_name = "autochem-w-ions"
#model_name = "qroc-methane-intel"



model_path= "models/"
collection_id = "empty"

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
    # "qroc-red-NHMC-nagfor",
    "qroc-methane-ion-nagfor"
    #"qroc-methane-intel"
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


# drop species in ignore list
ignore_list = [
    "HClS",
    "H2OS",
    "HONO2S",
]
df_species = df_species[[!(name ∈ ignore_list) for name ∈ df_species.varname],:]

# update species df to include species index:
df_species.idx_species = [i for i ∈ 1:nrow(df_species)]



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
@assert all(isfile.(bimol_dbs))

trimol_dbs
@assert all(isfile.(trimol_dbs))

photo_dbs
@assert all(isfile.(photo_dbs))


function rxn_in_list(rxn, rxn_list)
    for rxnᵢ ∈ rxn_list
        if rxnᵢ.reactants == rxn.reactants && rxnᵢ.products == rxn.products
            return true
        end
    end
    return false
end


function reactants_and_products_in_list(rxn, df_species)
    rxnspecs = vcat(rxn.reactants, rxn.products)
    all([e ∈ df_species.varname for e ∈ rxnspecs])
end


function find_first_rxn(rxn, rxn_list)
    for i ∈ 1:length(rxn_list)
        if rxn.reactants == rxn_list[i].reactants && rxn.products == rxn_list[i].products
            return i
        end
    end
    return -1
end



function pick_rxns(dbs, df_species; type=:bimol)
    type_dict = Dict(
        :bimol => (; type=BimolecularReaction, readdb=read_bimol),
        :trimol => (; type=TrimolecularReaction, readdb=read_trimol),
        :photolysis => (; type=PhotolysisReaction, readdb=read_photolysis),
    )

    @assert type ∈ keys(type_dict)

    rxns = type_dict[type].type[]

    for db ∈ dbs
        for rxn ∈ type_dict[type].readdb(db)
            if !("?" ∈ rxn.reactants) && !("?" ∈ rxn.products)
                if reactants_and_products_in_list(rxn, df_species)
                    if rxn_in_list(rxn, rxns)
                        i = find_first_rxn(rxn, rxns)
                        rxns[i] = rxn
                    else
                        push!(rxns, rxn)
                    end
                end
            end
        end
    end

    return rxns
end

# generate the picked reactions:
bimol_db = pick_rxns(bimol_dbs, df_species);
trimol_db = pick_rxns(trimol_dbs, df_species, type=:trimol);
photolysis_db = pick_rxns(photo_dbs, df_species, type=:photolysis);

photolysis_db

# create new versions of the databases with integers instead of variable names
function get_species_index(species, df_species)
    idx_row = findfirst(df_species.varname .== species)
    return df_species.idx_species[idx_row]
end


function convert_rxn_vars_to_idxs(rxn::BimolecularReaction, df_species)
    reactants = [get_species_index(r, df_species) for r ∈ rxn.reactants]
    products = [get_species_index(p, df_species) for p ∈ rxn.products]

    return BimolecularReaction(
        rxn.idx,
        rxn.source,
        reactants,
        products,
        rxn.prod_stoich,
        rxn.a1,
        rxn.a2,
        rxn.a3,
        rxn.a4,
        rxn.a5,
        rxn.contains_OH,
        rxn.contains_HONO2,
        rxn.contains_CO,
        rxn.all_reactants_HO2
    )
end

function convert_rxn_vars_to_idxs(rxn::TrimolecularReaction, df_species)
    reactants = [get_species_index(r, df_species) for r ∈ rxn.reactants]
    products = [get_species_index(p, df_species) for p ∈ rxn.products]

    return TrimolecularReaction(
        rxn.idx,
        rxn.source,
        reactants,
        products,
        rxn.prod_stoich,
        rxn.a1,
        rxn.a2,
        rxn.a3,
        rxn.a4,
        rxn.a5,
        rxn.b1,
        rxn.b2,
        rxn.b3,
        rxn.b4,
        rxn.b5,
        rxn.c1,
        rxn.c2,
        rxn.c3,
        rxn.c4,
        rxn.c5,
        rxn.contains_m,
        rxn.contains_ClCO,
        rxn.contains_ClCO3,
        rxn.contains_COCl2,
        rxn.contains_APINENE,
        rxn.contains_BPINENE,
        rxn.contains_ClCO_and_ClCO3,
        rxn.contains_ClCO_and_COCl2,
        rxn.contains_only_low
)
end

function convert_rxn_vars_to_idxs(rxn::PhotolysisReaction, df_species)
    reactants = [get_species_index(r, df_species) for r ∈ rxn.reactants]
    products = [get_species_index(p, df_species) for p ∈ rxn.products]

    return PhotolysisReaction(
        rxn.idx,
        rxn.source,
        reactants,
        products,
        rxn.prod_stoich,
        rxn.autochem_files,
        rxn.crosssection_files,
        rxn.quantumyield_files
    )
end

bimol_db_out = [convert_rxn_vars_to_idxs(rxn, df_species) for rxn ∈ bimol_db]
trimol_db_out = [convert_rxn_vars_to_idxs(rxn, df_species) for rxn ∈ trimol_db]
photo_db_out = [convert_rxn_vars_to_idxs(rxn, df_species) for rxn ∈ photolysis_db]

# save them to the model directory under /mechanism
open(joinpath(outpath, "mechanism", "bimol.json"), "w") do f
    JSON.print(f, bimol_db_out, 2)
end

open(joinpath(outpath, "mechanism", "trimol.json"), "w") do f
    JSON.print(f, trimol_db_out, 2)
end

open(joinpath(outpath, "mechanism", "photolysis.json"), "w") do f
    JSON.print(f, photo_db_out, 2)
end


T = 298.15
P = 996.0
d = M(T,P)

for i ∈ 1:length(bimol_db_out)
    println(i, "\t", bimol_db_out[i](T, P, d))
end


for i ∈ 1:length(trimol_db_out)
    println(i, "\t", trimol_db[i](T, P, d))
end



# create auto-documentation for databases

using LaTeXStrings

function get_tex(name::String, df_species::DataFrame)
    idx = findfirst(df_species.varname .== name)
    tex = df_species.printname[idx]
    return tex
end


docs_path = joinpath(outpath, "docs")
if !ispath(docs_path)
    mkpath(docs_path)
end

photo_path = joinpath(docs_path, "photolysis.qmd")


open(photo_path, "w") do f
    for rxn ∈ photo_db_out
        # reactants = [get_tex(r, df_species) for r ∈ rxn.reactants]
        # products = [get_tex(p, df_species) for p ∈ rxn.products]

        reactants = df_species.printname[rxn.reactants]
        products = df_species.printname[rxn.products]


        reactants = [ ("h\\nu" == r) ? "h\\nu" : "\\mathrm{$r}" for r ∈ reactants]
        products = ["\\mathrm{$p}" for p ∈ products]
        pstoich = Int.(rxn.prod_stoich)

        for i ∈ 1:length(products)
            if pstoich[i] > 1
                products[i] = "$(pstoich[i])" * products[i]
            end
        end

        reactants = join([r for r ∈ reactants], " + ")
        products = join([p for p ∈ products], " + ")

        out = """
    \\begin{equation}
        $reactants \\longrightarrow $products
    \\end{equation}
    """
        println(f, out)


    end
end
