ENV["GKSwstype"] = 100

@info "Setting up Julia Environment..."
using Pkg
Pkg.activate(".")
Pkg.instantiate()
@info "\t...finished"

using AutoChem
using CSV, DataFrames
using JSON
using ArgParse
using ProgressMeter
using LaTeXStrings, YAML

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
            default = "high_primed"
        "--unc_ext"
            help = "Extension for uncertainty files."
            arg_type = String
            default = "_std"
        "--qroc"
            help = "Autochem qroc used to select relevant species and reaction databases."
            arg_type = String
            #default = "qroc-methane"
            default = "qroc-methane-ion-nagfor"
        "--model_name"
            help = "Name for the resulting model used in output paths"
            arg_type = String
            #default = "methane"
            default = "autochem-w-ions"
        "--time_step"
            help = "The time step used during integration of mechanism (in minutes)."
            arg_type = Float64
            default = 15.0
    end


    parsed_args = parse_args(ARGS, s; as_symbols=true)

    @assert ispath(parsed_args[:data_basepath])
    @assert ispath(joinpath(parsed_args[:data_basepath], "number_densities", parsed_args[:collection_id]))
    @assert isfile(joinpath("src", "species", "species", parsed_args[:qroc] * ".csv"))

    return parsed_args
end





@info "Parsing command line flags..."

parsed_args = parse_commandline()

data_basepath = parsed_args[:data_basepath]
model_name = parsed_args[:model_name]
qroc = parsed_args[:qroc]
collection_id = parsed_args[:collection_id]
unc_ext = parsed_args[:unc_ext]
Δt_step = parsed_args[:time_step]  # time step in minutes

model_path= "models/"

@info "Setting up file paths..."

if !ispath(joinpath(model_path, model_name, "runs", collection_id))
    mkpath(joinpath(model_path, model_name, "runs", collection_id))
    mkpath(joinpath(model_path, model_name, "runs", collection_id, "mechanism"))
    mkpath(joinpath(model_path, model_name, "runs", collection_id, "figures"))
    mkpath(joinpath(model_path, model_name, "runs", collection_id, "4d-var"))
    mkpath(joinpath(model_path, model_name, "runs", collection_id, "EKF"))
end

outpath = joinpath(model_path, model_name, "runs", collection_id)

docs_path = joinpath(model_path, model_name, "docs")
if !ispath(docs_path)
    mkpath(docs_path)
end



# 1. Create Species List from list of qroc files
function join_species(qroc)
    df_species = CSV.File(joinpath("src", "species", "species", qroc * ".csv")) |> DataFrame
    return df_species
end


function join_totals(qroc)
    df_totals = CSV.File(joinpath("src", "species", "totals", qroc * ".csv")) |> DataFrame
    return df_totals
end


function join_ratios(qroc)
    df_ratios = CSV.File(joinpath("src", "species", "ratios", qroc * ".csv")) |> DataFrame
end



# generate the lists
@info "Joining species lists..."
df_species = join_species(qroc)
df_totals = join_totals(qroc)
df_ratios = join_ratios(qroc)


# drop species in ignore list
@info "Dropping species to ignore"
ignore_list = [
    "HClS",
    "H2OS",
    "HONO2S",
    "CO3-H2O",
    "CO3-H2OH2O",
    "NO2-H2O",
    "NO3-H2O",
    "NO3-H2OH2O",
]



# remove ignore species
@info "Removing ignored species from list..."
df_species = df_species[[!(name ∈ ignore_list) for name ∈ df_species.varname],:]


# update species df to include species index:
@info "Adding index to species df..."
df_species.idx_species = [i for i ∈ 1:nrow(df_species)]


# save them to the mechanism directory
@info "Saving lists..."
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


@info "Fetching reaction databases"

paths = [joinpath("src/autochem-databases/qroc-scripts", qroc)]
@assert ispath(paths)

bimol_dbs, trimol_dbs, photo_dbs = get_rxn_databases(paths);

@assert all(isfile.(bimol_dbs))
@assert all(isfile.(trimol_dbs))
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
@info "Picking reactions from rxn dbs"
bimol_db = pick_rxns(bimol_dbs, df_species);
trimol_db = pick_rxns(trimol_dbs, df_species, type=:trimol);
photolysis_db = pick_rxns(photo_dbs, df_species, type=:photolysis);


# create new versions of the databases with integers instead of variable names
@info "Converting variable names to indices"
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
@info "Saving finalized databases"

open(joinpath(outpath, "mechanism", "bimol.json"), "w") do f
    JSON.print(f, bimol_db_out, 2)
end

open(joinpath(outpath, "mechanism", "trimol.json"), "w") do f
    JSON.print(f, trimol_db_out, 2)
end

open(joinpath(outpath, "mechanism", "photolysis.json"), "w") do f
    JSON.print(f, photo_db_out, 2)
end


@info "Testing reaction rate coefficient calculations"

T = 298.15
P = 996.0
d = M(T,P)

for i ∈ 1:length(bimol_db_out)
    if bimol_db_out[i](T, P, d) == 0.0
        println(i)
        println(bimol_db_out[i])
    end
end

for i ∈ 1:length(trimol_db_out)
    trimol_db[i](T, P, d)
end


# create auto-documentation for databases

bimol_path = joinpath(docs_path, "bimol.qmd")
trimol_path = joinpath(docs_path, "trimol.qmd")
photo_path = joinpath(docs_path, "photolysis.qmd")



rxn = bimol_db_out[1]
get_reaction_tex(bimol_db_out[23])


bimol_db_out[6](T,P,d)
get_reaction_tex(bimol_db_out[6])

bimol_db_out[6]



get_reaction_tex(trimol_db_out[23])


@info "Creating _quarto.yml specification for docs"

quarto_dict = Dict(
    "project" => Dict(
        "type" => "website"
    ),
    "website" => Dict(
        "title" => "$(model_name)",
        "navbar" => Dict(
            "background" => "primary",
            "search" => true,
            "left" => [
                Dict(
                    "href" => "index.qmd",
                    "text" => "Overview"
                ),
                Dict(
                    "href" => "bimol.qmd",
                    "text" => "Bimolecular Reactions"
                ),
                Dict(
                    "href" => "trimol.qmd",
                    "text" => "Trimolecular Reactions"
                ),
                Dict(
                    "href" => "photolysis.qmd",
                    "text" => "Photolysis Reactions"
                ),
            ]
        ),
    ),
    "format" => Dict(
        "html" => Dict(
            "theme" => "cosmo",
            # "css" => "styles.css",
            "toc" => true,
        )
    ),
    "footnotes" => "margin",
    "references" => "margin",
)

YAML.write_file(joinpath(docs_path, "_quarto.yml"), quarto_dict)


@info "Creating index file"

open(joinpath(docs_path, "index.qmd"), "w") do f
    println(f, "This is is the homepage for $(model_name)\n\n")
    println(f, "| Index | Species Name | Variable Name | Is Integrated? |")
    println(f, "|:-:|:----:|:----:|:-:|")

    for row ∈ eachrow(df_species)
        is_int = true
        if row.is_integrated == 2
            is_int = false
        end

        out = "| $(row.idx_species) | \$\$ \\mathrm{$(row.printname)} \$\$ | $(row.varname) | $(is_int) |"
        println(f, out)
    end

    println(f, ": Model species list {.hover .bordered .striped}")
end


@info "Writing equations docs..."

open(bimol_path, "w") do f
    println(f, "| # | Bimolecular Reaction | Reaction Rate Coeff |")
    println(f, "|:-:|:-------:|:-------:|")

    for i ∈ 1:length(bimol_db_out)
        rxn = bimol_db_out[i]
        rrate = get_reaction_tex(rxn)
        out = "| $(i) | " * get_tex(rxn, df_species) * " | " * rrate * " |"

        println(f, out)
    end
    println(f, ": Bimolecular reaction definitions {.hover .bordered .striped}")
end

open(trimol_path, "w") do f
    println(f, "| # | Trimolecular Reaction | Reaction Rate Coeff |")
    println(f, "|:-:|:-------:|:--------:|")

    for i ∈ 1:length(trimol_db_out)
        rxn = trimol_db_out[i]
        rrate = get_reaction_tex(rxn)
        out = "| $(i) | " * get_tex(rxn, df_species) * " | " * rrate * " |"
        println(f, out)
    end
    println(f, ": Trimolecular reaction definitions {.hover .bordered .striped}")
end

open(photo_path, "w") do f
    println(f, "| # | Photolysis Reaction |")
    println(f, "|:-:|:--------------------:|")

    for i ∈ 1:length(photo_db_out)
        rxn = photo_db_out[i]
        out = "| $(i) | " * get_tex(rxn, df_species) * " |"

        println(f, out)
    end
    println(f, ": Photolysis reaction definitions {.hover .bordered .striped}")
end

