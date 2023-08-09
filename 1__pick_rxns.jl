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

# save them to the model directory under /mechanism
open(joinpath(outpath, "mechanism", "bimol.json"), "w") do f
    JSON.print(f, bimol_db, 2)
end

open(joinpath(outpath, "mechanism", "trimol.json"), "w") do f
    JSON.print(f, trimol_db, 2)
end

open(joinpath(outpath, "mechanism", "photolysis.json"), "w") do f
    JSON.print(f, photolysis_db, 2)
end


const kb = 1.380649e−23 # J/K

# number density
function M(T,P)
    press_pa = 100.0 * P  # 100 Pa / mbar
    # 1 Pa = 1 J / m3. We want cm^3 so convert:
    press_final = press_pa * 1.0e-6 # there are (10^2)^3 = 10^6 cm³/m³
    Mout = press_final/(kb*T)  # we now have a stand in for number density in molecules/cm³jk
    return Mout
end

# O2 and N2 based on typical relative abundance
O2(T,P) = 0.2095 * M(T,P)
N2(T,P) = 0.7809 * M(T,P)




T = 298.15
P = 996.0
d = M(T,P)

for i ∈ 1:length(bimol_db)
    println(i, "\t", bimol_db[i](T, P, d))
end


for i ∈ 1:length(trimol_db)
    println(i, "\t", trimol_db[i](T, P, d))
end



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
    for rxn ∈ photolysis_db
        reactants = [get_tex(r, df_species) for r ∈ rxn.reactants]
        products = [get_tex(p, df_species) for p ∈ rxn.products]

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




