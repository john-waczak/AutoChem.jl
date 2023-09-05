using AutoChem
using DelimitedFiles
using JSON


rxn_databases = "src/autochem-databases/rxn-databases"


# -------------------------------------------------------------
# Bimolecular Reaction Databases
# -------------------------------------------------------------

bimol_list = []
for f ∈ readdir(rxn_databases)
    if occursin("bi", f) && !occursin("bulk", f)
        push!(bimol_list, joinpath(rxn_databases, f))
    end
end

bimol_list
@assert all(isfile.(bimol_list))


outpath = "src/json-databases/bimol"
if !ispath(outpath)
    mkpath(outpath)
end


for bimol_file ∈ bimol_list
    println(bimol_file)
    fname = split(basename(bimol_file),".")[1]*".json"
    rxns, rxns_failed = parse_bimol_d(bimol_file)
    open(joinpath(outpath, fname), "w") do f
        JSON.print(f, rxns, 2)
    end


    if !isempty(rxns_failed)
        println("\tfailed!")
        println("\t", rxns_failed)
    end

end


test_bi = joinpath(outpath, "basic-bi.json")
typeof(read_bimol(test_bi))


# -------------------------------------------------------------
# Trimolecular Reaction Databases
# -------------------------------------------------------------

readdir(rxn_databases)
trimol_list = []
for f ∈ readdir(rxn_databases)
    if occursin("tri", f)
        push!(trimol_list, joinpath(rxn_databases, f))
    end
end

trimol_list
@assert all(isfile.(trimol_list))


parse_trimol_d(trimol_list[1])

outpath = "src/json-databases/trimol"
if !ispath(outpath)
    mkpath(outpath)
end


for trimol_file ∈ trimol_list
    println(trimol_file)
    fname = split(basename(trimol_file),".")[1]*".json"
    rxns, rxns_failed = parse_trimol_d(trimol_file)
    open(joinpath(outpath, fname), "w") do f
        JSON.print(f, rxns, 2)
    end


    if !isempty(rxns_failed)
        println("\tfailed!")
        println("\t", rxns_failed)
    end

end

test_tri = joinpath(outpath, "basic-tri.json")
typeof(read_trimol(test_tri))


# -------------------------------------------------------------
# Photolysis Reaction Databases
# -------------------------------------------------------------
using CSV, DataFrames
df_lookup = CSV.read("./src/species/master_species_list.csv", DataFrame)


readdir(rxn_databases)
photo_list = []
for f ∈ readdir(rxn_databases)
    if occursin("ph.d", f) || occursin("-ph", f) || occursin("photo", f)
        push!(photo_list, joinpath(rxn_databases, f))
    end
end

photo_list
@assert all(isfile.(photo_list))


photo_list[1]

db, failed = parse_photolysis_d(photo_list[1])

db[1].autochem_files
db[1].crosssection_files
db[1].quantumyield_files
db[1].reactants
db[1].products


function append_σ_info!(db)
    for rxn ∈ db
        for reactant ∈ rxn.reactants
            if reactant ∈ df_lookup.varname && reactant != "Photon"
                # get idx of reactant
                idx = findfirst(df_lookup.varname .== reactant)
                if !ismissing(df_lookup[idx, "mpi-mainz-uvviz"])
                    println("σ: ", df_lookup.varname[idx], "\t", df_lookup[idx, "mpi-mainz-uvviz"])
                    if length(rxn.crosssection_files) < 2
                        rxn.crosssection_files[1] = df_lookup[idx, "mpi-mainz-uvviz"]*".csv"
                    else
                        push!(rxn.crosssection_files, df_lookup[idx, "mpi-mainz-uvviz"]*".csv")
                    end
                end
            end
        end
    end
end


outpath = "src/json-databases/photolysis"
if !ispath(outpath)
    mkpath(outpath)
end


for photolysis_file ∈ photo_list
    println(photolysis_file)
    fname = split(basename(photolysis_file),".")[1]*".json"
    txtname = split(basename(photolysis_file),".")[1]*".txt"

    rxns, rxns_failed = parse_photolysis_d(photolysis_file)
    append_σ_info!(rxns)

    open(joinpath(outpath, fname), "w") do f
        JSON.print(f, rxns, 2)
    end

    if !isempty(rxns_failed)
        println("\tfailed!")
        println("\t", rxns_failed)
    end


end


test_photolysis = joinpath(outpath, "cri-photo.json")
typeof(read_photolysis(test_photolysis))




