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
readdir(rxn_databases)
photo_list = []
for f ∈ readdir(rxn_databases)
    if occursin("ph.d", f) || occursin("-ph", f) || occursin("photo", f)
        push!(photo_list, joinpath(rxn_databases, f))
    end
end

photo_list
@assert all(isfile.(photo_list))


parse_photolysis_d(photo_list[1])


outpath = "src/json-databases/photolysis"
if !ispath(outpath)
    mkpath(outpath)
end


for photolysis_file ∈ photo_list
    println(photolysis_file)
    fname = split(basename(photolysis_file),".")[1]*".json"
    rxns, rxns_failed = parse_photolysis_d(photolysis_file)
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





