using AutoChem
using DelimitedFiles
using JSON

rxn_databases = "src/autochem-databases/rxn-databases"
@assert ispath(rxn_databases)


# 1. read in reaction databases and parse to JSON object

rxns, rxns_failed = parse_bimol_d(joinpath(rxn_databases, "v5_ac-bi.d"))



open("test.json", "w")  do f
    JSON.print(f, rxns)
end


JSON.parsefile("test.json")
