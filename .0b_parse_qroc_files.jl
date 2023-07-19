using AutoChem
using DelimitedFiles
using CSV, DataFrames
using JSON

qroc_files = "src/autochem-databases/qroc-scripts"
@assert ispath(qroc_files)

qroc_list = joinpath.(qroc_files, readdir(qroc_files))



function get_end_of_block(lines, idx_start)
    for i ∈ idx_start+1:length(lines)
        if occursin("EOF", lines[i])
            return i
        end
    end
    return -1
end


function get_species_list(path, section)
    qroc_lines = readlines(path)

    species_idxs = findall([occursin(section, line) && occursin("<<EOF", line) for line ∈ qroc_lines])

    species_slices = [i+2:get_end_of_block(qroc_lines, i)-1 for i ∈ species_idxs]

    return [qroc_lines[slice] for slice ∈ species_slices]
end



function process_species(species_list)
    table = (; is_integrated=Int[], varname=String[], printname=String[])

    for i ∈ 1:length(species_list)
        res = [strip(val) for val ∈ split(species_list[i], "'") if strip(val) != ""]

        push!(table.is_integrated, parse(Int, res[1]))
        push!(table.varname, res[2])
        push!(table.printname, res[3])
    end

    df_spec = table |> DataFrame

    return df_spec
end



function process_total(total_list)
    table = (; species_in_total=Int[], name=String[])

    for i ∈ 1:length(total_list)
        res = [strip(val) for val ∈ split(total_list[i], "'") if strip(val) != ""]

        push!(table.species_in_total, parse(Int, res[1]))
        push!(table.name, res[2])
    end

    df_tot = table |> DataFrame

    return df_tot
end



function process_ratio(ratio_list)
    table = (; amount_a=Int[], species_a=String[], amount_b=Int[], species_b=String[])

    for i ∈ 1:length(ratio_list)
        res = [strip(val) for val ∈ split(ratio_list[i], "'") if strip(val) != ""]
        if !isempty(res)
            push!(table.amount_a, parse(Int, res[1]))
            push!(table.species_a, res[2])
            push!(table.amount_b, parse(Int, res[3]))
            push!(table.species_b, res[4])
        end
    end

    df_ratio = table |> DataFrame

    return df_ratio
end


df_species = DataFrame()
df_totals = DataFrame()
df_ratios = DataFrame()

for qroc_f ∈ qroc_list
    println(qroc_f)
    try
        species = get_species_list(qroc_f, "specie.d")
        totals = get_species_list(qroc_f, "total.ctl")
        ratios = get_species_list(qroc_f, "ratio.ctl")


        for specie ∈ species
            df_species = vcat(df_species, process_species(specie))
        end

        for total ∈ totals
            df_totals = vcat(df_totals, process_total(total))
        end

        for ratio ∈ ratios
            df_ratios = vcat(df_ratios, process_ratio(ratio))
        end
    catch e
        println(qroc_f, " failed!")
        println(e)
    end
end

# get only the unique rows from the list of names
idx_unique = [findfirst(x->x==name, df_species.varname) for name ∈ unique(df_species.varname)]
df_species = df_species[idx_unique, :]

idx_unique = [findfirst(x->x==name, df_totals.name) for name ∈ unique(df_totals.name)]
df_totals = df_totals[idx_unique, :]

idx_unique = [findfirst(x->x==name, df_ratios.species_a) for name ∈ unique(df_ratios.species_a)]
df_ratios = df_ratios[idx_unique, :]


# save the files
outpath = "src/json-databases/species"
if !ispath(outpath)
    mkpath(outpath)
end

CSV.write(joinpath(outpath, "species.csv"), df_species)
CSV.write(joinpath(outpath, "totals.csv"), df_species)
CSV.write(joinpath(outpath, "ratios.csv"), df_species)

