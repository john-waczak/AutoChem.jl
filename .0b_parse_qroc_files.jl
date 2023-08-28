using AutoChem
using DelimitedFiles
using CSV, DataFrames
using JSON


# save the files
outpath = "src/species"
if !ispath(outpath)
    mkpath(outpath)
    mkpath(joinpath(outpath, "species"))
    mkpath(joinpath(outpath, "totals"))
    mkpath(joinpath(outpath, "ratios"))
end



qroc_files = "src/autochem-databases/qroc-scripts"
@assert ispath(qroc_files)


qroc_list = joinpath.(qroc_files, readdir(qroc_files))
@assert all(isfile.(qroc_list))

qroc_list

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


master_species_dfs = []

for qroc_f ∈ qroc_list
    println(qroc_f)
    fname = split(split(qroc_f, "/")[end], ".")[1]

    try
        species = vcat(process_species.(get_species_list(qroc_f, "specie.d"))...)
        push!(master_species_dfs, species)
        totals = vcat(process_total.(get_species_list(qroc_f, "total.ctl"))...)
        ratios = vcat(process_ratio.(get_species_list(qroc_f, "ratio.ctl"))...)


        CSV.write(joinpath(outpath, "species", fname * ".csv"), species)
        CSV.write(joinpath(outpath, "totals", fname * ".csv"), totals)
        CSV.write(joinpath(outpath, "ratios", fname * ".csv"), ratios)

    catch e
        println("\t", qroc_f, " failed!")
        println(e)
    end


end


master_species_list = unique!(vcat(master_species_dfs...))
CSV.write(joinpath(outpath, "species", "master_species_list.csv"), master_species_list)
