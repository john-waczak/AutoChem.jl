using HTTP
using Gumbo
using Cascadia
using DelimitedFiles
using CSV, DataFrames


# 0. Set up paths

outpath_base = "../../assets/mpi-mainz-uvviz/raw"
if !isdir(outpath_base)
    mkpath(outpath_base)
end

outpath_crosssections = joinpath(outpath_base, "cross-sections")
if !isdir(outpath_crosssections)
    mkpath(outpath_crosssections)
end

outpath_quantumyields = joinpath(outpath_base, "quantum-yields")
if !isdir(outpath_quantumyields)
    mkpath(outpath_quantumyields)
end



base_url = "https://uv-vis-spectral-atlas-mainz.org/uvvis/"
cross_sections_url = joinpath(base_url, "cross_sections/")
quantum_yields_url = joinpath(base_url, "quantum_yields/")








# --------------------------------------------------------------------------------- #
# 1. Cross Sections
# --------------------------------------------------------------------------------- #

function getcategories(text)
    h = parsehtml(text)
    body = h.root[2]
    s = eachmatch(Selector("table"), body)

    table = s[4][1]
    nrows = length(table.children)

    categories = String[]
    for i ∈ 1:nrows
        push!(categories, table[i][1][1][1].text)
    end

    categories_nospace = replace.(categories, " "=>"%20")
    return categories, categories_nospace
end


function hassubcategories(url)
    r = HTTP.get(url)
    text = String(r.body)
    return occursin("Subcategories:", text), text
end

function getsubcategories(text)
    h = parsehtml(text)
    body = h.root[2]
    s = eachmatch(Selector("table"), body)

    table = s[4][1]
    nrows = length(table.children)

    subcategories = String[]

    for i ∈ 1:nrows
        push!(subcategories, table[i][1][1][1].text)
    end

    subcategories_nospace = replace.(subcategories, " "=>"%20")
    return subcategories, subcategories_nospace
end



function hascasnumber(url)
    r = HTTP.get(url)
    text = String(r.body)
    return occursin("CAS RN", text), text
end

function getspecies(text)
    h = parsehtml(text)
    body = h.root[2]
    s = eachmatch(Selector("table"), body)
    table = s[4][1]

    nrows = length(table.children)

    formulae = String[]
    for i ∈ 2:nrows
        push!(formulae, table[i][1][1][1].text)
    end

    return formulae
end


# loop through Formula and generate a csv with links to each dataset
function hasdataset(url)
    r = HTTP.get(url)
    text = String(r.body)
    return occursin("Data Sets:", text), text
end


function getdatasets(text)
    h = parsehtml(text)
    body = h.root[2]
    s = eachmatch(Selector("table"), body)

    if occursin("graphical representations not available", text)
        table = s[5][1]
        nrows = length(table.children)

        datasets = String[]
        for i ∈ 2:nrows
            new_url = table[i][1][1].attributes["href"]
            if split(new_url, "/")[1] != "https"
                new_url = joinpath(base_url, new_url)
            end
            push!(datasets,new_url)
        end
    else
        table = s[6][1]
        nrows = length(table.children)

        datasets = String[]
        for i ∈ 2:nrows
            new_url = table[i][1][1].attributes["href"]
            if split(new_url, "/")[1] != "https"
                new_url = joinpath(base_url, new_url)
            end
            push!(datasets,new_url)
        end
    end

    # check to make sure we start with the base_url
    # for i ∈ 1:length(datasets)
    #     if split(datasets[i], "/")[1] != "https"
    #         datasets[i] = joinpath(base_url, datasets[i])
    #     end

    # end

    return datasets
end


function getdatainfo(url)
    r = HTTP.get(url)
    text = String(r.body)
    h = parsehtml(text)
    body = h.root[2]
    s = eachmatch(Selector("table"), body)
    length(s)

    table = s[4][1]

    datafile = table[1][2][1].attributes["href"]
    datafileurl = replace("https://uv-vis-spectral-atlas-mainz.org/"*datafile, " "=>"%20")

    out_dict = Dict()


    out_dict[:download_url] = datafileurl
    out_dict[:name] = table[2][2][1].text
    out_dict[:formula] =  table[3][2][1].text


    splitname = split(replace(datafile, "/uvvis_data/cross_sections"=>"cross-sections"), "/")
    fname = joinpath(splitname[1:end-1]..., out_dict[:formula], splitname[end])
    fname = split(fname, ".")[1]*".csv"
    out_dict[:fname] = fname

    out_dict[:author_year] = table[4][2][1].text

    temp = lstrip(table[5][2][1].text)
    idxs = findfirst("K", temp)
    temp = temp[1:idxs[1]-1]
    if occursin("-", temp)
        out_dict[:T1] = parse(Float64, split(temp, "-")[1])
        out_dict[:T2] = parse(Float64, split(temp, "-")[2])
    elseif occursin(",", temp)
        out_dict[:T1] = parse(Float64, split(temp, ",")[1])
        out_dict[:T2] = parse(Float64, split(temp, ",")[2])
    elseif occursin("pm", temp)
        μₜ = parse(Float64, split(temp, "pm")[1])
        out_dict[:T1] = μₜ
        out_dict[:T2] = NaN
    else
        out_dict[:T1] = parse(Float64, temp)
        out_dict[:T2] = NaN
    end


    λ = lstrip(table[6][2][1].text)
    idxs = findfirst("nm", λ)
    λ = λ[1:idxs[1]-1]

    if occursin("-", λ)
        out_dict[:λ1] = parse(Float64, split(λ, "-")[1])
        out_dict[:λ2] = parse(Float64, split(λ, "-")[2])
    elseif occursin(",", λ)
        out_dict[:λ1] = parse(Float64, split(λ, ",")[1])
        out_dict[:λ2] = parse(Float64, split(λ, ",")[2])
    elseif occursin("pm", λ)
        μλ = parse(Float64, split(λ, "pm")[1])
        out_dict[:λ1] = μλ
        out_dict[:λ2] = NaN
    else
        out_dict[:λ1] = parse(Float64, λ)
        out_dict[:λ2] = NaN
    end

    try
        nitems = length(table[7][2].children)
        out_dict[:doi] = table[7][2][nitems].attributes["href"]
    catch
        out_dict[:doi] = ""
    end

    try
        out_dict[:comments] = table[8][2][1].text
    catch
        out_dict[:comments] = ""
    end

    return out_dict
end


getdatainfo("https://uv-vis-spectral-atlas-mainz.org/uvvis/cross_sections/Alcohols/H2N(CH2)2OH_Kuo(2020)_295.9-299.7K_185-250nm.txt")


# # this doesn't work for everything yet... everything seems to have it's own format...
# function getdata(url, fname, outpath)
#     tmp = download(url)

#      df = CSV.File(tmp) |> DataFrame

#     # figure out which columns actually have data
#     datacols = [eltype(col) != Missing for col ∈ eachcol(df)]

#     df = df[!, datacols]

#     # alternative strategy for one line files
#     if size(df) == (0,0)
#         dlm = readdlm(tmp)
#         df = DataFrame()
#         df.λ = [dlm[1,1]]
#         df.σ = [dlm[1,2]]
#     end

#     # check for weird edge case
#     if typeof(df[1,1]) <: String
#         if occursin("---", df[1,1])
#             df = df[2:end, :]
#         end
#     end


#     return_df = DataFrame()

#     return_df.λ = df[!, 1]
#     return_df.σ = df[!, 2]

#     if ncol(df) > 2
#         return_df.uncertainty = df[!,3]
#         if any(ismissing.(df[!,3])) && ncol(df) >= 4
#             # maybe it got misplaced to the left?
#             idxs = ismissing.(df[!,3])
#             return_df.uncertainty[idxs] = df[idxs, 4]
#         end
#     end

#     CSV.write(joinpath(outpath, fname), return_df)
# end



# generate list of categories
r = HTTP.get(cross_sections_url)
text = String(r.body)
categories, categories_nospace = getcategories(text)

# create subdirectories
for cat ∈ categories
    catpath = joinpath(outpath_crosssections, cat)
    if !isdir(catpath)
        mkpath(catpath)
    end
end


# loop through each category and look for subcategories
urllist = String[]
failure_list = String[]



for i ∈ 1:length(categories)
    url = joinpath(cross_sections_url, categories_nospace[i])
    hassubcat, text = hassubcategories(url);
    if hassubcat
        println(categories[i], " has subcategories!")
        subcats, subcats_nospace = getsubcategories(text)

        for j ∈ 1:length(subcats)
            subcatpath = joinpath(outpath_crosssections, categories[i], subcats[j])
            if !isdir(subcatpath)
                mkpath(subcatpath)
            end

            subcaturl = joinpath(url, subcats_nospace[j])
            hascas, text = hascasnumber(subcaturl);
            formulae = getspecies(text)

            for formula ∈ formulae
                formulapath = joinpath(outpath_crosssections, categories[i], subcats[j], formula)
                if !isdir(formulapath)
                    mkpath(formulapath)
                end

                formulaurl = joinpath(subcaturl, formula*".spc")
                push!(urllist, formulaurl)

                # for each formula, grab all datasets and generate csv with links.
                try
                    hasdata, text = hasdataset(formulaurl)
                    if hasdata
                        datasets = getdatasets(text)

                        summary_dict = Dict(
                            :T1 => Float64[],
                            :T2 => Float64[],
                            :λ1 => Float64[],
                            :λ2 => Float64[],
                            :author_year => String[],
                            :download_url => String[],
                            :comments => String[],
                            :formula => String[],
                            :name => String[],
                            :doi=> String[],
                            :fname => String[],
                        )

                        for dataset ∈ datasets
                            println("\t working on $(dataset)")
                            try
                                info = getdatainfo(replace(dataset, " "=>"%20"))

                                push!(summary_dict[:T1], info[:T1])
                                push!(summary_dict[:T2], info[:T2])
                                push!(summary_dict[:λ1], info[:λ1])
                                push!(summary_dict[:λ2], info[:λ2])
                                push!(summary_dict[:author_year], info[:author_year])
                                push!(summary_dict[:download_url], info[:download_url])
                                push!(summary_dict[:comments], info[:comments])
                                push!(summary_dict[:formula], info[:formula])
                                push!(summary_dict[:name], info[:name])
                                push!(summary_dict[:doi], info[:doi])
                                push!(summary_dict[:fname], info[:fname])

                                # getdata(info[:download_url], info[:fname], "../")

                            catch e
                                println("FAILURE: $(dataset)")
                                println(e)
                                push!(failure_list, dataset)
                            end
                        end

                        CSV.write(joinpath(formulapath, "summary.csv"), summary_dict)
                    end
                catch
                    println("FAILURE: Couldn't read $(formulaurl)")
                end

            end

        end
    else
        hascas, text = hascasnumber(url);
        formulae = getspecies(text)

        for formula ∈ formulae
            formulapath = joinpath(outpath_crosssections, categories[i], formula)
            if !isdir(formulapath)
                mkpath(formulapath)
            end

            formulaurl = joinpath(url, formula*".spc")
            push!(urllist, formulaurl)

            try
                # for each formula, grab all datasets and generate csv with links.
                hasdata, text = hasdataset(formulaurl)
                if hasdata
                    datasets = getdatasets(text)
                    summary_dict = Dict(
                        :T1 => Float64[],
                        :T2 => Float64[],
                        :λ1 => Float64[],
                        :λ2 => Float64[],
                        :author_year => String[],
                        :download_url => String[],
                        :comments => String[],
                        :formula => String[],
                        :name => String[],
                        :doi=> String[],
                        :fname => String[],
                    )


                    for dataset ∈ datasets
                        println("\t working on $(dataset)")
                        try
                            info = getdatainfo(replace(dataset, " "=>"%20"))


                            push!(summary_dict[:T1], info[:T1])
                            push!(summary_dict[:T2], info[:T2])
                            push!(summary_dict[:λ1], info[:λ1])
                            push!(summary_dict[:λ2], info[:λ2])
                            push!(summary_dict[:author_year], info[:author_year])
                            push!(summary_dict[:download_url], info[:download_url])
                            push!(summary_dict[:comments], info[:comments])
                            push!(summary_dict[:formula], info[:formula])
                            push!(summary_dict[:name], info[:name])
                            push!(summary_dict[:doi], info[:doi])
                            push!(summary_dict[:fname], info[:fname])

                            # getdata(info[:download_url], info[:fname], "../")

                        catch e
                            println("FAILURE: $(dataset)")
                            println(e)
                            push!(failure_list, dataset)
                        end
                    end

                    CSV.write(joinpath(formulapath, "summary.csv"), summary_dict)
                end

            catch
                println("FAILURE: Couldn't read $(formulaurl)")
            end

        end
    end
end

# save file with urls
writedlm(joinpath(outpath_base, "cross_section_urls.csv"), urllist, ",")
writedlm(joinpath(outpath_base, "failed_urls.csv"), failure_list, ",")



# --------------------------------------------------------------------------------- #
# 2. Quantum Yields
# --------------------------------------------------------------------------------- #

function hasreaction(url)
    r = HTTP.get(url)
    text = String(r.body)
    return occursin("Reaction:", text), text
end

function getreactions(text)
    h = parsehtml(text)
    body = h.root[2]
    s = eachmatch(Selector("table"), body)
    table = s[4][1]

    nrows = length(table.children)

    rxns = String[]
    rxns_url = String[]
    for i ∈ 2:nrows
        push!(rxns, table[i][1][1][1].text)
        push!(rxns_url, table[i][1][1].attributes["href"])
    end

    rxns_url = replace.(rxns_url, " "=>"%20", "{"=>"%7B", "}"=>"%7D")

    return rxns, rxns_url
end

function getrxninfo(url)
    r = HTTP.get(url)
    text = String(r.body)
    h = parsehtml(text)
    body = h.root[2]
    s = eachmatch(Selector("table"), body)
    length(s)

    table = s[4][1]


    datafile = table[1][2][1].attributes["href"]
    datafileurl = replace("https://uv-vis-spectral-atlas-mainz.org/"*datafile, " "=>"%20", "{"=>"%7B", "}"=>"%7D")

    out_dict = Dict()


    out_dict[:download_url] = datafileurl
    out_dict[:name] = table[2][2][1].text
    out_dict[:reaction] =  table[3][2][1].text

    splitname = split(replace(datafile, "/uvvis_data/quantum_yields"=>"quantum-yields"), "/")
    fname = joinpath(splitname[1:end-1]..., out_dict[:reaction], splitname[end])
    fname = split(fname, ".")[1]*".csv"
    out_dict[:fname] = fname

    out_dict[:author_year] = table[4][2][1].text


    temp = lstrip(table[5][2][1].text)
    idxs = findfirst("K", temp)
    temp = temp[1:idxs[1]-1]
    if occursin("-", temp)
        out_dict[:T1] = parse(Float64, split(temp, "-")[1])
        out_dict[:T2] = parse(Float64, split(temp, "-")[2])
    elseif occursin(",", temp)
        out_dict[:T1] = parse(Float64, split(temp, ",")[1])
        out_dict[:T2] = parse(Float64, split(temp, ",")[2])
    elseif occursin("pm", temp)
        μₜ = parse(Float64, split(temp, "pm")[1])
        out_dict[:T1] = μₜ
        out_dict[:T2] = NaN
    else
        out_dict[:T1] = parse(Float64, temp)
        out_dict[:T2] = NaN
    end


    λ = lstrip(table[6][2][1].text)
    idxs = findfirst("nm", λ)
    λ = λ[1:idxs[1]-1]

    if occursin("-", λ)
        out_dict[:λ1] = parse(Float64, split(λ, "-")[1])
        out_dict[:λ2] = parse(Float64, split(λ, "-")[2])
    elseif occursin(",", λ)
        out_dict[:λ1] = parse(Float64, split(λ, ",")[1])
        out_dict[:λ2] = parse(Float64, split(λ, ",")[2])
    elseif occursin("pm", λ)
        μλ = parse(Float64, split(λ, "pm")[1])
        out_dict[:λ1] = μλ
        out_dict[:λ2] = NaN
    else
        out_dict[:λ1] = parse(Float64, λ)
        out_dict[:λ2] = NaN
    end

    try
        nitems = length(table[7][2].children)
        out_dict[:doi] = table[7][2][nitems].attributes["href"]
    catch
        out_dict[:doi] = ""
    end

    try
        out_dict[:comments] = table[8][2][1].text
    catch
        out_dict[:comments] = ""
    end

    return out_dict
end



# generate list of categories
r = HTTP.get(quantum_yields_url)
text = String(r.body)
categories, categories_nospace = getcategories(text)

# create subdirectories
for cat ∈ categories
    catpath = joinpath(outpath_quantumyields, cat)
    if !isdir(catpath)
        mkpath(catpath)
    end
end

# loop through each category and look for subcategories
urllist = String[]
failure_list = String[]

for i ∈ 1:length(categories)
    url = joinpath(quantum_yields_url, categories_nospace[i])
    hassubcat, text = hassubcategories(url);
    if hassubcat
        println(categories[i], " has subcategories!")
        subcats, subcats_nospace = getsubcategories(text)

        for j ∈ 1:length(subcats)
            subcatpath = joinpath(outpath_quantumyields, categories[i], subcats[j])
            if !isdir(subcatpath)
                mkpath(subcatpath)
            end

            subcaturl = joinpath(url, subcats_nospace[j])


            hasrxn, text = hasreaction(subcaturl);
            rxns, rxns_url = getreactions(text)

            for k ∈ 1:length(rxns)
                rxnpath = joinpath(outpath_quantumyields, categories[i], subcats[j], rxns[k])
                if !isdir(rxnpath)
                    mkpath(rxnpath)
                end

                rxnurl = joinpath(base_url, rxns_url[k])
                push!(urllist, rxnurl)

                try
                    hasdata, text = hasdataset(rxnurl)
                    if hasdata
                        datasets = getdatasets(text)
                        datasets = replace.(datasets, " "=>"%20", "{"=>"%7B", "}"=>"%7D")

                        summary_dict = Dict(
                            :T1 => Float64[],
                            :T2 => Float64[],
                            :λ1 => Float64[],
                            :λ2 => Float64[],
                            :author_year => String[],
                            :download_url => String[],
                            :comments => String[],
                            :reaction => String[],
                            :name => String[],
                            :doi=> String[],
                            :fname => String[],
                        )

                        for dataset ∈ datasets
                            println("\t working on $(dataset)")
                            try
                                info = getrxninfo(dataset)

                                push!(summary_dict[:T1], info[:T1])
                                push!(summary_dict[:T2], info[:T2])
                                push!(summary_dict[:λ1], info[:λ1])
                                push!(summary_dict[:λ2], info[:λ2])
                                push!(summary_dict[:author_year], info[:author_year])
                                push!(summary_dict[:download_url], info[:download_url])
                                push!(summary_dict[:comments], info[:comments])
                                push!(summary_dict[:reaction], info[:reaction])
                                push!(summary_dict[:name], info[:name])
                                push!(summary_dict[:doi], info[:doi])
                                push!(summary_dict[:fname], info[:fname])

                                # getdata(info[:download_url], info[:fname], "../")

                            catch e
                                println("FAILURE: $(dataset)")
                                println(e)
                                push!(failure_list, dataset)
                            end
                        end

                        CSV.write(joinpath(rxnpath, "summary.csv"), summary_dict)
                    end
                catch
                    println("FAILURE: Couldn't read $(rxnurl)")
                end

            end

        end
    else

        hasrxn, text = hasreaction(url);
        rxns, rxns_url = getreactions(text)

        for k ∈ 1:length(rxns)
            rxnpath = joinpath(outpath_quantumyields, categories[i], rxns[k])
            if !isdir(rxnpath)
                mkpath(rxnpath)
            end

            rxnurl = joinpath(base_url, rxns_url[k])
            push!(urllist, rxnurl)

            try
                # for each formula, grab all datasets and generate csv with links.
                hasdata, text = hasdataset(rxnurl)
                if hasdata
                    datasets = getdatasets(text)
                    datasets = replace.(datasets, " "=>"%20", "{"=>"%7B", "}"=>"%7D")

                    summary_dict = Dict(
                        :T1 => Float64[],
                        :T2 => Float64[],
                        :λ1 => Float64[],
                        :λ2 => Float64[],
                        :author_year => String[],
                        :download_url => String[],
                        :comments => String[],
                        :reaction=> String[],
                        :name => String[],
                        :doi=> String[],
                        :fname => String[],
                    )


                    for dataset ∈ datasets
                        println("\t working on $(dataset)")
                        try
                            info = getrxninfo(dataset)

                            push!(summary_dict[:T1], info[:T1])
                            push!(summary_dict[:T2], info[:T2])
                            push!(summary_dict[:λ1], info[:λ1])
                            push!(summary_dict[:λ2], info[:λ2])
                            push!(summary_dict[:author_year], info[:author_year])
                            push!(summary_dict[:download_url], info[:download_url])
                            push!(summary_dict[:comments], info[:comments])
                            push!(summary_dict[:reaction], info[:reaction])
                            push!(summary_dict[:name], info[:name])
                            push!(summary_dict[:doi], info[:doi])
                            push!(summary_dict[:fname], info[:fname])

                            # getdata(info[:download_url], info[:fname], "../")

                        catch e
                            println("FAILURE: $(dataset)")
                            println(e)
                            push!(failure_list, dataset)
                        end
                    end

                    CSV.write(joinpath(rxnpath, "summary.csv"), summary_dict)
                end

            catch
                println("FAILURE: Couldn't read $(rxnurl)")
            end
        end
    end
end

# save file with urls
writedlm(joinpath(outpath_base, "quantum_yield_urls.csv"), urllist, ",")
writedlm(joinpath(outpath_base, "failed_urls_qy.csv"), failure_list, ",")
