using SparseArrays

function remove_duplicates(species, stoich)
    unique_species = unique(species)

    out_species = []
    out_stoich = []
    for uspec ∈ unique_species
        idxs = findall(x -> x == uspec, species)
        push!(out_species, uspec)
        push!(out_stoich, sum(stoich[idxs]))
    end
    return out_species, out_stoich
end

# function get_species_idxs(species, unique_species)
#     idxs = []
#     for spec ∈ species
#         idx = findfirst(x -> x == spec, unique_species)
#         push!(idxs, idx)
#     end

#     if idxs == [nothing]
#         idxs = nothing
#     end

#     return idxs
# end


function generate_stoich_mat(df_species, db_bimol, db_trimol, db_photo;model_name::String="autochem-w-ions")
    # is_integrated = df_species.is_integrated .== 1

    # df_integrated = df_species[is_integrated, :];

    outpath = joinpath("models", model_name, "mechanism", "N.csv")
    outpath2 = joinpath("models", model_name, "mechanism", "R.csv")
    if isfile(outpath)
        rm(outpath)
        rm(outpath2)
    end

    n_rxns = length(db_bimol) + length(db_trimol) + length(db_photo)
    N = zeros(Int, nrow(df_species), n_rxns)

    # loop over reactions
    i = 1
    for rxn ∈ db_bimol
        reactants_stoich = [1 for _ ∈ rxn.reactants]
        products_stoich = rxn.prod_stoich

        # check for duplicates and update
        reactants, reactants_stoich = remove_duplicates(rxn.reactants, reactants_stoich)
        products, products_stoich = remove_duplicates(rxn.products, products_stoich)


        # generate index lookups for reactants and products
        # idxs_reactants = get_species_idxs(reactants, df_integrated.varname)
        # idxs_products = get_species_idxs(products, df_integrated.varname)
        idxs_reactants = reactants
        idxs_products = products

        # update N
        if idxs_reactants != nothing
            for j ∈ 1:size(idxs_reactants,1)
                if idxs_reactants[j] != nothing
                    N[idxs_reactants[j],i] = -1*reactants_stoich[j]
                end
            end
        end

        if idxs_products != nothing
            for j ∈ 1:size(idxs_products,1)
                if idxs_products[j] != nothing
                    N[idxs_products[j], i] = products_stoich[j]
                end
            end
        end

        i+=1
    end

    for rxn ∈ db_trimol
        reactants_stoich = [1 for _ ∈ rxn.reactants]
        products_stoich = rxn.prod_stoich

        # check for duplicates and update
        reactants, reactants_stoich = remove_duplicates(rxn.reactants, reactants_stoich)
        products, products_stoich = remove_duplicates(rxn.products, products_stoich)

        # generate index lookups for reactants and products
        idxs_reactants = reactants
        idxs_products = products

        # update N
        if idxs_reactants != nothing
            for j ∈ 1:size(idxs_reactants,1)
                if idxs_reactants[j] != nothing
                    N[idxs_reactants[j],i] = -1*reactants_stoich[j]
                end
            end
        end

        if idxs_products != nothing
            for j ∈ 1:size(idxs_products,1)
                if idxs_products[j] != nothing
                    N[idxs_products[j], i] = products_stoich[j]
                end
            end
        end

        i+=1

    end

    for rxn ∈ db_photo
        reactants_stoich = [1 for _ ∈ rxn.reactants]
        products_stoich = rxn.prod_stoich

        # check for duplicates and update
        reactants, reactants_stoich = remove_duplicates(rxn.reactants, reactants_stoich)
        products, products_stoich = remove_duplicates(rxn.products, products_stoich)

        # generate index lookups for reactants and products
        idxs_reactants = reactants
        idxs_products = products

        # update N
        if idxs_reactants != nothing
            for j ∈ 1:size(idxs_reactants,1)
                if idxs_reactants[j] != nothing
                    N[idxs_reactants[j],i] = -1*reactants_stoich[j]
                end
            end
        end

        if idxs_products != nothing
            for j ∈ 1:size(idxs_products,1)
                if idxs_products[j] != nothing
                    N[idxs_products[j], i] = products_stoich[j]
                end
            end
        end

        i+=1

    end

    # sparsify
    N = sparse(N)
    I,J,V = findnz(N)

    # write sparse matrix to csv file for easy loading
    df = DataFrame([:I => I, :J => J, :V => V])
    CSV.write(outpath, df)

    # now we want to generate our reactants stoich matrix
    R = copy(N)
    R[R .> 0 ] .= 0
    R = -1 .* R
    I,J,V = findnz(R)

    df = DataFrame([:I => I, :J => J, :V => V])
    CSV.write(outpath2, df)

    return N
end


