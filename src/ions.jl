function get_positive_indices(df_species)
    idxs_positive = Int[]
    for i ∈ 1:nrow(df_species)
        if occursin("+", df_species.varname[i])
            push!(idxs_positive, i)
        end
    end

    return idxs_positive
end

function get_negative_indices(df_species)
    idxs_negative = Int[]
    for i ∈ 1:nrow(df_species)
        if occursin("-", df_species.varname[i])
            push!(idxs_negative, i)
        end
    end

    return idxs_negative
end

