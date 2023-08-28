struct PhotolysisReaction{T0<:Integer, T1<:AbstractString, T2<:AbstractString, T3<:Real, T4<:AbstractString, T5<:AbstractString, T6<:AbstractString}<: Reaction
    idx::T0
    source::String
    reactants::AbstractVector{T1}
    products::AbstractVector{T2}
    prod_stoich::AbstractVector{T3}
    autochem_files::AbstractVector{T4}
    crosssection_files::AbstractVector{T5}
    quantumyield_files::AbstractVector{T6}
end



PhotolysisReaction(rdict::Dict) = PhotolysisReaction(
    rdict["idx"],
    rdict["source"],
    String.(rdict["reactants"]),
    String.(rdict["products"]),
    convert(Vector{typeof(rdict["prod_stoich"][1])}, rdict["prod_stoich"]),
    String.(rdict["autochem_files"]),
    String.(rdict["crosssection_files"]),
    String.(rdict["quantumyield_files"])
)

# define how to convert Bimol reaction into JSON for parsing
JSON.lower(r::PhotolysisReaction) = (;
                                      idx=r.idx,
                                      source=r.source,
                                      reactants=r.reactants,
                                      products=r.products,
                                      prod_stoich=r.prod_stoich,
                                      autochem_files=r.autochem_files,
                                      crosssection_files=r.crosssection_files,
                                      quantumyield_files=r.quantumyield_files,
                                      )



"""
    parse_photolysis_d(path)

Parse an AutoChem photolysis reaction database, returning a vector of `PhotolysisReaction` objects and a second vector containing the indices of any reactions that failed to parse.
"""
function parse_photolysis_d(path)
    photo = readdlm(path)

    idx_start = findall([typeof(i)<:Integer for i ∈ photo[:,1]])
    N_rxns = length(idx_start)

    # generate reaction slices
    idx_slices = []
    for i ∈ 1:N_rxns-1
        push!(idx_slices, idx_start[i]:(idx_start[i+1]-1))
    end
    push!(idx_slices, idx_start[end]:size(photo,1))


    rxns = PhotolysisReaction[]
    rxns_failed = []

    for i ∈ 1:N_rxns
        try
            rxn_data = photo[idx_slices[i], :]

            push!(
                rxns,
                PhotolysisReaction(
                    rxn_data[1,1],
                    String(rxn_data[1,2]),
                    [String(reactant) for reactant ∈ rxn_data[2,2:end] if reactant != ""],
                    [String(product) for product ∈ rxn_data[3,2:end] if product != ""],
                    [stoich for stoich ∈ rxn_data[4,2:end] if stoich != ""],
                    String.(rxn_data[5:end,1]),
                    [""],
                    [""]
                )
            )
        catch e
            push!(rxns_failed, i)
        end

    end

    return rxns, rxns_failed
end



"""
    read_photolysis(path)

Parse a photolysis reaction database (in JSON format) returning a vector of `PhotolysisReaction` objects.
"""
function read_photolysis(path)
    rxn_dicts = JSON.parsefile(path)
    rxns = Vector{PhotolysisReaction}(undef, length(rxn_dicts))

    for i ∈ 1:length(rxn_dicts)
        rxns[i] = PhotolysisReaction(rxn_dicts[i])
    end

    return rxns
end

