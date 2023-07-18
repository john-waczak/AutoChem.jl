struct BimolecularReaction{T0<:Integer, T1<:AbstractString, T2<:AbstractString, T3<:Real, T4<:Real}<: Reaction
    idx::T0
    source::String
    reactants::AbstractVector{T1}
    products::AbstractVector{T2}
    prod_stoich::AbstractVector{T3}
    a1::T4
    a2::T4
    a3::T4
    a4::T4
    a5::T4
end


BimolecularReaction(rdict::Dict) = BimolecularReaction(
    rdict["idx"],
    rdict["source"],
    String.(rdict["reactants"]),
    String.(rdict["products"]),
    convert(Vector{typeof(rdict["prod_stoich"][1])}, rdict["prod_stoich"]),
    rdict["a1"],
    rdict["a2"],
    rdict["a3"],
    rdict["a4"],
    rdict["a5"]
)

# define how to convert Bimol reaction into JSON for parsing
JSON.lower(r::BimolecularReaction) = (;
                                      idx=r.idx,
                                      source=r.source,
                                      reactants=r.reactants,
                                      products=r.products,
                                      prod_stoich=r.prod_stoich,
                                      a1=r.a1,
                                      a2=r.a2,
                                      a3=r.a3,
                                      a4=r.a4,
                                      a5=r.a5
                                      )





"""
    parse_bimol_d(path)

Parse an AutoChem bimolecular reaction database, returning a vector of `BimolecularReaction` objects and a second vector containing the indices of any reactions that failed to parse.
"""
function parse_bimol_d(path)
    bi = readdlm(path)

    # bimol reactions come in blocks of 5 rows
    N_rxns = size(bi, 1)÷5

    # generate reaction slices
    idx_slices = [i:i+4 for i ∈ 1:5:N_rxns*5]

    rxns = BimolecularReaction[]
    rxns_failed = []

    for i ∈ 1:N_rxns
        try
            rxn_data = bi[idx_slices[i], :]
            coeffs = [val for val ∈ rxn_data[5,:] if typeof(val) <: Real]
            @assert length(coeffs) == 5
            push!(
                rxns,
                BimolecularReaction(
                    rxn_data[1,1],
                    String(rxn_data[1,2]),
                    [String(reactant) for reactant ∈ rxn_data[2,2:end] if reactant != ""],
                    [String(product) for product ∈ rxn_data[3,2:end] if product != ""],
                    [stoich for stoich ∈ rxn_data[4,2:end] if stoich != ""],
                    coeffs[1],
                    coeffs[2],
                    coeffs[3],
                    coeffs[4],
                    coeffs[5],
                )
            )
        catch e
            push!(rxns_failed, i)
        end

    end

    return rxns, rxns_failed
end



"""
    read_bimol(path)

Parse a bimolecular reaction database (in JSON format) returning a vector of `BimolecularReaction` objects.
"""
function read_bimol(path)
    rxn_dicts = JSON.parsefile(path)
    rxns = Vector{BimolecularReaction}(undef, length(rxn_dicts))

    for i ∈ 1:length(rxn_dicts)
        rxns[i] = BimolecularReaction(rxn_dicts[i])
    end

    return rxns
end

