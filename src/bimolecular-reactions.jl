struct BimolecularReaction{
    T0<:Integer,
    T1<:AbstractString,
    T2<:AbstractString,
    T3<:Real,
    T4<:Real}<: Reaction

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



function parse_bimol_d(path)
    bi = readdlm(path)

    # bimol reactions come in blocks of 5 rows
    N_rxns = Int(size(bi, 1)/5)

    # generate reaction slices
    idx_slices = [i:i+4 for i ∈ 1:5:N_rxns*5]

    rxns = BimolecularReaction[]
    rxns_failed = []

    for i ∈ 1:N_rxns
        try
            rxn_data = bi[idx_slices[i], :]
            push!(
                rxns,
                BimolecularReaction(
                    rxn_data[1,1],
                    String(rxn_data[1,2]),
                    [String(reactant) for reactant ∈ rxn_data[2,2:end] if reactant != ""],
                    [String(product) for product ∈ rxn_data[3,2:end] if product != ""],
                    [stoich for stoich ∈ rxn_data[4,2:end] if stoich != ""],
                    rxn_data[5,1],
                    rxn_data[5,2],
                    rxn_data[5,3],
                    rxn_data[5,4],
                    rxn_data[5,5],
                )
            )
        catch e
            push!(rxns_failed, i)
        end

    end

    return rxns, rxns_failed
end



# we'll want a JSON parser and a .d parser function
# we'll also want a JSON writer
