struct BimolecularReaction{T0<:Integer, T1<:Union{AbstractString, Int}, T2<:Union{AbstractString, Int}, T3<:Real, T4<:Real}<: Reaction
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
    contains_OH::Bool
    contains_HONO2::Bool
    contains_CO::Bool
    all_reactants_HO2::Bool
end


"""
    (rxn::BimolecularReaction)(T,P,M)

Given a reaction `rxn` of type `BimolecularReaction`, compute the reaction rate coefficient as a function of

- `T`: temperature in *Kelvin*
- `P`: pressure in *mbar*
- `M`: the total particle number density

"""
function (rxn::BimolecularReaction)(T,P,M)
    k = rxn.a1 * (T^rxn.a2) * exp(-rxn.a5/T)

    # handle extra temp dependence
    if rxn.a3 > 0.0 && rxn.a4 != 0.0
        k = k * (T/rxn.a3)^rxn.a4
    end

    # handle special cases
    if rxn.all_reactants_HO2 || (rxn.contains_CO && rxn.contains_OH)
        # note Patm = P/P₀ with P₀=1013.25 mbar/atm
        k = k * (1.0 + (0.6 * P/1013.25))
    elseif rxn.contains_HONO2 && rxn.contains_OH
        # we should double check w/ Dr. Lary for sources for this
        z₀ = 2.4e-14 * exp(460.0/T)
        z₂ = 2.7e-17 * exp(2199.0/T)
        z₃ = 6.50e-34 * exp(1335.0/T)
        k = z₀ + (z₃*M)/(1.0 + (z₃*M/z₂))
    end

    return k
end




BimolecularReaction(rdict::Dict) = BimolecularReaction(
    rdict["idx"],
    rdict["source"],
    convert(Vector{typeof(rdict["reactants"][1])}, rdict["reactants"]),
    convert(Vector{typeof(rdict["products"][1])}, rdict["products"]),
    convert(Vector{typeof(rdict["prod_stoich"][1])}, rdict["prod_stoich"]),
    rdict["a1"],
    rdict["a2"],
    rdict["a3"],
    rdict["a4"],
    rdict["a5"],
    rdict["contains_OH"],
    rdict["contains_HONO2"],
    rdict["contains_CO"],
    rdict["all_reactants_HO2"],
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
                                      a5=r.a5,
                                      contains_OH=r.contains_OH,
                                      contains_HONO2=r.contains_HONO2,
                                      contains_CO=r.contains_CO,
                                      all_reactants_HO2=r.all_reactants_HO2
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

            reactants = [String(reactant) for reactant ∈ rxn_data[2,2:end] if reactant != ""]
            products = [String(product) for product ∈ rxn_data[3,2:end] if product != ""]


            push!(
                rxns,
                BimolecularReaction(
                    rxn_data[1,1],
                    String(rxn_data[1,2]),
                    reactants,
                    products,
                    [stoich for stoich ∈ rxn_data[4,2:end] if stoich != ""],
                    coeffs[1],
                    coeffs[2],
                    coeffs[3],
                    coeffs[4],
                    coeffs[5],
                    any(reactants .== "OH"),
                    any(reactants .== "HONO2"),
                    any(reactants .== "CO"),
                    all(reactants .== "HO2")
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



