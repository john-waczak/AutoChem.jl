struct TrimolecularReaction{T0<:Integer, T1<:AbstractString, T2<:AbstractString, T3<:Real, T4<:Real}<: Reaction
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
    b1::T4
    b2::T4
    b3::T4
    b4::T4
    b5::T4
    c1::T4
    c2::T4
    c3::T4
    c4::T4
    c5::T4
    contains_m::Bool
    contains_ClCO::Bool
    contains_ClCO3::Bool
    contains_COCl2::Bool
    contains_APINENE::Bool
    contains_BPINENE::Bool
    contains_ClCO_and_ClCO3::Bool
    contains_ClCO_and_COCl2::Bool
    contains_only_low::Bool
end




"""
    (rxn::TrimolecularReaction)(T,P,M)

Given a reaction `rxn` of type `TrimolecularReaction`, compute the reaction rate coefficient as a function of

- `T`: temperature in *Kelvin*
- `P`: pressure in *mbar*
- `M`: the total particle number density


"""
function (rxn::TrimolecularReaction)(T,P,M)
    # pre-allocate
    k = 0.0
    k₀ = 0.0
    kᵢ = 0.0
    fc = rxn.c1

    k₀ = M*rxn.a1*(T^rxn.a2)*exp(-rxn.a5/T)
    if rxn.a3 > 0.0
        k₀ = k₀ * (T/rxn.a3)^rxn.a4
    end

    if !rxn.contains_only_low
        # if there is a high and low pressure term.
        kᵢ = rxn.b1*(T^rxn.b2)*exp(-rxn.b5/T)

        if rxn.b3 > 0.0
            kᵢ = kᵢ * (T/rxn.b3)^rxn.b4
        end

        if rxn.c3 > 0.0
            fc = fc + (1 - rxn.c2)*exp(-T/rxn.c3)
        end

        if rxn.c4 > 0.0
            fc = fc + rxn.c2 * exp(-T/rxn.c4)
        end

        if rxn.c5 > 0.0
            fc = fc + exp(-rxn.c5/T)
        end

        if k₀ < 1e-36
            k = kᵢ
        elseif kᵢ < 1e-36
            k = k₀
        else
            ratio = k₀/kᵢ
            f = 10.0^(log10(fc)/(1+log10(ratio)^2))
            k = kᵢ*(ratio/(1+ratio))*f
        end
    else
        # there is only a low-pressure term
        k = k₀
    end


    # NOTE: k should now have a non-zero value


    if rxn.contains_m
        # if the reaction already contains m explicitly, divide it out
        # from the reaction rate.
        k = k / M
    end


    # now we deal with special cases
    if rxn.contains_ClCO_and_ClCO3
        zrate = 2.83e-14
        ptorr = P * 100.0  # assume P is provided in mbar
        zratio = 0.0645 * (107.0 + ptorr)
        k = zrate * ratio
    end

    if rxn.contains_ClCO_and_COCl2
        k = 2.83e-14
    end

    if rxn.contains_BPINENE && rxn.contains_APINENE
        k = 1.0/M
    end

    return k
end



TrimolecularReaction(rdict::Dict) = TrimolecularReaction(
    rdict["idx"],
    rdict["source"],
    String.(rdict["reactants"]),
    String.(rdict["products"]),
    convert(Vector{typeof(rdict["prod_stoich"][1])}, rdict["prod_stoich"]),
    rdict["a1"],
    rdict["a2"],
    rdict["a3"],
    rdict["a4"],
    rdict["a5"],
    rdict["b1"],
    rdict["b2"],
    rdict["b3"],
    rdict["b4"],
    rdict["b5"],
    rdict["c1"],
    rdict["c2"],
    rdict["c3"],
    rdict["c4"],
    rdict["c5"],
    rdict["contains_m"],
    rdict["contains_ClCO"],
    rdict["contains_ClCO3"],
    rdict["contains_COCl2"],
    rdict["contains_APINENE"],
    rdict["contains_BPINENE"],
    rdict["contains_ClCO_and_ClCO3"],
    rdict["contains_ClCO_and_COCl2"],
    rdict["contains_only_low"]
)

# define how to convert Bimol reaction into JSON for parsing
JSON.lower(r::TrimolecularReaction) = (;
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
                                       b1=r.b1,
                                       b2=r.b2,
                                       b3=r.b3,
                                       b4=r.b4,
                                       b5=r.b5,
                                       c1=r.c1,
                                       c2=r.c2,
                                       c3=r.c3,
                                       c4=r.c4,
                                       c5=r.c5,
                                       contains_m=r.contains_m,
                                       contains_ClCO=r.contains_ClCO,
                                       contains_ClCO3=r.contains_ClCO3,
                                       contains_COCl2=r.contains_COCl2,
                                       contains_APINENE=r.contains_APINENE,
                                       contains_BPINENE=r.contains_BPINENE,
                                       contains_ClCO_and_ClCO3=r.contains_ClCO_and_ClCO3,
                                       contains_ClCO_and_COCl2=r.contains_ClCO_and_COCl2,
                                       contains_only_low=r.contains_only_low
                                      )



"""
    parse_trimol_d(path)

Parse an AutoChem trimolecular reaction database, returning a vector of `TrimolecularReaction` objects and a second vector containing the indices of any reactions that failed to parse.
"""
function parse_trimol_d(path)
    tri = readdlm(path)

    # trimol reactions come in blocks of 7 rows
    N_rxns = size(tri, 1)÷7

    # generate reaction slices
    idx_slices = [i:i+6 for i ∈ 1:7:N_rxns*7]

    rxns = TrimolecularReaction[]
    rxns_failed = []

    for i ∈ 1:N_rxns
        try
            rxn_data = tri[idx_slices[i], :]
            coeffs_a = [val for val ∈ rxn_data[5,:] if typeof(val) <: Real]
            coeffs_b = [val for val ∈ rxn_data[6,:] if typeof(val) <: Real]
            coeffs_c = [val for val ∈ rxn_data[7,:] if typeof(val) <: Real]

            @assert length(coeffs_a) == 5
            @assert length(coeffs_b) == 5
            @assert length(coeffs_c) == 5


            reactants = [String(reactant) for reactant ∈ rxn_data[2,2:end] if reactant != ""]
            products = [String(product) for product ∈ rxn_data[3,2:end] if product != ""]


            push!(
                rxns,
                TrimolecularReaction(
                    rxn_data[1,1],
                    String(rxn_data[1,2]),
                    reactants,
                    products,
                    [stoich for stoich ∈ rxn_data[4,2:end] if stoich != ""],
                    coeffs_a[1],
                    coeffs_a[2],
                    coeffs_a[3],
                    coeffs_a[4],
                    coeffs_a[5],
                    coeffs_b[1],
                    coeffs_b[2],
                    coeffs_b[3],
                    coeffs_b[4],
                    coeffs_b[5],
                    coeffs_c[1],
                    coeffs_c[2],
                    coeffs_c[3],
                    coeffs_c[4],
                    coeffs_c[5],
                    any(reactants .== "m"),
                    any(reactants .== "ClCO"),
                    any(reactants .== "ClCO3"),
                    any(reactants .== "COCl2"),
                    any(reactants .== "APINENE"),
                    any(reactants .== "BPINENE"),
                    any(reactants .== "ClCO") && any(reactants .== "ClCO3"),
                    any(reactants .== "ClCO") && any(reactants .== "COCl2"),
                    all(coeffs_b .== 0.0)
                )
            )
        catch e
            push!(rxns_failed, i)
        end

    end

    return rxns, rxns_failed
end



"""
    read_trimol(path)

Parse a trimolecular reaction database (in JSON format) returning a vector of `TrimolecularReaction` objects.
"""
function read_trimol(path)
    rxn_dicts = JSON.parsefile(path)
    rxns = Vector{TrimolecularReaction}(undef, length(rxn_dicts))

    for i ∈ 1:length(rxn_dicts)
        rxns[i] = TrimolecularReaction(rxn_dicts[i])
    end

    return rxns
end
