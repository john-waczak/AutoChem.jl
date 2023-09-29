
function get_tex(rxn, df_species)
    reactants = df_species.printname[rxn.reactants]
    products = df_species.printname[rxn.products]

    reactants = [ ("h\\nu" == r) ? "h\\nu" : "\\mathrm{$r}" for r ∈ reactants]
    products = ["\\mathrm{$p}" for p ∈ products]
    pstoich = Int.(rxn.prod_stoich)

    for i ∈ 1:length(products)
        if pstoich[i] > 1
            products[i] = "$(pstoich[i])" * products[i]
        end
    end

    reactants = join([r for r ∈ reactants], " + ")
    products = join([p for p ∈ products], " + ")

    out = "\$\$ " * reactants *  "\\longrightarrow " * products  *" \$\$"
    println(out)

    return out
end



function get_reaction_tex(rxn::BimolecularReaction)
    out = "\$\$k = "

    if rxn.a1 != 0.0
        out *= "(\\textrm{$(rxn.a1)})"
    end

    if rxn.a2 != 0.0
        out *= "T^{\\textrm{$(rxn.a2)}}"
    end

    if rxn.a5 != 0.0
        out *= "\\exp(\\textrm{$(-rxn.a5)}/T)"
    end

    if rxn.a3 > 0.0 && rxn.a4 != 0.0
        out *= "(T/\\textrm{$(rxn.a3)})^{\\textrm{$(rxn.a4)}}"
    end


    if rxn.all_reactants_HO2 || (rxn.contains_CO && rxn.contains_OH)
        out *= "(1.0 + (0.6 * P/1013.25))"
    elseif rxn.contains_HONO2 && rxn.contains_OH
        out  = "\$\$ k = \\textrm{2.4e-14} \\cdot \\exp(460.0/T) + \\frac{\\textrm{6.50e-34} \\cdot M \\cdot \\exp(1335.0/T)}{1 + \\frac{\\textrm{6.50e-34} \\cdot M \\cdot \\exp(1335.0/T)}{ \\textrm{2.7e-17} \\cdot exp(2199.0/T)}}"

    end


    if out == "\$\$k = "
        out = "\$\$k = 0"
    end


    # add closing $$
    out *= " \$\$"



    return out
end


function get_reaction_tex(rxn::TrimolecularReaction)
    k₀ = "k_0 &= "
    kᵢ = "k_{\\infty} &= "
    fc = "f_c &= "

    # update kₒ
    k₀ *= "M"
    if rxn.a1 != 0.0
        k₀ *= "(\\textrm{$(rxn.a1)})"
    end

    if rxn.a2 != 0.0
        k₀ *= "T^{\\textrm{$(rxn.a2)}}"
    end

    if rxn.a5 != 0.0
        k₀ *= "\\exp(\\textrm{$(-rxn.a5)}/T)"
    end

    if rxn.a3 > 0.0
        k₀ *= "(T/\\textrm{$(rxn.a3)})^{\\textrm{$(rxn.a4)}}"
    end


    # update kᵢ
    if rxn.b1 != 0.0
        kᵢ *="(\\textrm{$(rxn.b1)})"
    end

    if rxn.b2 != 0.0
        kᵢ *= "T^{\\textrm{$(rxn.b2)}}"
    end

    if rxn.b5 != 0.0
        kᵢ *= "\\exp(\\textrm{$(-rxn.b5)}/T)"
    end

    if rxn.b3 > 0.0
        kᵢ *= "(T/\\textrm{$(rxn.b3)})^{\\textrm{$(rxn.b4)}}"
    end



    # update fc
    fc *= "(\\textrm{$(rxn.c1)}) "

    if rxn.c3 > 0.0
        fc *= "+ (\\textrm{$(rxn.c2)})\\cdot\\exp(\\textrm{-T/\\textrm{$(rxn.c3)}}) "
    end

    if rxn.c4 > 0.0
        fc *= "+ (\\textrm{$(rxn.c2)})\\cdot\\exp(\\textrm{-T/\\textrm{$(rxn.c4)}}) "
    end

    if rxn.c5 > 0.0
        fc *= "+ \\exp(\textrm{$(-rxn.c5)}/T)"
    end


    # if there's been no change, set to 0
    if k₀ == "k_0 &= "
        k₀ = "k_0 &= 0"
    end

    if kᵢ == "k_{\\infty} &= "
        kᵢ = "k_{\\infty} &= 0"
    end

    if fc == "f_c &= "
        fc = "f_c &= 0"
    end

    out = """
\\begin{align}
    $(k₀) \\\\
    $(kᵢ) \\\\
    $(fc) \\\\
\\end{align}
"""

    return out
end


