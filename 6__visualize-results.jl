using AutoChem
using CSV, DataFrames
using CairoMakie
using MintsMakieRecipes
using LaTeXStrings
using ProgressMeter

set_theme!(mints_theme)
update_theme!(
    figure_padding=30,
    Axis = (
        xticklabelsize=20,
        yticklabelsize=20,
        xlabelsize=22,
        ylabelsize=22,
        titlesize=25,
    ),
    Colorbar = (
        ticklabelsize=20,
        labelsize=22
    )
)




# set up model output directory
@info "Setting up file paths..."
collection_id = "empty"
unc_ext = "_std"
#model_name = "autochem-w-ions"
model_name = "methane"
model_path= "models"
outpath = joinpath(model_path, model_name, "runs", collection_id)
docs_path = joinpath(model_path, model_name, "docs")
figs_path = joinpath(outpath, "figures")
mkpath(joinpath(figs_path, "concentrations"))
mkpath(joinpath(figs_path, "lifetimes"))
mkpath(joinpath(figs_path, "production-rates"))
mkpath(joinpath(figs_path, "loss-rates"))

ekf_path = joinpath(outpath, "EKF")

@assert ispath(docs_path)
@assert ispath(figs_path)
@assert ispath(ekf_path)



# read in dataframes
@info "Loading data into DataFrames"

df_species = CSV.read(joinpath(outpath, "mechanism", "species.csv"), DataFrame);
df_params= CSV.read(joinpath(outpath, "mechanism", "state_parameters.csv"), DataFrame);
df_params_ϵ = CSV.read(joinpath(outpath, "mechanism", "state_parameters_ϵ.csv"), DataFrame);

df_ua = CSV.read(joinpath(ekf_path, "ekf_output.csv"), DataFrame)
df_ua_ϵ = CSV.read(joinpath(ekf_path, "ekf_ϵ_output.csv"), DataFrame)
df_w = CSV.read(joinpath(ekf_path, "ekf_measurements.csv"), DataFrame)
df_w_ϵ = CSV.read(joinpath(ekf_path, "ekf_measurements_ϵ.csv"), DataFrame)
df_τ = CSV.read(joinpath(ekf_path, "lifetimes.csv"), DataFrame)

# generate measurement indices
idx_meas = []
for meas ∈ names(df_w)
    if meas != "times"
        println(meas)
        idx = findfirst(df_species.printname.== meas)
        push!(idx_meas, idx)
    end
end
unique!(idx_meas)



ts = df_ua.times
ua_mr_vals = Matrix(df_ua[:, Not([:times])])'
ua_mr_ϵ = Matrix(df_ua_ϵ[:, Not([:times])])'
W_mr_val = Matrix(df_w[:, Not([:times])])'
W_mr_ϵ = Matrix(df_w_ϵ[:, Not([:times])])'




# generate parameter plots:
names(df_params)
idxs = [t > ts[1] for t ∈ df_params.t]

f, ax, l = lines(df_params.t[idxs]./60, df_params.temperature[idxs], axis=(; ylabel="Temperature (K)", xlabel="Time (hours)"), linewidth=3, color=mints_colors[3])
xlims!(ax, ts[1]/60, ts[end]/60)

save(joinpath(figs_path, "temperature.png"), f)
save(joinpath(figs_path, "temperature.svg"), f)
save(joinpath(figs_path, "temperature.pdf"), f)

f, ax, l = lines(df_params.t[idxs]./60, df_params.pressure[idxs], axis=(; ylabel="Pressure (mbar)", xlabel="Time (hours)"), linewidth=3, color=mints_colors[3])
xlims!(ax, ts[1]/60, ts[end]/60)

save(joinpath(figs_path, "pressure.png"), f)
save(joinpath(figs_path, "pressure.svg"), f)
save(joinpath(figs_path, "pressure.pdf"), f)






# Generate time series plots
ts_plot = ts ./ 60
ts_label = "time (hours)"

@showprogress for i ∈ 1:nrow(df_species)
    if df_species.is_integrated[i] == 1

        units, unit_mult = get_reasonable_mr_units(ua_mr_vals[i,:])
        print_name  = df_species[i, :varname]
        spec_name = df_species[i, :printname]
        latex_label= latexstring("\\mathrm{"*spec_name*"}\\text{ "*"($(units))}")


        fig = Figure();
        gl = fig[1,1] = GridLayout();
        ax_main = Axis(gl[1,1], xlabel=ts_label, ylabel=latex_label, title="Concentration Time Series");
        ax_right = Axis(gl[1,2],
                        leftspinevisible = false,
                        rightspinevisible = false,
                        bottomspinevisible = false,
                        topspinevisible = false,
                        );

        linkyaxes!(ax_main, ax_right);
        hidedecorations!(ax_right);
        colgap!(gl, 5);
        colsize!(gl, 2, Relative(0.15));
        b = band!(ax_main, ts_plot, (ua_mr_vals[i,:] .- ua_mr_ϵ[i,:]).* unit_mult, (ua_mr_vals[i,:] .+ ua_mr_ϵ[i,:]).* unit_mult, color=(mints_colors[1], 0.25))
        l = lines!(ax_main, ts_plot, ua_mr_vals[i,:] .* unit_mult, linewidth=3, color=mints_colors[1])

        d = density!(ax_right, ua_mr_vals[i,:] .* unit_mult, direction=:y, color=(mints_colors[1], 0.25), strokecolor=mints_colors[1], strokewidth=2)

        if i ∈ idx_meas
            idx_to_use = findfirst(x->x==i, idx_meas)
            s = scatter!(ax_main, ts_plot[2:end], W_mr_val[idx_to_use,2:end] .* unit_mult, color=mints_colors[2], markerwidth=2)
            err = errorbars!(ax_main, ts_plot[2:end], W_mr_val[idx_to_use,2:end] .* unit_mult, W_mr_ϵ[idx_to_use,2:end] .* unit_mult, color=mints_colors[2], whiskerwidth=5)

            # if we have measurements, add a legend
            #leg = Legend(fig[1,2], [l, s], ["EKF", "Measurements"])
            leg = axislegend(ax_main, [l, s], ["EKF", "Measurements"]; position = :lt, framevisible=false)
        end

        fig
        xlims!(ax_right, 0, nothing)
        xlims!(ax_main, ts_plot[1], ts_plot[end])
        fig

        save(joinpath(figs_path, "concentrations", "$(print_name).png"), fig)
        save(joinpath(figs_path, "concentrations", "$(print_name).pdf"), fig)
        save(joinpath(figs_path, "concentrations", "$(print_name).svg"), fig)
    end
end



# plot negative ion concentrations
df_species[idx_meas, :]

df_ions = CSV.read(joinpath(ekf_path, "ions.csv"), DataFrame)
df_ions_ϵ = CSV.read(joinpath(ekf_path, "ions_ϵ.csv"), DataFrame)


fig = Figure();
ax = Axis(fig[1,1], xlabel="time (hours)", ylabel="Negative Ions (cm⁻³)", title="Concentration Time Series");
b = band!(ax, ts ./ 60, df_ions[:, 3] .- df_ions_ϵ[:,3], df_ions[:,3] .+ df_ions_ϵ[:,3], color=(mints_colors[1], 0.25))
l = lines!(ax, ts ./ 60, df_ions[:,3], color=mints_colors[1], linewidth=3)
eb = errorbars!(ax, ts ./60.0, df_ions[:,4], df_ions_ϵ[:,4], color=mints_colors[2], whiskerwidth=5)
s = scatter!(ax, ts ./60.0, df_ions[:,4], color=mints_colors[2], markerwidth=2)
leg = axislegend(ax, [l ,s], ["EKF", "Measured"], framevisible=false)

save(joinpath(figs_path, "negative-ions.png"), fig)
save(joinpath(figs_path, "negative-ions.svg"), fig)
save(joinpath(figs_path, "negative-ions.pdf"), fig)

fig = Figure();
ax = Axis(fig[1,1], xlabel="time (hours)", ylabel="Positive Ions (cm⁻³)", title="Concentration Time Series");
b = band!(ax, ts ./ 60, df_ions[:, 1] .- df_ions_ϵ[:,1], df_ions[:,1] .+ df_ions_ϵ[:,1], color=(mints_colors[1], 0.25))
l = lines!(ax, ts ./ 60, df_ions[:,1], color=mints_colors[1], linewidth=3)
eb = errorbars!(ax, ts ./60.0, df_ions[:,2], df_ions_ϵ[:,2], color=mints_colors[2], whiskerwidth=5)
s = scatter!(ax, ts ./60.0, df_ions[:,2], color=mints_colors[2], markerwidth=2)
leg = axislegend(ax, [l ,s], ["EKF", "Measured"], framevisible=false)

save(joinpath(figs_path, "positive-ions.png"), fig)
save(joinpath(figs_path, "positive-ions.svg"), fig)
save(joinpath(figs_path, "positive-ions.pdf"), fig)





# Generate  lifetime plots
τs = Matrix(df_τ[!, Not([:times])])'

ts_plot = ts ./ 60
ts_label = "time (hours)"

@showprogress for i ∈ 1:nrow(df_species)
    if df_species.is_integrated[i] == 1
        units, unit_mult = get_reasonable_time_units(τs[i,:])
        print_name  = df_species[i, :varname]
        spec_name = df_species[i, :printname]
        latex_label= latexstring("\\mathrm{"*spec_name*"}\\text{ Lifetime "*"($(units))}")

        try
            fig = Figure();
            gl = fig[1,1] = GridLayout();
            ax_main = Axis(gl[1,1], xlabel=ts_label, ylabel=latex_label, title="Lifetime Time Series");
            ax_right = Axis(gl[1,2],
                            leftspinevisible = false,
                            rightspinevisible = false,
                            bottomspinevisible = false,
                            topspinevisible = false,
                            );

            linkyaxes!(ax_main, ax_right);
            hidedecorations!(ax_right);
            colgap!(gl, 5);
            colsize!(gl, 2, Relative(0.15));

            l = lines!(ax_main, ts_plot, τs[i,:] .* unit_mult, linewidth=3, color=mints_colors[1])
            d = density!(ax_right, τs[i,:] .* unit_mult, direction=:y, color=(mints_colors[1], 0.25), strokecolor=mints_colors[1], strokewidth=2)


            fig
            xlims!(ax_right, 0, nothing)
            xlims!(ax_main, ts_plot[1], ts_plot[end])
            fig

            save(joinpath(figs_path, "lifetimes", "$(print_name).png"), fig)
            save(joinpath(figs_path, "lifetimes", "$(print_name).pdf"), fig)
            save(joinpath(figs_path, "lifetimes", "$(print_name).svg"), fig)

        catch e
            println("$(print_name) failed to plot lifetime!")
            println(e)
        end
    end
end





# Generate Relative Contribution plots as a feature importance ranking
using DelimitedFiles
using Statistics
using LaTeXStrings


const db_bimol = read_bimol(joinpath(outpath, "mechanism", "bimol.json"));
const db_trimol = read_trimol(joinpath(outpath, "mechanism", "trimol.json"));
const db_photo = read_fitted_photolysis(joinpath(outpath, "mechanism", "fitted_photolysis.json"));


const derivatives_bimol = get_bimolecular_derivatives(db_bimol, df_species)
const derivatives_trimol = get_trimolecular_derivatives(db_trimol, df_species)
const derivatives_photo = get_photolysis_derivatives(db_photo, df_species)


du_bi = readdlm(joinpath(outpath, "EKF", "du_bi.csv"), ',')
du_tri = readdlm(joinpath(outpath, "EKF", "du_tri.csv"), ',')
du_photo = readdlm(joinpath(outpath, "EKF", "du_photo.csv"), ',')


# take averages
du_bi_μ = zeros(size(du_bi, 1))

minimum(du_bi)

du_tri_μ = zeros(size(du_tri, 1))
du_photo_μ = zeros(size(du_photo, 1))

for i ∈ axes(du_bi,1)
    du_bi_μ[i] = mean([du for du ∈ du_bi[i,:] if !isnan(du) && !isinf(du)])
end

for i ∈ axes(du_tri,1)
    du_tri_μ[i] = mean([du for du ∈ du_tri[i,:] if !isnan(du) && !isinf(du)])
end

for i ∈ axes(du_photo,1)
    du_photo_μ[i] = mean([du for du ∈ du_photo[i,:] if !isnan(du) && !isinf(du)])
end


# idx_min = argmin(du_bi[:,1])


# fig,ax,h = hist(du_bi_μ)
# ylims!(ax, -1, 1)
# fig


function get_reactants_tex(rxn, df_species)
    reactants = df_species.printname[rxn.reactants]
    products = df_species.printname[rxn.products]

    reactants = [ ("h\\nu" == r) ? "h\\nu" : "\\mathrm{$r}" for r ∈ reactants]
    products = ["\\mathrm{$p}" for p ∈ products]
    pstoich = Int.(rxn.prod_stoich)

    n_involved = length(reactants) + length(products)

    for i ∈ 1:length(products)
        if pstoich[i] > 1
            products[i] = "$(pstoich[i])" * products[i]
        end
    end

    reactants = join([r for r ∈ reactants], " + ")
    products = join([p for p ∈ products], " + ")

    return LaTeXString("\$"* reactants * "\\longrightarrow" * products *"\$")
end

size(derivatives_bimol)
size(db_bimol)

@showprogress for j ∈ 1:nrow(df_species)
    if df_species.is_integrated[j] == 1
        idx_species = df_species.idx_species[j]
        print_name = df_species.printname[j]
        save_name = df_species.varname[j]

        rxn_strings = LaTeXString[]
        rxn_freq = Float64[]

        for i ∈ axes(derivatives_bimol,1)
            deriv = derivatives_bimol[i]
            if deriv.idx_du == idx_species
                push!(rxn_strings, get_reactants_tex(db_bimol[deriv.idx_k], df_species))
                push!(rxn_freq, du_bi_μ[i])
            end
        end

        for i ∈ axes(derivatives_trimol,1)
            deriv = derivatives_trimol[i]
            if deriv.idx_du == idx_species
                push!(rxn_strings, get_reactants_tex(db_trimol[deriv.idx_k], df_species))
                push!(rxn_freq, du_tri_μ[i])
            end
        end

        for i ∈ axes(derivatives_photo,1)
            deriv = derivatives_photo[i]
            if deriv.idx_du == idx_species
                push!(rxn_strings, get_reactants_tex(db_photo[deriv.idx_k], df_species))
                push!(rxn_freq, du_photo_μ[i])
            end
        end



        try

            idx_prod = findall(rxn_freq .> 0.0)
            idx_loss = findall(rxn_freq .< 0.0)

            rxns_prod = rxn_strings[idx_prod]
            prod_rates = rxn_freq[idx_prod]

            rxns_loss = rxn_strings[idx_loss]
            loss_rates = rxn_freq[idx_loss]


            fig = Figure();
            if length(rxns_prod) > 10
                fig = Figure(; resolution=(1000,1000));
            end
            ax = Axis(fig[1,1],
                      yticks=(1:length(rxns_prod), rxns_prod),
                      xlabel="Production Rate (cm⁻³⋅s⁻¹)",
                      yminorgridvisible=false,
                      );
            b = barplot!(
                ax,
                1:length(prod_rates),
                prod_rates,
                direction=:x
            )

            save(joinpath(figs_path, "production-rates", "$(save_name).png"), fig)
            save(joinpath(figs_path, "production-rates", "$(save_name).pdf"), fig)
            save(joinpath(figs_path, "production-rates", "$(save_name).svg"), fig)

            fig = Figure();
            if length(rxns_loss) > 10
                fig = Figure(; resolution=(1000,1000));
            end
            ax = Axis(fig[1,1],
                      yticks=(1:length(rxns_loss), rxns_loss),
                      xlabel="Loss Rate (cm⁻³⋅s⁻¹)",
                      yminorgridvisible=false,
                      );
            b = barplot!(
                ax,
                1:length(loss_rates),
                loss_rates,
                direction=:x
            )

            save(joinpath(figs_path, "loss-rates", "$(save_name).png"), fig)
            save(joinpath(figs_path, "loss-rates", "$(save_name).pdf"), fig)
            save(joinpath(figs_path, "loss-rates", "$(save_name).svg"), fig)

        catch e
            println("$(print_name) failed to plot lifetime!")
            println(e)
        end
    end
end





@info "generating analysis doc"
# set up javascript output
qmd_path = joinpath(outpath, "output.qmd")

open(qmd_path, "w") do f
    header = """---
title: "Assimilation Results"
subtitle: "Mechanism: $(model_name), Collection: $(collection_id)"
author: "John Waczak"
date: today
format:
    html:
        page-layout: full
execute:
    echo: false
    output: true
---

"""
    println(f, header)

    println(f, "```{ojs}")
    println(f, "df_species =  FileAttachment(\"./mechanism/species.csv\").csv({typed: true});")
    println(f, "df_ekf =  FileAttachment(\"./EKF/ekf_output.csv\").csv({typed: true});")
    println(f, "df_ekf_unc =  FileAttachment(\"./EKF/ekf_ϵ_output.csv\").csv({typed: true});")
    println(f, "df_meas =  FileAttachment(\"./EKF/ekf_measurements.csv\").csv({typed: true});")
    println(f, "df_meas_unc =  FileAttachment(\"./EKF/ekf_measurements_ϵ.csv\").csv({typed: true});")

    println(f, "concentration_plots = ({")
    for row ∈ eachrow(df_species)
        if row.is_integrated == 1
            println(f, "\t\"$(row.printname)\": FileAttachment(\"./figures/concentrations/$(row.varname).svg\"),")
        end
    end
    println(f, "});")
    println(f, "\n")


    println(f, "lifetime_plots = ({")
    for row ∈ eachrow(df_species)
        if row.is_integrated == 1
            println(f, "\t\"$(row.printname)\": FileAttachment(\"./figures/lifetimes/$(row.varname).svg\"),")
        end
    end
    println(f, "});")
    println(f, "\n")


    println(f, "production_plots = ({")
    for row ∈ eachrow(df_species)
        if row.is_integrated == 1
            println(f, "\t\"$(row.printname)\": FileAttachment(\"./figures/production-rates/$(row.varname).svg\"),")
        end
    end
    println(f, "});")
    println(f, "\n")


    println(f, "loss_plots = ({")
    for row ∈ eachrow(df_species)
        if row.is_integrated == 1
            println(f, "\t\"$(row.printname)\": FileAttachment(\"./figures/loss-rates/$(row.varname).svg\"),")
        end
    end
    println(f, "});")
    println(f, "\n")




    println(f, "```")

    println(f, "```{ojs}")
    println(f, "filtered_species = df_species.filter(r => r.is_integrated == 1)");
    println(f, "viewof species = Inputs.select(filtered_species, {label: \"Species:\", format: x => x.varname, value: filtered_species[0]})");
    println(f, "```")

    println(f, "::: {layout-ncol=2}")
    println(f, "```{ojs}")
    println(f, "concentration_plots[species.printname].image()")
    println(f, "```")

    println(f, "```{ojs}")
    println(f, "lifetime_plots[species.printname].image()")
    println(f, "```\n")
    println(f, ":::\n")


    println(f, "::: {layout-ncol=2}")
    println(f, "```{ojs}")
    println(f, "production_plots[species.printname].image()")
    println(f, "```")

    println(f, "```{ojs}")
    println(f, "loss_plots[species.printname].image()")
    println(f, "```\n")
    println(f, ":::\n")


    plot_test = """
```{ojs}
data_to_plot = df_ekf.map((r) => ({time: r.times/60, concentration : r[species.printname]}));
```

    ```{ojs}
Plot.plot({
    marks: [
        Plot.lineY(data_to_plot, {x: "time", y: "concentration"})
    ]
})
```
"""
    println(f, plot_test)
end
