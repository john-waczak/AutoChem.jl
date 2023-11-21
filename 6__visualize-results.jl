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
    println(f, ":::")


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
