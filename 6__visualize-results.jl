using AutoChem
using DelimitedFiles, CSV, DataFrames
using JSON
using LinearAlgebra
using ParameterHandling, EarlyStopping
using StableRNGs
using BenchmarkTools
using ProgressMeter
using CairoMakie
using MintsMakieRecipes

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





# observation: ~1000 negative ions per cc is what we expect for outside... This definitely doesn't apply to our chamber.
df_species[idx_meas, :]


idx_to_use = 6
f1, ax, b1 = band(ts[idx_0:end] ./ 60.0, ua_mr_vals[idx_meas[idx_to_use], :] .- ua_mr_ϵ[idx_meas[idx_to_use], :], ua_mr_vals[idx_meas[idx_to_use],:] .+ ua_mr_ϵ[idx_meas[idx_to_use],:]; color=(mints_colors[2],0.2))
l1 = lines!(ax, ts[idx_0:end] ./ 60.0, ua_mr_vals[idx_meas[idx_to_use],:], lw=3, color=mints_colors[2])

eb1 = errorbars!(ax, ts[idx_0:end] ./ 60.0, W_mr_val[idx_to_use,:], W_mr_ϵ[idx_to_use,:], color=mints_colors[2])
sc1 = scatter!(ax, ts[idx_0:end] ./ 60.0, W_mr_val[idx_to_use,:], color=mints_colors[2])
f1


fig, ax, b = band(ts[idx_0:end] ./60.0, pos_tot_val[:] .- pos_tot_ϵ[:], pos_tot_val[:] .+ pos_tot_ϵ[:]; color=(mints_colors[1],0.2), axis=(xlabel="time (hours)", ylabel="Ions (cm⁻³)"))
l = lines!(ax, ts[idx_0:end] ./60.0, pos_tot_val[:], lw=3)
eb = errorbars!(ax, ts[idx_0:end] ./60.0, df_ions[:,"Total Positive Ions (measured)"], df_ions_ϵ[:,"Total Positive Ions (measured)"])
sc = scatter!(ax, ts[idx_0:end] ./60.0, df_ions[:,"Total Positive Ions (measured)"])
fig


b2 = band!(ax, ts[idx_0:end] ./ 60.0, neg_tot_val[:] .- neg_tot_ϵ[:], neg_tot_val[:] .+ neg_tot_ϵ[:]; color=(mints_colors[2],0.2))
l2 = lines!(ax, ts[idx_0:end] ./ 60.0, neg_tot_val[:], lw=3, color=mints_colors[2])
eb2 = errorbars!(ax, ts[idx_0:end] ./ 60.0, df_ions[:,"Total Negative Ions (measured)"], df_ions_ϵ[:,"Total Negative Ions (measured)"], color=mints_colors[2])
sc2 = scatter!(ax, ts[idx_0:end] ./ 60.0, df_ions[:,"Total Negative Ions (measured)"], color=mints_colors[2])

axislegend(ax, [l, l2, sc, sc2], ["Positive (modeled)", "Negative (modeled)", "Positive (observed)", "Negative (observed)"])

fig

save(joinpath(outpath, "EKF", "ion_totals.png"), fig)
save(joinpath(outpath, "EKF", "ion_totals.svg"), fig)





fig = Figure();
ax = Axis(fig[1,1], xlabel="time with ActivePure (hours)", ylabel="Concentration (ppqv)", title="NO₃")
ax2 = Axis(fig[1,2], xlabel="time with ActivePure (hours)", ylabel= "Concentration (ppqv)", title="NO₃⁻")


b1 = band!(ax, ts[idx_0:end] ./ 60.0, (df_ekf[:, "NO_3"]  .- df_ekf_ϵ[:, "NO_3"])./1e-15,  (df_ekf[:, "NO_3"]  .+ df_ekf_ϵ[:, "NO_3"])./1e-15; color=(mints_colors[1],0.2))
l1 = lines!(ax, ts[idx_0:end] ./ 60.0, df_ekf[:, "NO_3"] ./ 1e-15, linewidth=3)

b2 = band!(ax2, ts[idx_0:end] ./ 60.0, (df_ekf[:, "NO_3^-"]  .- df_ekf_ϵ[:, "NO_3^-"])./1e-15,  (df_ekf[:, "NO_3^-"]  .+ df_ekf_ϵ[:, "NO_3^-"])./1e-15; color=(mints_colors[2],0.2))
l2 = lines!(ax2, ts[idx_0:end] ./ 60.0, df_ekf[:, "NO_3^-"] ./ 1e-15, color=mints_colors[2], linewidth=3)

ylims!(ax, 0, 150)
ylims!(ax2, 0, 1)

xlims!(ax, ts[idx_0]./60, ts[end]./60)
xlims!(ax2, ts[idx_0]./60, ts[end]./60)

fig

save("NO3_comparison.png", fig)

f, ax, l = lines(ts./60, df_params.temperature, axis=(; ylabel="Temperature (°C)", xlabel="Time (hours)"), linewidth=3, color=mints_colors[3])
save("Temperature.png", f)



# 10_000 / U_noint[end,1]
# [296, 552] ./ U_noint[end,1]

fig = Figure();
ax = Axis(fig[1,1], xlabel="t (hours)", ylabel="NO₃⁻ Lifetime (s)");
ax2 = Axis(fig[1,1], ylabel="Temperature (°C)", yaxisposition=:right);
linkxaxes!(ax, ax2)
lines!(ax, ts[idx_0:end] ./60.0, df_τs[:,"NO_3^-"], linewidth=3)
lines!(ax2,ts[idx_0:end]./60.0, df_params.temperature[idx_0:end], color=mints_colors[2], linewidth=3)

fig



fig = Figure();
ax = Axis(fig[1,1], xlabel="t (hours)", ylabel="NO₃ Lifetime (s)");
ax2 = Axis(fig[1,1], ylabel="NO₃⁻ Lifetime (s)", yaxisposition=:right);
linkxaxes!(ax, ax2)
lines!(ax, ts[idx_0:end] ./60.0, df_τs[:,"NO_3"], linewidth=3)
lines!(ax2,ts[idx_0:end]./60.0, df_τs[:,"NO_3^-"], color=mints_colors[2], linewidth=3)

fig


save("NO3-comparison.pdf", fig)



fig = Figure();
ax = Axis(fig[1,1], xlabel="t (hours)", ylabel="NO₃ Lifetime");
ax2 = Axis(fig[1,1], ylabel="NO₃⁻ Lifetime", yaxisposition=:right);
linkxaxes!(ax, ax2)
lines!(ax, ts[idx_0:end] ./60.0, df_τs[:,"NO_3"], linewidth=3)
lines!(ax2,ts[idx_0:end]./60.0, df_τs[:,"NO_3^-"], color=mints_colors[2], linewidth=3)

fig



