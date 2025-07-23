using Plots
using LaTeXStrings
using JLD
using Statistics
include("../FitEntEntScal.jl")


function get_data_filename(L, τ, D)
    return "exm/data/Bulk_measure/Ising/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(D).jld"
end

function bulk_dynamics_model(x, p)
    return p[1] .- π/6 .* p[2] ./(x .^2)
end

function get_system_params(τ)
    if τ == log(1 + √2)
        D = 20
    else
        D = 30
    end
    return D
end

function load_data(τ, L, with_time_series=true)
    D = get_system_params(τ)

    data_path = get_data_filename(L, τ, D)
    
    if with_time_series
        return load(data_path, 
            "average_EE_tlis", "stderr_EE_tlis", "bulk_meanEElis", 
            "bulk_stderr_EElis", "time_average_free_energy", 
            "bulk_FE", "bulk_FE_stderr", "time_FEstderr", "time_FElis")
    else
        return load(data_path, 
            "average_EE_tlis", "stderr_EE_tlis", "bulk_meanEElis", 
            "bulk_stderr_EElis", "time_average_free_energy", 
            "bulk_FE", "bulk_FE_stderr")
    end
end

function fit_casimir_central_charge(L_list=collect(8:2:24))
    FE_lis = zeros(length(L_list))
    FEstderr_lis = zeros(length(L_list))
    
    for (idx, L) in enumerate(L_list)
        _, _, _, _, _, bulk_FE, bulk_FE_stderr = load_data(τ, L, false)
        FE_lis[idx] = bulk_FE
        FEstderr_lis[idx] = bulk_FE_stderr
    end

    function bulk_dynamics_model(x, p)
        return p[1] .- π/6 .* p[2] ./(x .^2)
    end

    fit = curve_fit(bulk_dynamics_model, L_list[1:end], (FE_lis ./L_list)[1:end], [0.5, 0.5])
    fit_paras = fit.param
    fit_err = stderror(fit)
    # decimal_places = ceil(Int, -log10(fit_err[2]))
    decimal_places = 2  # Set a fixed number of decimal places for the Casimir coefficient
    @show fit_paras
    x_fit = collect(8:10:1000)
    y_fit = bulk_dynamics_model(x_fit, fit_paras)
    
    fig = plot(
        vcat(0, 1 ./ x_fit .^2), 
        vcat(fit_paras[1], y_fit),
        # ylim=(0.3042, 0.3075), 
        # ylim=(0.304, 0.308),
        color=palette(:tab10)[2], 
        linewidth=2, 
        label=latexstring("c_{\\mathrm{Casimir}} = $(round(fit_paras[2], digits=decimal_places))"), 
        framestyle=:box
    )

    # tickslis = [8, 10, 12, 16, 20, 24]
    tickslis = [8, 10, 12]
    scatter!(fig, 
        1 ./ L_list .^2, FE_lis ./ L_list, 
        yerror = FEstderr_lis ./ L_list, 
        xticks = (
            vcat(0, [1 ./ tickslis .^2...]),
            vcat([L"0"], [latexstring("\\frac{1}{$(i)^2}") for i in tickslis])
        ),
        xlabel=L"1/L^2", ylabel=L"S_c/(TL)", label=false,
        marker=:circle, markersize=6, linewidth=2, 
        color=RGB(12/255, 159/255, 250/255),
        legend_background_color=nothing, legend_foreground_color=nothing, legend=:topright
    )

    return fig
end

fig_FE_scaling = fit_casimir_central_charge(collect(8:2:12))
savefig(fig_FE_scaling, "figs/Bulk_measure/Casimir_coef_scaling.pdf")

function fit_casimir_coef()
    cclis = zeros(length(γlis))
    ccstderr_lis = zeros(length(γlis))
    cCasimir_coef = zeros(length(γlis))
    cCasimir_coef_stderr = zeros(length(γlis))

    for (idx, γ) in enumerate(γlis)
        L_list = collect(8:2:20) 

        FE_lis = zeros(length(L_list))
        FEstderr_lis = zeros(length(L_list))
        for (idy, L) in enumerate(L_list)
            data_path = get_data_filename(L, τlis[idx], get_system_params(τlis[idx]))
    

            if isfile(data_path)
                bulk_meanEElis, bulk_FE, bulk_FE_stderr = load(data_path, "bulk_meanEElis", "bulk_FE", "bulk_FE_stderr")
                FE_lis[idy] = bulk_FE
                FEstderr_lis[idy] = bulk_FE_stderr


                if L == L_list[end] && γ != 1.0
                    
                    cent, fig = fitCCEntEntScal(bulk_meanEElis, mincut=3, pbc=true)
                    cclis[idx] = cent[1]  # Extract the central charge coefficient for the last L
                    ccstderr_lis[idx] = cent[2]  # Extract the error
                end

                if γ == 1.0 && L == L_list[end]
                    cent, fig = fitpart(bulk_meanEElis, mincut=2, pbc=true, part=:odd)
                    cclis[idx] = cent[1]  # Extract the central charge coefficient for the last L
                    ccstderr_lis[idx] = cent[2]
                end
            else
                @warn "Data file not found: $data_path"
            end

        end
        
        # Fit the Casimir coefficient
        fit = curve_fit(bulk_dynamics_model, L_list[2:end], (FE_lis ./L_list)[2:end], [0.5, 0.5])
        fit_paras = fit.param
        fit_err = stderror(fit)
        cCasimir_coef[idx] = fit_paras[2]  # Extract the Casimir coefficient
        cCasimir_coef_stderr[idx] = fit_err[2]  # Extract the error
    end

    return cclis, ccstderr_lis, cCasimir_coef, cCasimir_coef_stderr
end


# Plot the central charge, Casimir coefficient and their error
function plot_cc_Casimir()
    cclis, ccstderr_lis, cCasimir_coef, cCasimir_coef_stderr = fit_casimir_coef()

    # Plot the central charge coefficients
    fig = plot(
        γlis, cclis, 
        yerror=ccstderr_lis, 
        xlabel=L"γ", ylabel="Conformal data", label=latexstring("c_{\\mathrm{ent}}"),
        xticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 1],
        linewidth=2, markersize=5, marker=:circle,
        color=RGB(12/255, 159/255, 250/255),
        legend_background_color=nothing, legend_foreground_color=nothing
    )

    # Add the Casimir coefficient to the plot
    plot!(fig, 
        γlis, cCasimir_coef, 
        yerror=cCasimir_coef_stderr, 
        lw=2, markersize=5, marker=:circle,
        color=palette(:tab10)[2], label=latexstring("c_{\\mathrm{Casimir}}"),
    )

    return fig
end

fig_cc_Casimir = plot_cc_Casimir()
savefig(fig_cc_Casimir, "figs/Bulk_measure/Casimir_coef_gamma.pdf")