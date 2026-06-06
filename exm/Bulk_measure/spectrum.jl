using LinearAlgebra
using Plots
using LaTeXStrings
using Measurements
using LsqFit

Llis = collect(8:2:16)
GSE = zeros(length(Llis))
GSEerr= zeros(length(Llis))
for (i, L) in enumerate(Llis)
    data = load("exm/data/Bulk_measure/tf_spectrum_Born/L$(L)/ensemble_spectrum_L$(L)_gammaind10.jld2")
    avg_spectrum = data["avg_spectrum"]   
    stderr_spectrum = data["stderr_spectrum"]
    GSE[i] = avg_spectrum[end]
    GSEerr[i] = stderr_spectrum[end]
end

scatter(1 ./Llis .^2, GSE./Llis, yerror=GSEerr ./Llis)

function plot_scaling_dimension_Born_level(Llis, inds, λlength=9)
    α_lis = [0.00808169631522918 ± 0.00027268011041641855,
    0.0336858093655925 ± 0.0010069937182724683, 
    0.07677253264064307 ± 0.0024462134009264357,
    0.1340252511327597 ± 0.0008904660607937714,
    0.23672155605486694 ± 0.005984075300077113,
    0.35629717664541705 ± 0.00776122153413414,
    0.5648689542913156 ± 0.029405819478422905, 
    0.8943680189331076 ± 0.02930979355013548,
    1.508584182807505 ±  0.03482208339523906,
    2.096879171922712 ± 0.07620256174216884,
    3.9901816100072516 ± 0.5112930314418107,
    4.85864043594799 ± 0.44470089488871173,
    ]

    λLlis = zeros(length(Llis), λlength+1)
    λerrLlis = zeros(length(Llis), λlength+1)
    c = cgrad(:blues, λlength, categorical=true)
    fit_model(x, p) = p[1] .*x
    x_fit = range(0, stop = 2π ./minimum(Llis), length = 100)
    fig = plot(xticks = (vcat(0, 2π ./Llis), latexstring.(vcat("0", ["\\frac{2π}{$(i)}" for i in Llis]))), xlabel=L"2π/L", ylabel=L"ΔF/αT")

    for (idx, L) in enumerate(Llis)
        data = load("exm/data/Bulk_measure/tf_spectrum_Born/L$(L)/ensemble_spectrum_L$(L)_gammaind10.jld2")
        avg_spectrum = data["avg_spectrum"]   
        stderr_spectrum = data["stderr_spectrum"]
        λLlis[idx, :] = avg_spectrum
        λerrLlis[idx, :] = stderr_spectrum
    end
    @show λLlis

    for i in 1+1:λlength+1
        fit_inds = collect(2:length(Llis))
        fit = curve_fit(fit_model, 2π ./Llis[fit_inds], (λLlis[fit_inds, i] .- λLlis[fit_inds, 1])/α_lis[inds].val,  [0.5])
        y_fit = fit_model(x_fit, fit.param)
        plot!(fig, x_fit, y_fit, label = latexstring("Δ_{$(i-1)}=$(round(fit.param[1], digits=3))"), linewidth=2, color=c[i-1])
        scatter!(fig, 2π ./Llis, (λLlis[:, i] .- λLlis[:, 1])/α_lis[inds].val, yerr = λerrLlis[:, 2], label =false, markershape=:circle, color=c[i-1])
    end
    
    return fig
end
fig = plot_scaling_dimension_Born_level(Llis, 10)