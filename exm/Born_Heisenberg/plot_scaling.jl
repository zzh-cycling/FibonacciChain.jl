# Central-charge and free-energy conformal scaling for the Heisenberg Born
# dynamics data. Uses the collected Observables files (ensemble-averaged EE
# profiles and steady-state bulk free energies).
#
# Fits (PBC, Calabrese-Cardy):  S(l) = (c/3) ln[(L/π) sin(π l/L)] + a
#   - even/odd bond cuts are fitted separately (strong dimerized oscillation)
# Free energy:                    F(L)  = f L + b + d/L
#   - d = -π c v / 12 (bulk_FE is the record free energy per 2 layers, cf.
#     process_data), so define cv_eff = -12 d/π.
#
# Usage: julia --project plot_scaling.jl [τ_idx]   (default ind = 10, γ=0.95)

using JLD2
using Statistics
using Printf
using LsqFit
using Plots
using LaTeXStrings

include(joinpath(@__DIR__, "config.jl"))  # only for τlis/Δlis

ind = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
τ = τlis[ind]
γ = γlis[ind]
Llis = [8, 10, 12, 14, 16]

obsdir(L) = "exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)"

function load_obs(L, Δ)
    dir = obsdir(L)
    fn = filter(f -> startswith(f, "Observables_L$(L)_τ$(τ)_Δ$(Δ)_t"), readdir(dir))
    length(fn) == 1 || error("expected 1 Observables file for L=$L Δ=$Δ, got $(length(fn))")
    return load(joinpath(dir, only(fn)))
end

# ---------- load data ----------
Sprof = Dict{Tuple{Int,Float64},Vector{Float64}}()  # (L, Δ) -> S(l), l = 1..L-1
FE = Dict{Tuple{Int,Float64},Float64}()             # (L, Δ) -> bulk_FE
for L in Llis, Δ in Δlis
    d = load_obs(L, Δ)
    Sprof[(L, Δ)] = Float64.(d["bulk_meanEElis"])
    FE[(L, Δ)] = Float64(d["bulk_FE"])
end

# ---------- linear fits ----------
# OLS slope/intercept with residual-based slope error
function linfit(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    x̄, ȳ = mean(x), mean(y)
    sxx = sum(abs2, x .- x̄)
    a = sum((x .- x̄) .* (y .- ȳ)) / sxx
    r = y .- (a .* x .+ (ȳ - a * x̄))
    return a, ȳ - a * x̄, sqrt(sum(abs2, r) / (n - 2) / sxx)
end

# shared slope across L (each L gets its own intercept)
function joint_slope(xs, ys, gs)
    num, den = 0.0, 0.0
    for g in unique(gs)
        m = gs .== g
        num += sum((xs[m] .- mean(xs[m])) .* (ys[m] .- mean(ys[m])))
        den += sum(abs2, xs[m] .- mean(xs[m]))
    end
    a = num / den
    r = Float64[]
    for g in unique(gs)
        m = gs .== g
        append!(r, ys[m] .- (a .* xs[m] .+ (mean(ys[m]) - a * mean(xs[m]))))
    end
    σ² = sum(abs2, r) / (length(xs) - length(unique(gs)) - 1)
    return a, sqrt(σ² / den)
end

xscal(l, L) = log((L / π) * sin(π * l / L)) / 3  # slope of S vs xscal = c

# cuts by branch: :all, :even, :odd
function collect_cuts(Δ, branch)
    xs, ys, gs = Float64[], Float64[], Int[]
    for L in Llis
        llist = branch == :even ? (2:2:L-2) : branch == :odd ? (1:2:L-1) : (1:L-1)
        append!(xs, xscal.(collect(llist), L))
        append!(ys, Sprof[(L, Δ)][collect(llist)])
        append!(gs, fill(L, length(llist)))
    end
    return xs, ys, gs
end

cfit = Dict{Tuple{Float64,Symbol},Tuple{Float64,Float64}}()  # (Δ, branch) -> (c, err)
for Δ in Δlis, branch in (:all, :even, :odd)
    xs, ys, gs = collect_cuts(Δ, branch)
    cfit[(Δ, branch)] = joint_slope(xs, ys, gs)
end

# ---------- free energy F(L) = f L + e + d/L ----------
# e is a non-universal constant (e.g. the -ln 4 topological constant in the
# trivial FM phase); the conformal Casimir term is d/L with
# d = -π c v / 12 (bulk_FE is the record free energy per 2 layers, cf.
# process_data), so cv_eff = -12 d/π.
fe_model(L, p) = @. p[1] * L + p[2] + p[3] / L
fefit = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()  # Δ -> (params, errs)
for Δ in Δlis
    Ls = Float64.(Llis)
    Fs = [FE[(L, Δ)] for L in Llis]
    fit = curve_fit(fe_model, Ls, Fs, [-2.0, 0.0, 0.0])
    fefit[Δ] = (fit.param, stderror(fit))
end

# ---------- fit table ----------
println("\n Δ      c_all          c_even         c_odd          const e        cv_eff = -12d/π")
for Δ in Δlis
    ca, ce, co = cfit[(Δ, :all)], cfit[(Δ, :even)], cfit[(Δ, :odd)]
    e, d = fefit[Δ][1][2], fefit[Δ][1][3]
    derr = fefit[Δ][2][3]
    @printf(
        "%5.1f  %6.3f ± %.3f  %6.3f ± %.3f  %6.3f ± %.3f  %8.4f   %7.3f ± %.3f\n",
        Δ, ca[1], ca[2], ce[1], ce[2], co[1], co[2], e, -12 * d / π, 12 * derr / π
    )
end

# ---------- figures ----------
mkpath(joinpath(@__DIR__, "figs"))
default(frame = :box, grid = false, guidefontsize = 11, tickfontsize = 9, legendfontsize = 8)
Lcolor = Dict(L => cgrad(:viridis, length(Llis), categorical = true)[i] for (i, L) in enumerate(Llis))

# fig 1: EE scaling collapse for representative Δ
showΔ = [-1.0, 0.0, 0.5, 0.8, 0.9, 1.5]
p1 = plot(layout = (2, 3), size = (1200, 700))
for (k, Δ) in enumerate(showΔ)
    for branch in (:even, :odd), L in Llis
        llist = collect(branch == :even ? (2:2:L-2) : (1:2:L-1))
        scatter!(
            p1[k], xscal.(llist, L), Sprof[(L, Δ)][llist];
            color = Lcolor[L], marker = branch == :even ? :circle : :square,
            markersize = 4, label = k == 1 && branch == :even ? "L=$L" : "",
        )
    end
    for (branch, ls) in ((:even, :solid), (:odd, :dash))
        xs, ys, gs = collect_cuts(Δ, branch)
        c = cfit[(Δ, branch)][1]
        b̄ = mean(mean(ys[gs .== L]) - c * mean(xs[gs .== L]) for L in Llis)
        plot!(
            p1[k], sort(xs), c .* sort(xs) .+ b̄;
            color = :black, linestyle = ls, lw = 1.5, label = "",
        )
    end
    ce, co = cfit[(Δ, :even)], cfit[(Δ, :odd)]
    alldata = vcat([Sprof[(L, Δ)] for L in Llis]...)
    lo, hi = extrema(alldata)
    pad = max(0.08 * (hi - lo), 0.05)
    plot!(
        p1[k]; title = latexstring("\\Delta = $Δ,\\ c_{even} = $(round(ce[1], digits=3)),\\ c_{odd} = $(round(co[1], digits=3))"),
        titlefontsize = 9, ylims = (lo - pad, hi + pad),
        xlabel = L"\frac{1}{3}\ln[(L/\pi)\sin(\pi l/L)]", ylabel = L"S(l)",
    )
end
savefig(p1, joinpath(@__DIR__, "figs", "ee_collapse_ind$(ind).png"))

# fig 2: c(Δ) and cv_eff(Δ)
p2l = plot(xlabel = L"\Delta", ylabel = L"c", size = (550, 400))
for (branch, m, lab) in ((:even, :utriangle, "c (even cuts)"), (:odd, :dtriangle, "c (odd cuts)"))
    cs = [cfit[(Δ, branch)] for Δ in Δlis]
    plot!(
        p2l, Δlis, first.(cs); yerror = last.(cs), marker = m, markersize = 4,
        label = lab,
    )
end
hline!(p2l, [0.0]; color = :gray, linestyle = :dot, label = "")

p2r = plot(xlabel = L"\Delta", ylabel = L"c", size = (550, 400))
cvs = [-12 * fefit[Δ][1][3] / π for Δ in Δlis]
cverrs = [12 * fefit[Δ][2][3] / π for Δ in Δlis]
plot!(p2r, Δlis, cvs; yerror = cverrs, marker = :diamond, markersize = 4, label = L"cv_{\rm eff}\ (-12d/\pi)")
cs = [cfit[(Δ, :even)] for Δ in Δlis]
plot!(p2r, Δlis, first.(cs); yerror = last.(cs), marker = :circle, markersize = 4, label = L"c\ (\rm entanglement,\ even)")
hline!(p2r, [0.0]; color = :gray, linestyle = :dot, label = "")
p2 = plot(p2l, p2r; layout = (1, 2), size = (1100, 400), left_margin = 5Plots.mm)
savefig(p2, joinpath(@__DIR__, "figs", "central_charge_ind$(ind).png"))

# fig 3: free energy conformal scaling
p3l = plot(xlabel = L"L", ylabel = L"F(L)", size = (550, 400), legend = :bottomright)
p3r = plot(xlabel = L"1/L", ylabel = L"F - fL - e", size = (550, 400), legend = :bottomleft)
Lfine = range(minimum(Llis) - 0.5, maximum(Llis) + 0.5, length = 100)
for Δ in showΔ
    Fs = [FE[(L, Δ)] for L in Llis]
    plot!(p3l, Llis, Fs; marker = :circle, markersize = 4, label = "Δ=$Δ")
    f, e, d = fefit[Δ][1]
    plot!(p3l, Lfine, fe_model(Lfine, [f, e, d]); label = "", linestyle = :dash)
    res = Fs .- (f .* Llis .+ e)
    plot!(p3r, 1 ./ Llis, res; marker = :circle, markersize = 4, label = "Δ=$Δ")
    xinv = range(0, 1 / minimum(Llis) + 0.01, length = 50)
    plot!(p3r, xinv, d .* xinv; label = "", linestyle = :dash)
end
p3 = plot(p3l, p3r; layout = (1, 2), size = (1100, 400), left_margin = 5Plots.mm)
savefig(p3, joinpath(@__DIR__, "figs", "free_energy_scaling_ind$(ind).png"))

println("\nFigures saved to $(joinpath(@__DIR__, "figs"))")
