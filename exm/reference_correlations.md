# Reference-qubit correlations and anisotropy

`Born_OBF/corr_calculate.jl` and `HybridEvolution/corr_calculate.jl` use the shared
exact-state implementation in `reference_correlations.jl`. They implement the
reference-copy probe used by `Bulk_measure/corr_calculate.jl`, with model-specific
layer schedules and fresh conditional Born sampling after probe insertion.

## Protocol

Start from the existing circuit initial state: the uniform `|+>` product state
for OBF, or the all-zero coherent fusion path for Fibonacci hybrid circuits.
Evolve `warmup` complete periods. For a spatial sample (`separation=0`), copy
sites `1` and `L/2+1` into two references at that boundary. For a temporal sample,
copy site `L/2+1` twice, separated by `separation` complete periods. Evolve another
`final_periods` after the second insertion, then compute

$$
I(R_1:R_2)=S(R_1)+S(R_2)-S(R_1R_2).
$$

This is reference mutual information in nats, not an ordinary operator two-point
expectation. Copying preserves the original branch amplitudes; it does not reset
the system site or enforce equal Schmidt weights. References subsequently remain
passive. The Born probability is the sum over all reference branches, and the
same sampled Kraus operator acts on every branch. Normalizing or sampling each
branch independently would give a different circuit.

A period has **14 layers for OBF** and **2 layers for Fibonacci hybrid**. All
insertions occur at period boundaries. OBF uses the existing coupling-dependent
`XZZ`, `ZZX`, `X`, and `ZZ` schedule, including the `lambda >= 10` convention in
`Born_OBF/config.jl`. Hybrid uses projective measurements with probability `p`
and independent uniform unitary angles otherwise. The callable hybrid wrapper
also accepts finite `τ`. There is no final half-strength layer and no unitary-only
scrambling stage. These choices are stored in the trajectory metadata.

## Run

From the package root, small smoke runs are:

```sh
julia --project=. exm/Born_OBF/corr_calculate.jl 6 0.4 7 6 6 1,2,3,4 1 10 tmp/obf_reference_scan
julia --project=. exm/HybridEvolution/corr_calculate.jl 8 0.5 8 8 1,2,3,4 1 10 tmp/hybrid_reference_scan
```

OBF arguments are `L lambda tau_idx warmup final_periods dt_csv seed_start
seed_end output_dir`; hybrid arguments omit `lambda tau_idx` and instead take
`p` after `L`. Times are integer periods, **not multiples of L**. The temporal
list must be sorted, unique, and positive; the spatial sample is added
separately. At least two seeds are required for SEM. Use a new output directory;
existing data are never overwritten. These example sizes/times/sample counts
are not convergence recommendations. Production outputs belong under
`exm/data/OBF/` or `exm/data/HybridEvolution/`.

Each `dt*_seed*.jld2` stores probe geometry, initialization, circuit parameters,
reference entropies, MI, `active_mask`, `measurement_mask`, `outcomes`,
`unitary_angles`, and layer-wise negative log Born probabilities. Inactive OBF
columns have `active_mask=false`: they are neither unitary gates nor measurement
outcomes. Angles are only defined at active unitary locations. The full record
includes warmup. `averaged.jld2` stores paired seeds, raw MI arrays, means, SEM,
and the crossing/bootstrap results. Temporal array rows are seeds and columns
are separations. No large final state is saved.

For deterministic replay, load one trajectory as a NamedTuple and pass it as
`record` to the corresponding wrapper with the same parameters. The returned
`state` is available for checks in memory. Merely including either entry script
defines functions and does not start a scan. Run the two entry scripts in separate
Julia sessions to avoid repeatedly including the shared definitions.

## Extract alpha

Under the cylinder-scaling assumption for the probe observable,

$$
\frac{\overline{I(0,\delta t)}}{\overline{I(L/2,0)}}
=\left[\sinh\left(\frac{\pi\alpha\delta t}{L}\right)\right]^{-2\Delta}.
$$

Matching the ensemble means at a crossing $t_*$ eliminates $\Delta$:

$$
\alpha=\frac{\operatorname{asinh}(1)L}{\pi t_*}.
$$

The code uses bracketed linear interpolation, rather than the old unconstrained
exponential-plus-offset fit. It does not extrapolate. Absent or multiple crossings
return `NaN`; extend or refine the temporal grid and improve sampling. Whole-seed
bootstrap preserves covariance across all separations and the shared spatial
baseline. `alpha_sem` and percentile `alpha_ci` are reported only when at least
95% of bootstrap replicates have an unambiguous crossing. All bootstrap alpha
values, including failures, are saved.

Alpha is in sites per complete period; divide by the number of layers per period
when comparing with a coefficient expressed per layer. A crossing alone does not
establish CFT scaling. Check warmup and final-time convergence, system-size
trends, temporal-grid interpolation bias, and the scaling form near the regime
of interest. Bootstrap errors do not include these systematics. In OBF, also
check insertion-origin dependence of the three-site Trotter grouping; lengths
and antipodal separations compatible with that grouping (e.g. L divisible by 6)
make comparisons easier. For hybrid circuits, do not assume the same CFT form
holds throughout the p scan.

The state contains at most four system vectors. This backend is exact and
exponential in system size; it is not an MPS implementation.
