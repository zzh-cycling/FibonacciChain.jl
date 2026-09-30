# Hybrid Evolution Research Workflow

Keep reusable scripts in this directory as top-level Julia scripts. Put one-off sweep and validation scripts under the ignored `tmp/` directory. The package module is defined only in `src/FibonacciChain.jl`; do not introduce another `module` here. Do not add target-specific tests under `test/` or register them in `test/runtests.jl`.

## Paper Sweep Rules

Use periodic Fibonacci chains with $L=8,12,16,20$. The $p$ grid is the sorted set of keys in this directory's `config.jl`: $0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.95$. Evolve every point for exactly $t=4L$ complete periods and record observables at every integer period, including $t=0$. The values returned by `get_hybrid_time_in_L` are separate convergence guidance; this sweep uses the keys as its $p$ grid.

Run 1000 Born trajectories at every $(L,p)$ point. Reuse the identical seed list `1:1000` at every point. Use the exact backend, the all-zero coherent fusion path, projective measurements with $\tau=\infty$, independent measurement locations with probability $p$, independent random unitary angles, and the even-then-odd layer order.

Save the full measurement record by default in every trajectory file: `measurement_mask`, `outcomes`, and `unitary_angles`. Read `outcomes` only where `measurement_mask` is true; false mask entries denote unitary gates, not a measurement outcome. Keep these raw records in trajectory files; the small averaged file contains only statistics.

Save each completed trajectory under `exm/data/HybridEvolution/`, retaining the existing keys and adding `reference_entropy_final`. This is the final ancilla von Neumann entropy for the topological-charge reference-qubit construction in `../Bulk_measure/topological_charge_sharpening.jl`. Sector-preserving hybrid gates keep the reference branches orthogonal, so $q_1=(\langle Y\rangle+1/\phi)/(\phi+1/\phi)$ and $S_{\mathrm{ref}}=-q_1\log q_1-(1-q_1)\log(1-q_1)$. Collect the per-point means, SEM, and sharpening results into one small `averaged.jld2` file for the paper repository. Before using it, verify every point has exactly the seed list `1:1000`, $t=4L$, and 1000 samples. Interpret the output as finite-time data.
