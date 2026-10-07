# Fibonacci string-net PEPS

The main entry point now builds a two-dimensional PEPS directly. It does not
enumerate configurations or construct a dense wavefunction.

```julia
include("exm/StringNetCondensate/peps.jl")

psi = fibonacci_stringnet_peps(6, 4)                    # OBC
cylx = fibonacci_stringnet_peps(6, 4; bc=:cylinder_x)   # x PBC, y OBC
cyly = fibonacci_stringnet_peps(6, 4; bc=:cylinder_y)   # x OBC, y PBC

psi.tensors         # Vector{ITensor}: local vertex tensors
psi.links           # virtual Index for each honeycomb edge
psi.physical_sites  # one dimension-two physical Index per edge
psi.physical_edges  # physical edge IDs owned by each tensor
psi.lattice.edges   # endpoint vertex IDs for every edge
```

`Lx` and `Ly` count **hexagonal plaquettes**, in the two axial directions of a
honeycomb rhombus. Thus `(1,1)` is one hexagon. Construction time and storage
are $O(L_xL_y)$ at fixed bond dimension. No tensor-network package beyond the
repository's existing `ITensors` dependency is needed. The result is a finite
honeycomb PEPS with explicit graph connectivity, rather than a rectangular
matrix of square-lattice tensors. It is not a `PEPSKit` object.

The physical Hilbert space has one qubit on each honeycomb edge:
$|0\rangle$ is vacuum and $|1\rangle$ is $\tau$. To avoid duplicating physical
qubits, all edges are assigned to their endpoint on the A sublattice. An A
vertex tensor has up to three physical qubit indices; a B vertex tensor has
only virtual indices (equivalently physical dimension one). `physical_edges`
makes this convention explicit. Each physical Index appears on exactly one
tensor, and each virtual Index appears on exactly two.

The five bulk virtual states are the admissible triples

$$
(i,a,b)\in\{(0,0,0),(0,1,1),(1,0,1),(1,1,0),(1,1,1)\},
$$

where $i$ is the physical edge label and $a,b$ are the auxiliary labels of the
faces on the left and right of the edge. Thus the bulk bond dimension is
$D=5$. Fixing exterior face labels to vacuum reduces boundary bonds to $D=2$.
`bond_basis[e]` records the actual ordering for each edge, directed from the
first to the second endpoint in `lattice.edges[e]`.

For three counterclockwise edges $i,j,k$, write the adjacent face labels
$a,b,c$ in the wedges $(k,i),(i,j),(j,k)$. The local tensor weight is

$$
T_v=(d_i d_j d_k)^{1/4}(d_a d_b d_c)^{1/6}
G(i,j,k;a,b,c),\qquad
G(i,j,k;a,b,c)=\frac{[F^{ijc}_{a}]_{kb}}{\sqrt{d_kd_b}},
$$

with $d_0=1$ and $d_1=\phi$. The shared virtual indices enforce agreement of
physical labels and face labels. Every internal hexagon contributes six
factors $d_a^{1/6}$, producing the required loop weight $d_a$.
This is the triple-line construction of
[Schotte et al., arXiv:1909.06284, Appendix A](https://arxiv.org/html/1909.06284#A1),
with each physical edge copied out only once. Missing edges at degree-two
boundary vertices have label zero; exterior auxiliary labels are also zero.

For both OBC and the cylinder closures implemented here, the unnormalized state
is

$$
|\Psi_{\mathrm{raw}}\rangle=
\prod_p(B_p^0+\phi B_p^\tau)|0\cdots0\rangle,
\qquad \langle\Psi_{\mathrm{raw}}|\Psi_{\mathrm{raw}}\rangle
=(1+\phi^2)^{L_xL_y}.
$$

The default `normalize=true` distributes the corresponding normalization among
local tensors. It requires no contraction. Pass `normalize=false` to obtain
amplitudes in the convention where the empty configuration has amplitude one.
On a cylinder, this selects the state obtained by applying the bulk plaquette
projectors to the vacuum, with both exterior auxiliary boundaries fixed to
zero and no noncontractible MPO insertion. It is **not** the spherical closure
used by the older ladder evaluator and does not return all topological sectors.
The periodic direction requires at least two plaquettes; the open direction
may have length one. Doubly periodic torus geometry is not implemented.

The integer vertex coordinates embed as $(u/2,\sqrt{3}v/2)$. The two plaquette
translations are $(3,1)$ and $(0,2)$ in these coordinates. `cylinder_x` identifies
$(u,v)\sim(u+3L_x,v+L_x)$, and `cylinder_y` identifies
$(u,v)\sim(u,v+2L_y)$. The metadata retain the local edge vectors across each seam.

Small-system contractions are provided for verification:

```julia
small = fibonacci_stringnet_peps(2, 1)
config = zeros(Int, length(small.lattice.edges))
a = peps_amplitude(small, config)
state = contract_peps(small)  # ITensor, open physical indices in physical_sites
```

Both helpers use exact contraction and guard intermediate tensor sizes through
`max_elements`. `contract_peps` also guards the final dense state size. They are
not scalable two-dimensional contraction algorithms; large PEPS construction
remains cheap even when exact contraction is impractical. ITensors may warn
about the number of open indices when explicitly contracting a dense state.

Run the dedicated PEPS checks with:

```sh
julia --startup-file=no --project=. exm/StringNetCondensate/checks/test_StringNetPEPS.jl
```

The tests compare full OBC and cylinder states with an independent implementation
of plaquette loop insertion, verify $B_p|\Psi\rangle=|\Psi\rangle$, check
normalization and the unique ownership of physical indices, and compare OBC
amplitudes with the diagram evaluator. Larger patches are checked using selected
amplitudes and reversed plaquette-projector order. These checks live in
`exm/StringNetCondensate/checks/` and run separately from the core package test suite.

# Comparing Born trajectories with a deformed PEPS

`trajectory_equivalence.jl` implements the finite-tension and strict-projection
tests in `fibonacci_peps_stringnet_equivalence_lecture_note.md`. It contains no
module and introduces no package dependency. Given a normalized base reference,

$$
A_J(g)=Z_J^{-1/2}e^{\sum_{e\in Z}J_e\delta_{g_e,0}}A_{\rm SN}(g),\qquad
Z_J=\sum_g |A_{\rm SN}(g)|^2 e^{2\sum_{e\in Z}J_e\delta_{g_e,0}}.
$$

Probability weights contain **twice** the amplitude tension. With all selected
$J_e=-\infty$, `p_zigzag_tau` is the joint probability of the selected edges
being tau, not the product of their individual probabilities.

```julia
include("exm/StringNetCondensate/trajectory_equivalence.jl")
psi = fibonacci_stringnet_peps(1, 1)
reference = exact_stringnet_reference(psi)  # reusable across samples and J
Z = [1, 2]                               # illustrative, user-selected edges
deformed = deform_stringnet_peps(psi; edges=Z, J=-0.4)
deformed.peps                            # actual locally filtered PEPS

# Executable smoke example; replace these arrays with independent Born records.
records = [ones(Int, length(psi.links)) for _ in 1:20]
result = compare_trajectories(reference, records; edges=Z, J=-Inf,
                              rng=MersenneTwister(42))
result.summary
result.rows                              # counts, empirical p, target q, log q

fit = fit_string_tension(reference, records; edges=Z,
                         J_grid=vcat(-Inf, collect(-3:0.1:0)))
```

Local filtering scales with PEPS size and does not contract the network. Exact
reference construction stores $2^E$ amplitudes for $E$ physical edges, with
`max_elements=1<<20` by default; contraction intermediates are also guarded.
The reference is built once. All trajectories are then counted with equal
weight, and the tension fit uses counts of distinct records. These routines
are small-system checks; they do not provide a boundary-MPS contraction for a
full production spacetime lattice.

`TrajectoryEdgeMap(psi, record_edges; fixed_edges, outcome_labels)` explicitly
maps each raw record entry to one PEPS edge. The integer array `record_edges`
must have the same shape as the record. Entries refer to `psi.lattice.edges`;
unrecorded, unfixed edges are summed over in $q(s)$. `fixed_edges` instead
conditions on an event, whose target probability is returned separately.
The sampled source ensemble must have the same conditioning. A zigzag projector
can act on unrecorded edges through `edges=Z, J=-Inf`; it is part of the target
deformation, distinct from subsequent boundary conditioning.

The summary reports TV, Jensen–Shannon divergence, empirical classical fidelity,
forbidden records, known local branching violations, and unseen target mass.
`tv_null_pvalue` calibrates the empirical TV against independent samples of the
same size drawn from the exact target. Histogram fidelity has finite-sample
bias and is not quantum fidelity. Use independent held-out trajectories for
this calibration after fitting $J$. For correlated windows from one trajectory,
set `iid=false`: the IID calibration, standard errors, and unbiased-fidelity
estimates are then disabled.

If the **normalized probability of exactly the compared record** is available,
pass `log_probabilities`. The interface also estimates

$$
\sum_s\sqrt{p(s)q(s)}
=\mathbb E_{s\sim p}\sqrt{q(s)/p(s)},\qquad
D_{\rm KL}(p\Vert q)=\mathbb E_{s\sim p}\log[p(s)/q(s)].
$$

`sampled_target_mass` estimates $\sum_{s:p(s)>0}q(s)$, not an unknown PEPS
normalization. Importance estimates can have large variance if $p$ poorly
covers $q$; ESS and standard errors help assess this but cannot reveal every
unsampled rare event. Do not multiply Born sample counts by $p(s)$ again.

For complete physical records with known scalar coherent amplitudes, additionally
pass their unit complex `phases` to estimate the overlap. Alternatively,
`compare_exact_amplitudes(reference, amplitude_dict; edges=Z, J=...)` compares
full vectors, including relative phases and a global-phase-aligned error.
The dictionary uses physical-label tuples; missing keys assert zero amplitude.
Sampling alone supplies no phases. In particular, tracing the retained final
chain in a monitored circuit generally produces a mixed record state; assigning
$\sqrt{p(s)}$ to every record does not reconstruct its coherent wavefunction.

## Existing CFT data on hpc2ust

Read-only inspection of
`exact/dense_eigensolver/L8/gammaind10/periods32_trajectory_seed4695.jld2`
under `exm/data/Bulk_measure/monitored_dynamics_tci/` found:

* `sample`: `BitMatrix`, shape `(64,4)`;
* `sample_free_energy`: 64 `Float32` layer values;
* `initial_state="TCI_GS"`, `gamma=0.95`, `tau=1.8317808230648227`;
* initial/final Y expectation and entropies, but no final state or amplitude phases.

The generating circuit uses PBC, even sites in odd layers and odd sites in even
layers. Raw `false` means tau and `true` means vacuum: use
`outcome_labels=(1,0)`, reversing the PEPS encoding. At finite measurement
strength these are Kraus outcome labels. The derived spacetime mapping places
them on horizontal rungs with a non-diagonal physical filter; see
[the exact recoupling](honeycomb_recoupling.md).

`cft_trajectory_record` accepts the dictionary returned by `JLD2.load`, so the
comparison implementation does not itself depend on JLD2 or access the cluster:

```julia
using JLD2
# paths: local existing JLD2 files from a single ensemble, with distinct seeds.
records = [cft_trajectory_record(JLD2.load(path); layers=1:2) for path in paths]

# Supply a derived embedding; shape must be (2,L/2) for the prefix above.
mapping = TrajectoryEdgeMap(psi, record_edges; outcome_labels=(1,0),
                           fixed_edges=boundary_labels)
result = compare_trajectories(reference, records; mapping, edges=zigzag_edges,
                              J=-0.4)
```

For all columns of a chronological prefix beginning at layer 1,
`record.log_probability = -sum(sample_free_energy[prefix])` is available.
A spacetime crop, a suffix, or a subset of columns returns `nothing`: summing
their layer free energies would give the wrong marginal probability. If the
source is further postselected to match fixed boundary labels, these saved
probabilities must also be renormalized before use as `log_probabilities`.
Their numerical accuracy is limited by Float32 storage and, for MPS data,
truncation. `terminal_half_strength=true` exposes the generating code's current
last-layer $\tau/2$ convention; the files do not save that switch, so the adapter
does not assume its value silently.

The production TCI initial boundary and summed final-chain boundary differ from
the vacuum-capped PEPS above. The exact correspondence is now implemented in
`honeycomb_recoupling.jl`: outcomes map to horizontal rungs, zigzag edges are
pinned tau, and explicit boundary maps attach the TCI state and decode the final
chain. Crucially, the projective circuit also carries a diagonal rung fugacity;
the bare zigzag projection alone does not reproduce it in this edge basis.
See [the derivation and executable interfaces](honeycomb_recoupling.md).

## Periodic cylinder reference and topological sectors

Our convention is **space around the circumference, measurement time along the
axis**, with an initial boundary vector and an open final fusion-path index.
Spatial PBC does not close time into a torus. `cylinder_transfer.jl` builds an
independent F-symbol transfer calculation in this geometry. Its local projector
at fusion-path site $i$ is

$$
[P_i^{\mathbf1}]_{b'b}
=[F^{a\tau\tau}_d]_{b',\mathbf1}
 [F^{a\tau\tau}_d]^*_{b,\mathbf1},\qquad
a=x_{i-1},\ b=x_i,\ d=x_{i+1},
$$

with periodic neighbors. It uses the repository Hamiltonian only to prepare
the TCI boundary vector; the F-symbol projectors and loop operator are built
independently and checked against the existing measurement implementation.

```julia
include("exm/StringNetCondensate/cylinder_transfer.jl")
boundary = fibonacci_cylinder_reference(8; tau=atanh(0.95))
cylinder_sector_weights(boundary)           # Y=phi, vacuum_flux ≈ 1
# datasets = [JLD2.load(path) for path in paths]
# audit = audit_cft_trajectories(boundary, datasets)

# Actual inspected prefix; tau/2 belongs to layer 64, not the prefix's layer 2.
prefix = Bool[1 0 0 0; 1 0 0 1]
replayed = cylinder_record_probability(boundary, prefix; total_layers=64)
replayed.layer_log_probabilities

tiny = fibonacci_cylinder_reference(4; tau=atanh(0.95))
joint = cylinder_record_state(tiny, 2)
joint.record_purity                        # diagnoses a mixed record marginal
postselected = cylinder_record_state(tiny, 2; final_state=tiny.initial_state)
postselected.amplitudes                    # coherent amplitudes, with phases
postselected.postselection_probability
```

The inspected L=8 directory contains 10,000 trajectories. For seed 4695 the
first two saved log probabilities are `[-4.3262553215, -3.1630616188]`;
independent F contraction gives `[-4.3262555119, -3.1630616562]`, within the
Float32 storage precision. The boundary energy agrees to $3\times10^{-15}$.
This is a check of the circuit-to-F-network implementation, **not** a fidelity
measurement against the honeycomb PEPS. No production ensemble was recomputed.
The transfer calculation is linear in trajectory length at fixed width; sparse
local projectors are reused. Exact boundary ED and the dense Y operator are
guarded by `max_basis=512` before basis enumeration. Full coherent record
enumeration has a separate `max_elements` guard.

The loop around the spatial circle measures flux through the cylinder. Its
Fibonacci fusion algebra gives

$$
Y^2=I+Y,\qquad
P_{\mathbf1}=\frac{Y+\phi^{-1}I}{\phi+\phi^{-1}},\qquad
P_\tau=\frac{\phi I-Y}{\phi+\phi^{-1}}.
$$

These formulas and the eigenvalues $\phi,-\phi^{-1}$ agree with
[Buican and Gromov, arXiv:1701.02800, Sec. 5.1](https://arxiv.org/html/1701.02800#S5.SS1).
`fibonacci_flux_projectors` and `cylinder_sector_weights` resolve these two
boundary sectors; `sector=:tau_flux` selects the lowest TCI-Hamiltonian state
in the other sector. Every local Kraus operator commutes with Y, so an input
in a single sector remains there in each trajectory. A superposition of sectors
can have its sector weights updated by conditioning on an outcome.

Two Y eigenvalues do not by themselves specify a logical qubit or fully label
the doubled theory. The doubled-Fibonacci bulk sectors are
$(\mathbf1,\mathbf1)$, $(\tau,\mathbf1)$, $(\mathbf1,\bar\tau)$,
and $(\tau,\bar\tau)$, resolved by the tube algebra's central idempotents;
see [Bultinck et al., arXiv:1511.08090, Appendix D.1](https://arxiv.org/html/1511.08090#A4.SS1).
Which sectors survive on an open cylinder depends on its end boundaries.
The two polynomials in Y implemented here do not construct those four PEPS
idempotents. `honeycomb_boundary_symmetry` now transports them to the PEPS
virtual cut, including its boundary metric. The TCI input with Y=phi differs
from the alternating fusion-path caps reproducing the original
`fibonacci_stringnet_peps` closure.

There is also a necessary **physical-basis check** at finite strength. Writing
$h=(1+e^{-2\tau})^{-1/2}$ and $\ell=e^{-\tau}h$, the actual circuit uses

$$
K_{\mathbf1}=hP_{\mathbf1}+\ell P_\tau,\qquad
K_\tau=\ell P_{\mathbf1}+hP_\tau.
$$

Consequently its coherent outcome tensor is obtained from the projective
one by the non-diagonal physical filter

$$
R(\tau)=\begin{pmatrix}h&\ell\\\ell&h\end{pmatrix}
=(h+\ell)U\begin{pmatrix}e^J&0\\0&1\end{pmatrix}U^\dagger,
\qquad U=\frac1{\sqrt2}\begin{pmatrix}1&1\\-1&1\end{pmatrix},
\qquad J=\log\tanh\frac\tau2.
$$

`fibonacci_outcome_filter` and `fibonacci_outcome_deformation` implement this
identity. The checks compare the **whole coherent record state**, retaining the
final chain, with the product of these physical filters. The honeycomb
recoupling fixes the full physical rung filter to
`R(tau)*Diagonal([1,sqrt(phi)])/phi`. Diagonalizing R also rotates the underlying
projective state and the measurement basis; that rotated J is not the note's
zigzag-edge tension. Raw Born samples in the original basis cannot be converted
to samples in the rotated basis by a deterministic label reversal.

Run the deformation, statistics, phase-sensitivity, and data-schema checks with:

```sh
julia --startup-file=no --project=. exm/StringNetCondensate/checks/test_TrajectoryEquivalence.jl
```

# Small-system periodic-ladder reference

`stringnet.jl` constructs the normalized state in the orthonormal edge-label
basis, with labels `0` (vacuum) and `1` ($\tau$). The cyclic ordering at each
trivalent vertex is part of the diagram data.

The amplitude is the spherical evaluation of the labeled diagram. Equivalently,
the two holes of the annular ladder are capped by vacuum disks: a winding $\tau$
loop has amplitude $\phi$ relative to the empty configuration. This specifies the
closure; it does not implement arbitrary annular boundary conditions or other
flux sectors. The local data are the unitary Fibonacci F-symbols, a $\tau$-loop
value $\phi$, and a $\tau$-digon factor $\sqrt{\phi}$. In particular, the one-dimensional
block satisfies $[F^{\tau\tau\tau}_{0}]_{\tau\tau}=1$.

References: Lin, Levin and Burnell, [arXiv:2012.14424, Eq. (103)](https://arxiv.org/html/2012.14424#S7.SS4);
Minev et al., [arXiv:2406.12820, Appendix C](https://arxiv.org/html/2406.12820v2#A3).
The latter supplies the independent check

$$\operatorname{Eval}(G)^2 = \frac{\chi(\widehat{G},\phi+2)}{\phi+2}.$$

```julia
include("exm/StringNetCondensate/stringnet.jl")
lad = Ladder(8)
configs, amplitudes = stringnet_ground_state(lad)
S = arc_entropy(configs, amplitudes, lad, 4)
```

`arc_entropy` uses natural logarithms and the edge tensor-product partition
$A=\{t_i,b_i,v_i:1\leq i\leq\ell\}$. This is not a union of complete
plaquettes. The one-column cuts have overlapping boundary vertices; they need
not share the entropy plateau of longer arcs. For $L=4,6,8$ the one-column
entropy is 1.7195478808 and the interior plateau is 2.3004047708. The function
accepts normalized real or complex amplitudes and sums duplicate configurations
coherently. It works for general edge states, not just fusion-admissible ones.

Configuration enumeration uses four transfer states with closure pruning and
copies only completed configurations. Ordering is deterministic but differs
from the earlier dictionary-based enumeration; always keep configurations
paired with their returned amplitudes. Diagram evaluation reuses the ladder
geometry and binary canonical keys. Memo keys are internal and should not be
persisted across implementation versions. Entropy uses independent blocks of
the sparse Schmidt matrix and SVD, with the smaller subsystem on the row side.
There is still exponential cost in enumerating the full state; this is a
small-system exact calculation, not an MPS solver.

Run the regression suite and larger-ladder checks with:

```sh
julia --startup-file=no --project=. exm/StringNetCondensate/checks/test_StringNetCondensate.jl
```

The tests cover F-block unitarity, the pentagon identity, brute-force basis
counts, independent chromatic-polynomial weights, every edge F-move for L=3,
reduction-path independence, sphere normalization, graph relabeling, and
entropy against dense SVD for complex random states. They live in
`exm/StringNetCondensate/checks/` and run separately from the core package test
suite. A failing test makes the script exit unsuccessfully.

Local Julia 1.12.2 measurements, including fresh memo construction each run:

| Operation | Before | After | Allocated before / after |
| --- | ---: | ---: | ---: |
| L=8 enumeration, 29,375 configurations | 30.42 ms | 1.46 ms | 89.2 / 8.06 MiB |
| L=8 ground state | 1972.53 ms | 136.02 ms | 3328.8 / 215.7 MiB |
| L=6 entropy, ell=5 | 110.98 ms | 0.29 ms | 14.50 / 0.40 MiB |
| L=8 entropy, ell=4 | 11.44 ms | 13.56 ms | 5.87 / 3.18 MiB |

The balanced-cut SVD trades some time for lower allocation and avoids forming
a density matrix. The improvement is largest for unequal partitions. These
allocations are cumulative Julia allocations, not peak resident memory.
Matching configurations gives a maximum normalized-amplitude difference of
$2.6\times10^{-15}$ for $L=6,8$; the F-symbol bug was bypassed in the old ladder evaluator by
its vacuum-absorption preprocessing, but affected the F-symbol API itself.
