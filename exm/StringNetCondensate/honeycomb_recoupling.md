# Exact recoupling on a spatial cylinder

The operational comparison is between the coherent trajectory tensor
$W_{x_f,s}=\langle x_f|K_s|\psi_0\rangle$ and a honeycomb PEPS with the same
initial boundary and an open final fusion path. Spatial PBC closes the
circumference; measurement time runs along the cylinder axis. The Born weight
is $p(s)=\sum_{x_f}|W_{x_f,s}|^2$. A final vector instead defines a postselected
scalar amplitude. These two boundary contractions must not be interchanged.

`honeycomb_recoupling.jl` implements this correspondence with actual local
ITensors. Its independent sparse transfer implementation contracts paired
G vertices for many trajectories without enumerating the record Hilbert space.
No additional module or package dependency is introduced.

## Local identity and the physical rung filter

Use the existing PEPS convention $0=\mathbf1$, $1=\tau$, $d_0=1$, $d_1=\phi$.
At a measurement site let $a=x_{i-1}$, $b=x_i$, and $d=x_{i+1}$ be fusion-path
labels. The projective fusion outcome $s$ has matrix element

$$
[P_i^s]_{b'b}=[F^{a\tau\tau}_d]_{b',s}[F^{a\tau\tau}_d]^*_{b,s}.
$$

In our real unitary gauge, a honeycomb vertex with physical labels
$(\tau,\tau,s)$ and adjacent face labels $(d,b,a)$ is

$$
v_s(a,b,d)=(\phi^2d_s)^{1/4}(d_a d_b d_d)^{1/6}
G(\tau,\tau,s;d,b,a),\qquad
G(\tau,\tau,s;d,b,a)=\frac{[F^{a\tau\tau}_d]_{b,s}}{\sqrt{d_b d_s}}.
$$

Resolving one measurement into its lower and upper trivalent vertices gives

$$
v_s(a,b',d)v_s(a,b,d)
=\frac{\phi}{\sqrt{d_s}}
\frac{(d_a d_d)^{1/3}}{(d_b d_{b'})^{1/3}}
[P_i^s]_{b'b}.
$$

This identity is checked for every admissible local label. The quantum-dimension
factor depends on the **physical outcome** and cannot be removed as a common
normalization or a virtual gauge. Reversing the identity requires the rung
filter

$$
C=\frac{1}{\phi}\operatorname{diag}(1,\sqrt\phi)
=\phi^{-1/2}\operatorname{diag}(e^{J_{\rm rung}},1),\qquad
J_{\rm rung}=-\frac12\log\phi.
$$

At finite strength, the actual Kraus operators additionally mix projective
outcomes through

$$
R(\tau)=\frac{1}{\sqrt{1+e^{-2\tau}}}
\begin{pmatrix}1&e^{-\tau}\\e^{-\tau}&1\end{pmatrix},\qquad
M(\tau)=R(\tau)C.
$$

Thus the monitored state is the tau-pinned honeycomb G network with **M on each
rung**, including its two boundary rung rows. It is not, in this physical edge
basis, merely the bare zigzag-projected state. The exact projective limit is
$M(\infty)=C$. At finite strength M is non-diagonal; diagonalizing R alone also
changes the projective reference state and does not turn raw Born records into
samples in that rotated basis. The terminal layer uses $M(\tau/2)$ when the
generating code's half-strength convention is enabled.

## Boundary maps and flux symmetry

In a complete layer, each unmeasured face appears as an outer neighbor twice.
Define a diagonal matrix on the periodic fusion-path basis by

$$
q_r(x)=\prod_{i\notin I_r}d_{x_i}^{1/3}
       \prod_{i\in I_r}d_{x_i}^{-1/3},\qquad
Q_r=\operatorname{diag}q_r(x),
$$

where $I_r$ contains even sites in odd layers and odd sites in even layers.
The paired-G layer with the physical filter M obeys

$$
\mathcal G_r(s_r)=Q_r K_r(s_r)Q_r.
$$

The next layer swaps measured and unmeasured sites, so $Q_{r+1}Q_r=I$.
Consequently all interior factors cancel:

$$
\mathcal G_T(s_T)\cdots\mathcal G_1(s_1)
=Q_T K_s Q_1.
$$

The input PEPS boundary is therefore $Q_1^{-1}|\psi_0\rangle$; the decoded
output is $Q_T^{-1}|v_{\rm out}\rangle$. A final boundary bra is
$\langle f|Q_T^{-1}$. Summing over final states means using the boundary metric
$Q_T^{-2}$, not the Euclidean norm of the undecoded virtual tensor. The input
and output are encoded on the cut zigzag bonds by
$x\mapsto\{(x_{i-1},x_i)\}_{i=1}^L$, with periodic indices. This construction
accepts arbitrary complex boundary vectors, including superpositions of sectors.

The same maps transport the chain loop operator onto the PEPS cut:

$$
Y_{\rm in}=Q_1^{-1}YQ_1,\qquad
Y_{\rm out}=Q_TYQ_T^{-1}.
$$

They satisfy $Y_{\rm cut}^2=I+Y_{\rm cut}$ and are Hermitian in the respective
metrics $Q_1^2$ and $Q_T^{-2}$. The projectors

$$
P_{\mathbf1}=\frac{Y_{\rm cut}+\phi^{-1}I}{\phi+\phi^{-1}},\qquad
P_\tau=\frac{\phi I-Y_{\rm cut}}{\phi+\phi^{-1}}
$$

therefore explicitly identify the chain's two Y sectors on the virtual PEPS
boundary. This does not add the missing tube-algebra operators needed to resolve
all four doubled-Fibonacci bulk anyons. It also does not identify the two
sector spaces with a two-dimensional logical Hilbert space.

## Relation to the existing finite honeycomb patch

For T measurement layers on L worldlines, the interior graph is
`HoneycombPatch(T-2,L÷2;bc=:cylinder_y)`, provided $T\geq3$. Here T is the number
of **layers**, not periods. The first and last measurement rows supply dangling
boundary rungs. Removing their outer trivalent vertices leaves exactly the
degree-two boundary vertices of the original patch.

For $2\leq r\leq T-1$, outcome $(r,c)$ maps to the horizontal edge centered at

$$
u=3(r-2),\qquad v=i\pmod L,\qquad
i=\begin{cases}2c,&r\text{ odd},\\2c-1,&r\text{ even}.\end{cases}
$$

All slanted physical edges are zigzag worldlines and have label tau.
`honeycomb_record_geometry` returns these exact edge IDs; its first/last rows
contain zero because those outcomes belong to boundary rungs, not bulk edges.

The original PEPS fixes the exterior auxiliary faces and the missing boundary
physical edges to vacuum. To recover that particular closure, take the
alternating fusion path with $x_i=\tau$ on measured sites and $x_i=0$ on
unmeasured sites at each end, and fix both boundary outcome rows to vacuum.
`honeycomb_vacuum_cap` constructs these vectors. With the Q maps above,

$$
A_{\rm raw\ G}(s_{\rm boundary}=0,s_{\rm interior})
=\phi^{L/2}\,
A_{\rm original\ PEPS}(s_{\rm interior},Z=\tau),
$$

where the original PEPS uses `normalize=false`. The scalar comes from the
removed end vertices and boundary maps. Tests compare **every** interior
amplitude, including its sign, for $(L,T)=(4,3),(4,4),(6,3)$. Applying C on all
rungs then gives the circuit amplitude exactly. Omitting C yields fidelity
strictly below one even with these matched caps. A TCI ground-state boundary
is a different vector and is attached through Q directly, without replacing
it by the alternating cap. For L=8 the alternating cap has Y expectation zero,
with sector weights 0.2763932023 and 0.7236067977, whereas the TCI input has
Y=phi. This explicitly distinguishes the old vacuum-cap convention from the
chain's vacuum-flux eigenstate.

With matched alternating caps, the normalized fidelities between the bare
zigzag-projected state and the projective circuit state are 0.9748801606 for
$(L,T)=(4,3)$ and 0.9446600726 for $(4,4)$. Including C restores the complete
complex amplitude vector, without fitting a phase or normalization factor.

For a data-level audit, ten existing L=8, gamma=0.95, 32-period trajectories
were read from hpc2ust (seeds 1, 10, 100, 1000, 10000, 1001--1005). Each has
64 layers and 256 outcomes. Paired-G transfer and independent circuit replay
agree in whole-record log probability to $1.8\times10^{-13}$; the maximum
difference from the saved Float32 free-energy sums is $1.2\times10^{-6}$.
The per-record results are saved in
`exm/data/StringNetCondensate/recoupling_audit/L8_gamma095_exact_10_records.csv`.
This is an implementation audit of the derived identity, not an inference
from ten samples about arbitrary PEPS states. No trajectories were regenerated.

## Using the implementation

```julia
include("exm/StringNetCondensate/honeycomb_recoupling.jl")

# Full coherent amplitude comparison on a small spatial cylinder.
boundary = fibonacci_cylinder_reference(4; tau=atanh(0.95), sector=:vacuum_flux)
psi = recoupled_honeycomb_peps(4, 3; tau=boundary.tau)
g = contract_recoupled_honeycomb(psi, boundary)
f = cylinder_record_state(boundary, 3)
@assert isapprox(g.joint_amplitudes, f.joint_amplitudes; atol=1e-12)

# The bare G network, with the same TCI boundary and zigzag projection.
bare = recoupled_honeycomb_peps(4, 3; representation=:stringnet)
q = contract_recoupled_honeycomb(bare, boundary)
q.probabilities ./ q.norm2  # normalize the bare candidate explicitly

geometry = honeycomb_record_geometry(4, 3)
geometry.record_edges
geometry.zigzag_edges
honeycomb_boundary_symmetry(boundary, 3; side=:output)

# Many existing records at fixed L: reuse sparse G-vertex operators.
using JLD2
boundary8 = fibonacci_cylinder_reference(8; tau=atanh(0.95))
cache = HoneycombTransferCache(boundary8)
# paths must identify records from the same initial-state/parameter ensemble.
# rows = compare_honeycomb_trajectories(cache, (JLD2.load(p) for p in paths))
```

The last interface compares the actual G-network log probability with saved
`sample_free_energy`; it does not estimate a wavefunction phase from sample
frequencies. For a prefix, pass the original `total_layers` to
`honeycomb_record_probability` so that the terminal half-strength filter is not
applied prematurely. `representation=:stringnet` returns an unnormalized
weight, deliberately distinguished from a normalized probability.

PEPS construction has linear size at fixed local dimensions (zigzag bond
dimension 3, rung bond dimension 5). Exact open-boundary contraction and coherent
record enumeration remain guarded exponential calculations. The cached sparse
G transfer is linear in record length at fixed width; its fusion-path dimension
still grows exponentially with circumference. The present work does not supply
an MPS approximation for very wide cylinders.

Run all recoupling checks with:

```sh
julia --startup-file=no --project=. exm/StringNetCondensate/checks/test_HoneycombRecoupling.jl
```

The new identities follow algebraically from the local tensor convention and
are numerically tested, rather than inferred from a similarity of probability
histograms. Background on the boundary Fibonacci loop and its two eigenvalues
is given in [arXiv:1701.02800, Sec. 5.1](https://arxiv.org/html/1701.02800#S5.SS1);
the distinction from the four doubled-Fibonacci tube sectors is explained in
[arXiv:1511.08090, Appendix D.1](https://arxiv.org/html/1511.08090#A4.SS1).
