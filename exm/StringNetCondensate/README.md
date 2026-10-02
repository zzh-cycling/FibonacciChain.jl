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
