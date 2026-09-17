# Ising anyon chain

This implementation follows [Shao, arXiv:2308.00747v2, Sections 3.4.1–3.4.2](https://arxiv.org/html/2308.00747v2#S3.SS4).
`IsingAnyon()` uses periodic fusion paths. `SpinHalf()` with `model_type=:Ising`
uses the existing spin Hilbert space. The two constructors count different sites.

## 1. Hilbert space

Write the paper's objects ``I,\eta,\mathcal D`` as ``I,\eta,\sigma``:

```math
\eta\times\eta=I,\qquad \eta\times\sigma=\sigma,\qquad
\sigma\times\sigma=I+\eta.
```

For ``L`` physical σ anyons, define

```math
\mathcal H_L=\operatorname{span}\{\lvert x_1\cdots x_L\rangle:
N_{x_i\sigma}^{x_{i+1}}=1,\quad x_{L+1}=x_1\}.
```

For even ``L=2n``, the paths split into σ on odd links or σ on even
links; the remaining labels independently take values ``I,\eta``. Thus

```math
\mathcal H_{2n}\simeq(\mathbb C^2)^{\otimes n}\oplus
(\mathbb C^2)^{\otimes n},\qquad \dim\mathcal H_{2n}=2^{n+1}.
```

Odd ``L`` gives an empty space. `anyon_basis(model)` returns sorted
`NTuple{L,UInt8}` paths, with `0 = I`, `1 = η`, `2 = σ`; tuple index is
the link index. Both sectors are retained. Open boundaries are rejected
because their fusion paths require a choice of boundary charges.

## 2. Hamiltonian

In the ``I,\eta`` order the nontrivial F move is

```math
F^{\sigma\sigma\sigma}_{\sigma}=\frac1{\sqrt2}
\begin{pmatrix}1&1\\1&-1\end{pmatrix}.
```

Let ``P_i^I`` project the two σ anyons adjacent to link ``i`` onto vacuum.
Our local action, obtained by applying this F move, is

```math
P_i^I=\begin{cases}
(1+X_i)/2,&x_i\in\{I,\eta\},\\
(1+Z_{i-1}Z_{i+1})/2,&x_i=\sigma,
\end{cases}\qquad P_i^\eta=1-P_i^I.
```

Here ``X`` exchanges ``I,\eta`` and ``Z`` has eigenvalues ``+1,-1``.
Indices wrap around the ring. `ising_fusion_projector(model, i;
channel=:I)` returns a sparse matrix; use `channel=:eta` for the complement.

The code retains the exact projector normalization of Eq. (3.69):

```math
H=-J\sum_{i\ \mathrm{odd}}P_i^I-h\sum_{i\ \mathrm{even}}P_i^I.
```

`J=h=1` is the uniform chain. On the sector with σ on odd links,
our formula gives, with ``n=L/2`` effective spins,

```math
H\big|_{\sigma\ \mathrm{odd}}=
-\frac J2\sum_{j=1}^n Z_jZ_{j+1}
-\frac h2\sum_{j=1}^n X_j-\frac{n(J+h)}2\,1.
```

The other sector exchanges ``J,h``. Consequently the uniform model consists
of two copies of the critical spin Hamiltonian divided by two and shifted
by ``-L/2``. This normalization matters when comparing numerical energies.
`anyon_ham` and `anyon_ham_sparse` construct dense and sparse matrices.

## 3. Measurement operator

The following measurement protocol is our extension of the repository's
spin-chain convention, using the paper's fusion projectors. Define
``A_i=P_i^I-P_i^\eta=2P_i^I-1``, so ``A_i^2=1``. For strength ``\tau\ge0``
and outcome ``s=\pm1``, use

```math
M_{i,s}(\tau)=\frac{\exp(s\tau A_i/2)}{\sqrt{2\cosh\tau}}
=\frac{e^{s\tau/2}P_i^I+e^{-s\tau/2}P_i^\eta}
{\sqrt{2\cosh\tau}}.
```

These satisfy ``\sum_s M_{i,s}^\dagger M_{i,s}=1``. A normalized input
``|\psi\rangle`` gives probability ``p_s=\|M_{i,s}\psi\|^2`` and conditional
state ``M_{i,s}\psi/\sqrt{p_s}`` for a nonzero-probability outcome.

`sign=false` means ``s=+1`` (vacuum in the projective limit);
`sign=true` means ``s=-1`` (η). At zero strength both operators equal
``1/\sqrt2``. At infinite strength they are ``P_i^I,P_i^\eta``.
The implementation uses bounded eigen-amplitudes so `τ=Inf` is supported.
`measuremap` returns the **unnormalized** state.

```julia
using FibonacciChain, LinearAlgebra

model = AnyonModel(IsingAnyon(), 6; J=1.0, h=1.0)
basis = anyon_basis(model)                   # 16 paths, two 8-dimensional sectors
H = anyon_ham_sparse(model)
P = ising_fusion_projector(model, 2; channel=:I)
ψ = zeros(ComplexF64, length(basis)); ψ[1] = 1
ϕ = measuremap(model, 0.7, ψ, 2, false)
p = norm(ϕ)^2
ψ_conditional = ϕ / sqrt(p)
M = FibonacciChain.measure_matrix(model, 0.7, 2, false)
@assert M * ψ ≈ ϕ
```

This API covers fusion-path bases, Hamiltonians, and local measurements.
MPS, entanglement, symmetry-sector and multilayer evolution interfaces for
this new basis are not implemented.
