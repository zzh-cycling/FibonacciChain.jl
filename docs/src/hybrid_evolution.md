# Hybrid Monitored Evolution

The hybrid evolution routines implement monitored Fibonacci circuits in which
each valid spacetime location contains either a weak measurement or a unitary
fusion-channel gate. Exact-state and Matrix Product State (MPS) backends use the
same circuit geometry and can replay the same stored circuit realization.

## Circuit Setup

The current implementation acts on an even-length Fibonacci chain. One period
contains two staggered layers:

- odd global layers act on sites `2, 4, ..., N`;
- even global layers act on sites `1, 3, ..., N-1`.

Each layer therefore contains `N ÷ 2` events. For

```math
\Delta t=t_2-t_1+1,
```

a complete schedule has size

```math
(2\Delta t)\times (N/2).
```

Rows label staggered layers, while columns enumerate the valid event locations
within a layer. For example, at ``N=6``, columns `1:3` correspond to sites
`2, 4, 6` in the first layer and sites `1, 3, 5` in the second layer.

At each event, the circuit applies either a measurement Kraus operator or

```math
U(\theta)=\Pi_1+e^{i\theta}\Pi_\tau.
```

The alternatives are exclusive: a location is a measurement brick or a
unitary brick, not a unitary followed by a probabilistic measurement.

The Fibonacci measurement outcomes use the convention

- `false`: the ``\tau`` fusion channel;
- `true`: the trivial, or ``1``, fusion channel.

For measurement strength ``\tau_m``, the corresponding Kraus operators are

```math
M_\tau(\tau_m)
=\frac{e^{\tau_m}\Pi_\tau+\Pi_1}{\sqrt{e^{2\tau_m}+1}},
```

```math
M_1(\tau_m)
=\frac{\Pi_\tau+e^{\tau_m}\Pi_1}{\sqrt{e^{2\tau_m}+1}}.
```

They obey

```math
M_\tau^\dagger M_\tau+M_1^\dagger M_1=I.
```

The limit ``\tau_m\to\infty`` gives the projectors ``\Pi_\tau`` and ``\Pi_1``.

## Configuration

A `HybridConfig` embeds a [`MeasureConfig`](@ref) and adds the hybrid
gate probability and unitary-angle settings:

```julia
using FibonacciChain, Random

config = HybridConfig(
    τ = 1.0,
    t₁ = 1,
    t₂ = 20,
    mode = :Born,
    rng = MersenneTwister(42),
    enable_τ_eff = false,
    p = 0.35,
    θ = π,
    random_angles = false,
)
```

The hybrid-specific fields are:

- `p`: probability that an event is a measurement, with ``0\leq p\leq1``;
- `θ`: fixed angle used by unitary bricks when `random_angles=false`;
- `random_angles`: whether each unitary brick draws an independent angle.

The embedded `MeasureConfig` controls the measurement strength, time interval,
sampling mode, random-number generator, and MPS truncation settings. If
`enable_τ_eff=true`, the final layer of the run uses measurement strength
`τ / 2`; all other layers use `τ`.

The same configuration can be constructed from an existing measurement
configuration:

```julia
measure_config = MeasureConfig(
    τ = 1.0,
    t₂ = 20,
    mode = :Born,
    rng = MersenneTwister(42),
)

config = HybridConfig(
    measure_config;
    p = 0.35,
    θ = π,
    random_angles = false,
)
```

## Independent Gate Selection

In `mode=:Born`, every valid event draws a fresh uniform random number

```math
r_{x,t}\sim\operatorname{Uniform}(0,1)
```

and selects

```math
g_{x,t}=
\begin{cases}
\text{measurement}, & r_{x,t}<p,\\
\text{unitary}, & r_{x,t}\geq p.
\end{cases}
```

Consequently, the gate indicators are independent Bernoulli variables,

```math
g_{x,t}\overset{\mathrm{i.i.d.}}{\sim}\operatorname{Bernoulli}(p).
```

The gate choice is independent of the current quantum state, previous
measurement outcomes, and local fusion channel. The implementation is
equivalent to

```julia
is_measurement = p == 1 ? true :
                 p == 0 ? false : rand(rng) < p
```

The cases `p == 1` and `p == 0` are handled without consuming a gate-choice
random number. They give the measurement-only and unitary-only limits,
respectively.

If an event is a measurement, its outcome is then sampled from the current
state. The probability of the ``\tau`` outcome is

```math
q_\tau=\left\|M_\tau|\psi\rangle\right\|^2,
```

and a second random draw selects between ``M_\tau`` and ``M_1``. Gate locations are
therefore state-independent, but measurement outcomes generally are neither
independent nor identically distributed: earlier events change the state on
which later Born probabilities are evaluated.

In `mode=:sample`, no gate location is drawn. The routine reads the Boolean
gate mask from a supplied [`HybridGateSchedule`](@ref), and `config.p` does not
affect the replayed circuit.

## Unitary-Angle Sampling

For `random_angles=false`, every unitary brick uses the configured angle,

```math
\theta_{x,t}=\theta.
```

For `random_angles=true`, an angle is drawn only after the event has been
selected as a unitary:

```math
\theta_{x,t}=2\pi u_{x,t},
\qquad
u_{x,t}\sim\operatorname{Uniform}(0,1).
```

Conditional on the unitary locations, the stored angles are therefore
independent and uniformly distributed on ``[0,2\pi)``. Measurement locations do
not draw an angle and store `NaN` in the angle array.

The exact-state backend applies the gate without constructing a dense matrix:

```math
U(\theta)|\psi\rangle
=\left[I+\left(e^{i\theta}-1\right)\Pi_\tau\right]|\psi\rangle.
```

The MPS backend constructs the same local operator from the two fusion-channel
projectors. In particular,

```math
U(0)=I,
```

while

```math
U(\pi)=\Pi_1-\Pi_\tau.
```

## RNG Draw Order

Gate choices, measurement outcomes, and random unitary angles use the same
`MeasureConfig.rng`. Within each event, the draw order is:

1. draw the gate type when ``0<p<1``;
2. for a measurement, draw its Born outcome;
3. for a random-angle unitary, draw its angle.

A fixed-angle unitary does not consume the third draw. Changing `p` or
`random_angles` can therefore shift the random-number stream seen by later
events, even when the initial seed is unchanged. A seed reproduces a run only
when the complete algorithm and configuration are unchanged.

Use a stored schedule, rather than only a seed, when deterministic replay is
required.

## Schedules and Deterministic Replay

A [`HybridGateSchedule`](@ref) stores three arrays of equal size:

- `measurement_mask[layer, column]`: `true` for a measurement and `false` for
  a unitary;
- `outcomes[layer, column]`: the measurement outcome, meaningful only where
  the mask is `true`;
- `unitary_angles[layer, column]`: the unitary angle, meaningful only where
  the mask is `false`.

The two-argument constructor

```julia
HybridGateSchedule(measurement_mask, outcomes)
```

creates a legacy fixed-angle schedule with `NaN` angle entries. During replay,
a finite stored angle takes precedence; a `NaN` angle at a unitary location
falls back to `config.θ`.

The following example samples and then exactly replays an exact-state
trajectory:

```julia
using FibonacciChain, LinearAlgebra, Random

N = 6
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
state = zeros(ComplexF64, length(anyon_basis(model)))
state[1] = 1

born_config = HybridConfig(
    τ = 1.0,
    t₂ = 4,
    mode = :Born,
    rng = MersenneTwister(7),
    enable_τ_eff = false,
    p = 0.5,
    random_angles = true,
)

sampled = bulk_evolution(model, state, born_config)

replay_config = HybridConfig(
    τ = 1.0,
    t₂ = 4,
    mode = :sample,
    enable_τ_eff = false,
    p = 0.5,  # Ignored during replay.
    random_angles = true,
)

replayed = bulk_evolution(
    model,
    state,
    replay_config,
    sampled.schedule,
)

@assert isapprox(replayed.state, sampled.state; atol=1e-13)
```

In `mode=:Born`, a schedule must not be supplied. In `mode=:sample`, a schedule
is required and must have the dimensions implied by `t₁`, `t₂`, and `N`.

## Evolution Results

Both the exact-state and MPS methods return a
[`HybridMeasurementOutcome`](@ref) containing:

- `state`: the final exact vector or MPS;
- `schedule`: the generated or replayed circuit realization;
- `free_energys`: one measurement surprisal per layer;
- `entanglement_entropys`: one half-chain entropy per period.

For exact-state evolution, `track_y_expectation=true` in the embedded
`MeasureConfig` also records `y_expectation_values` after every period. It uses
the same normalized charge-expectation helper as measurement-only evolution
and constructs a dense charge matrix. The vector is empty when tracking is
disabled. Hybrid MPS evolution supports the same option by contracting
`topological_charge_mpo` after each period, without constructing a dense charge
matrix. MPS charge tracking requires periodic boundaries; finite truncation
can affect the measured charge and should be checked for convergence.

For normalized evolution, the layer free energy is

```math
F_{\mathrm{layer}}
=-\sum_{m\text{ in layer}}\log q_m,
```

where the sum includes measurement events only. Unitary events do not
contribute. Exact evolution computes the entropy from the anyonic reduced
density matrix, while MPS evolution uses the central bond entropy.

The MPS method additionally uses `cutoff`, `mindim`, `maxdim`, and
`truncate_every_events` from the embedded `MeasureConfig`. Setting
`enforce_fibonacci_constraint=true` reprojects the MPS onto legal periodic
Fibonacci fusion paths after every layer.

## Topological-Sector Diagnostics

Two exact-state routines use the same stored hybrid schedules:

- `hybrid_lyapunov_spectrum` computes the finite-time Lyapunov spectrum
  of the unnormalized transfer-matrix product by periodic QR factorization;
- `hybrid_bayesian_evolution` evolves trivial- and ``\tau``-sector
  hypotheses under a common measurement record and returns likelihood,
  posterior, record-fidelity, conditional-entropy, mutual-information, and
  Bayes-error estimators.

For periodic Fibonacci chains, the topological charge operator used by these
routines has eigenvalues

```math
y_{\mathrm{trivial}}=\varphi,
\qquad
y_\tau=-\varphi^{-1},
\qquad
\varphi=\frac{1+\sqrt{5}}{2}.
```

Both topological-sector diagnostics require periodic boundary conditions. The
hybrid `bulk_evolution` methods do not impose this check, but the current
staggered schedule and regression tests are designed for even periodic chains.
