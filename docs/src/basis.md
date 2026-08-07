# Basis Functions

This section covers the fundamental basis generation and manipulation functions for anyon chains.

## Anyon Basis Generation

The `anyon_basis` function generates the computational basis states for various anyon models in different symmetry sector:

### Fibonacci Anyons
For Fibonacci anyons, where quantum dimension is $\phi = \frac{1+\sqrt{5}}{2}$, the basis consists of configurations that satisfy the fusion constraints. Each site can be in state `0` ($\tau$ particle) or `1` (vacuum), with the constraint that no two adjacent $1$ particles can fuse (forbidden configurations `11`are excluded). If you are familiar with PXP model, you will quickly get it.

### Ising Anyons  
For Ising anyons (Majorana fermions), where quantum dimension $d=\sqrt{2}$, the basis includes all possible $\sigma$ and 1 configurations, representing the two types of Ising anyons.


### Nomral Spin 
quantum dimension $d=2$.

### Symmetry
The symmetry you can choose nowadays is:

- Translational symmetry $T$
- Topological symmetry $Y$

Potential symmetry:
- Inversion symmetry $I$

### Usage Examples

```julia
using FibonacciChain

# Generate Fibonacci basis for 6 sites with PBC
N = 6
model_fibo = AnyonModel(FibonacciAnyon(), N; pbc=true)
basis_pbc = anyon_basis(model_fibo)

# Generate Ising anyon basis  
model_ising = AnyonModel(SpinHalf(), N; model_type=:Ising, pbc=true)
basis_ising = anyon_basis(model_ising)

# Generate basis in momentum sector k=0
basis_k0, rep_dict = anyon_basis(model_fibo, 0)
```

## Hamiltonian Construction

Now we support 

- Ising interaction: $H = -\sum_i Z_iZ_{i+1} - \sum_i X_i$
- Ferromagnetic Fibonacci anyon interaction: $H = -\sum_i \Pi_{i}^0$
- Antiferromagnetic Fibonacci anyon interaction: $H = -\sum_i \Pi_{i}^1$
  
## Topological Symmetries


The topological symmetry operations are crucial for understanding the anyonic nature of the system. They encode how the fusion outcomes transform under topological operations.

## Reduced Density Matrices

For a given model and subsystems, compute reduced density matrices:

- `anyon_rdm(model, subsystems, state)`
- Sector version: `anyon_rdm_sec(model, subsystems, state, k)`

And we also have `disjoint_rdm` for two parallel chains with two different anyon types.