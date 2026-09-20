using FibonacciChain
using ITensors
using ITensorMPS
using Random

"""Generate a coherent-state MPS trajectory using hybrid bulk_evolution."""
function samples_generate_mps(L::Int, p::Float64, periods::Int, seed::Int;
    stride::Int = 1, save_schedule::Bool = false,
    cutoff::Float64 = 1e-12, mindim::Int = 1, maxdim::Int = 256,
    truncate_every_events::Int = 1, enforce_fibonacci_constraint::Bool = false,
)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and >= 4"))
    periods >= 1 && stride >= 1 && seed >= 0 ||
        throw(ArgumentError("periods/stride must be positive; seed nonnegative"))
    isfinite(cutoff) && cutoff >= 0 && 1 <= mindim <= maxdim &&
        truncate_every_events >= 1 || throw(ArgumentError("Invalid MPS truncation settings"))
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    state, sites = initial_mps(L)
    Y = topological_charge_mpo(sites; pbc = true)
    initial_y = real(inner(prime(state), Y, state)) / real(inner(state, state))
    phi = (1 + sqrt(5.0)) / 2
    initial_weight = (initial_y + inv(phi)) / (phi + inv(phi))
    0 < initial_weight < 1 || error("Initial state must have both sector components")
    config = HybridConfig(
        τ = Inf,
        t₂ = periods,
        mode = :Born,
        rng = MersenneTwister(seed),
        enable_τ_eff = false,
        track_y_expectation = true,
        p = p,
        random_angles = true,
        cutoff = cutoff,
        mindim = mindim,
        maxdim = maxdim,
        truncate_every_events = truncate_every_events,
        enforce_fibonacci_constraint = enforce_fibonacci_constraint,
    )
    outcome = bulk_evolution(model, sites, state, config)
    times = sort!(unique!([0; collect(stride:stride:periods); periods]))
    entropy = [ee_mps(state, L ÷ 2); Float64.(outcome.entanglement_entropys)][times .+ 1]
    y_expectation = [initial_y; Float64.(outcome.y_expectation_values)][times .+ 1]
    measurement_count = count(outcome.schedule.measurement_mask)
    schedule = save_schedule ? outcome.schedule : nothing
    final_bond_dimension = maxlinkdim(outcome.state)
    return (; L, p, periods, seed, times, entropy, y_expectation,
        initial_weight, measurement_count, schedule, final_bond_dimension)
end

function process_task_mps(task)
    L, p, periods, seed, stride, save_schedule, settings = task
    return samples_generate_mps(L, p, periods, seed; stride, save_schedule, settings...)
end
