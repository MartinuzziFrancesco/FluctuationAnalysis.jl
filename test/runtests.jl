using ReTestItems
using FluctuationAnalysis

# Each test file holds one or more standalone `@testitem`s that ReTestItems runs
# in isolated modules and, with `nworkers > 0`, in parallel. The worker count is
# taken from `RETESTITEMS_NWORKERS` (set in CI) and otherwise defaults to a small
# share of the local cores.
const NWORKERS = parse(
    Int, get(ENV, "RETESTITEMS_NWORKERS", string(min(Sys.CPU_THREADS, 4)))
)

runtests(FluctuationAnalysis; nworkers = NWORKERS)
