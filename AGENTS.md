# AGENTS.md

Rules for AI coding agents in **FluctuationAnalysis.jl** — a type-generic Julia package for fluctuation/scaling analysis of time series (`dfa`, `mfdfa`, `dcca`, `dma`, `mfdma`, `hurst`). Priorities, in order: **correctness, clean API, type stability, extensibility, docs.**

## Commands

Always pass `--startup-file=no` (the user's `startup.jl` breaks precompilation on Julia 1.12).

```bash
# Test (ReTestItems, parallel)
julia --startup-file=no --project=test -e 'using Pkg; Pkg.test("FluctuationAnalysis")'
# Format (Runic — CI enforces it; run before committing)
julia --startup-file=no --project=@runic -m Runic --inplace src/ test/ docs/make.jl
# Format check only
julia --startup-file=no --project=@runic -e 'using Runic; exit(Runic.main(["--check", "--diff", "."]))'
# Build docs
julia --startup-file=no --project=docs docs/make.jl
```

## Workflow rules

- **Do not add new analysis methods without asking the user first.** When approved, verify the math against the cited primary source before implementing.
- **Test each function as you add it** — never batch testing to the end.
- Run the test suite **and** the Runic check before considering a change done.
- Snapshot a working state before big changes; commit each logical step.
- Reuse existing primitives (`integrated_profile`, `logarithmic_scales`,
  `__segment_views`, `__segment_variance(s)`, `__segment_covariance(s)`,
  `__fluctuation_curve`, `__q_order_fluctuation`, `loglog_fit`,
  `__scale_selection`) instead of duplicating logic. One concern per file.
- Extend via the existing seams: subtype `AbstractDetrender`, `AbstractHurstEstimator`, or `AbstractFluctuationResult`.

## Code rules

- **Follow the [SciML Style Guide](https://github.com/SciML/SciMLStyle).** Format with **Runic only** (never Blue/JuliaFormatter).
- **No single-letter or symbol variable names** — use `series`, `profile`, `scale`, `exponent`, never `x`, `s`, `α`.
- **Package-private helpers start with `__`** — exported API keeps ordinary
  `snake_case`; every unexported implementation function uses a double-underscore
  prefix (for example, `__segment_views`).
- **No unnecessary comments.** Let code and docstrings speak.
- **Each function does a single thing.** Keep functions small and composable.
- **Type-generic, never hardcode `Float64`.** Preallocate with
  `eltype`/`float`/`similar`; results preserve the input scalar type, including
  ForwardDiff dual numbers (Float32 → Float32, BigFloat → BigFloat).
- **Normalize inputs at exported boundaries** to a 1-based plain array; accept `Array`, views, `OffsetArray`, ranges, and `Int`/`Float32`/`Float64`/`BigFloat`/`Rational` eltypes.
- **Use `Int.(collect(scales))`, never `collect(Int, scales)`** — only 1-arg `collect` resets offset axes to `Base.OneTo`; the 2-arg form, comprehensions, and generators preserve them.
- **Keep per-segment statistics on `mean` (1/s).** This preserves the exact reduction identities (`dcca(x, x) == dfa(x)`, MFDFA `h(2)` → DFA, DMA = MFDMA at `q = 2`) — verify they still hold after touching the fluctuation machinery.
- **Structs use ConcreteStructs `@concrete`** (untyped fields). Exception: those with validating inner constructors (`PolynomialDetrender`, `MovingAverage`) stay plain.
- **Explicit imports** (`using Pkg: names` / `import Pkg`); blank line between `using` and `import` groups.
- **Every dependency** (incl. test deps) has a `[compat]` bound — no `^` or `>=`.

## Test rules

- **ReTestItems.** Put each concern in `test/<concern>_tests.jl` (the `_tests.jl` suffix is required) as standalone `@testitem "name" begin ... end` blocks; each declares its own imports and runs isolated/parallel.
- Split expensive statistical-convergence tests into their own testitems so they parallelize.
- When touching the public boundary, add coverage to `input_types_tests.jl` (every container/eltype variant must match the `Vector{Float64}` baseline).
- Keep `ad_tests.jl` passing: exported analyses support ForwardDiff derivatives
  with respect to series values while discrete configuration remains fixed.
- Aqua + JET live in `quality_tests.jl`.

## Doc rules

- Docstrings follow the **SciML standard**: signature line with `-> ReturnType`, summary, then `# Arguments` / `# Keywords` / `# Returns` / `# Throws` (`# Fields` for types). Document the abstract interface, not every subtype.
- For math use `@doc doc"""..."""` before the definition (`doc"..."` is raw — write `\tau`, `\alpha` directly); inline math uses Documenter double-backticks.
- Docs site (`docs/make.jl`, `checkdocs=:exports`): Home → Tutorials → API reference; tutorials precede the reference, all tutorial code in runnable `@example` blocks.
