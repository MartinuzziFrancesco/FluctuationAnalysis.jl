I’d stop expanding AD support and focus on making the package release-ready.

Recommended order:

1. Finish the `fm/ad` change

- Review and commit the current ForwardDiff work.
- Document the narrow guarantee: differentiation with respect to series values under fixed configuration.
- Avoid promising compatibility with every AD backend.

2. Address the deferred correctness edge cases

- Define behavior for zero-variance and constant series.
- Handle zero or negative fluctuation quantities consistently, especially DCCA signed covariance.
- Validate logarithm domains before log-log fitting.
- Check insufficient segments, invalid scale ranges, and rank-deficient polynomial fits.
- Add explicit tests for every defined failure convention.

The principal formulas and reduction identities have already been corrected; these numerical/domain cases are the remaining correctness priority.

3. Strengthen numerical reliability

- Test very short and very long series.
- Add nearly constant and badly scaled inputs.
- Exercise `Float32`, `Float64`, and `BigFloat`.
- Add property tests for:
  - `dcca(x, x) == dfa(x)`
  - MFDFA at `q = 2` reducing to DFA
  - MFDMA at `q = 2` reducing to DMA
  - invariance under additive offsets
  - expected scaling under multiplying the signal

4. Stabilize the public API

- Decide which result fields and accessor functions are public.
- Confirm consistent keyword names across DFA/MFDFA and DMA/MFDMA.
- Ensure internal `__` functions stay unexported and undocumented as public API.
- Review exception types and validation messages.
- Avoid API expansion until these contracts settle.

5. Improve practical documentation

- One opinionated tutorial for selecting scales and detrending order.
- Guidance on interpreting DFA, DCCA, DMA, and multifractal outputs.
- Examples showing common invalid analyses: too few scales, overfitting trends, zero fluctuations.
- A reproducible comparison against published/reference data.

6. Benchmark and profile

- Measure allocations and runtime by series length, scale count, and overlap mode.
- Optimize only demonstrated bottlenecks.
- Pay particular attention to repeated polynomial design matrices and segment allocations.

7. Prepare a release

- Add CI across supported Julia versions.
- Verify compatibility bounds and package metadata.
- Run Aqua, JET, Runic, tests, and docs in CI.
- Add changelog/release notes describing the corrected reduction identities and ForwardDiff support.
- Tag a release only after correctness edge cases have explicit contracts.

The immediate next task should be item 2: resolve and test the remaining zero/degenerate-domain correctness behavior. That has much more scientific value than adding additional AD backends.
