# Long Properly Colored Paths in the Union of Two Hamiltonian Paths

Independent research note (May 2026) on an extremal question about
two-edge-colored multigraphs whose color classes are Hamiltonian paths. The
paper gives a reduction lemma connecting the problem to a step in the
Norin–Steiner–Thomassé–Wollan approach to the Lovász vertex-transitive path
problem, proves a tight extremal value for a specific block-reversal family,
and refutes the natural half-bound conjecture with an explicit counterexample
at n = 18.

The full writeup is in [paper.pdf](paper.pdf).

## Background

For a permutation σ ∈ S_n, let H_σ be the two-edge-colored multigraph on
{0, 1, ..., n−1} whose blue edges form the natural-order Hamiltonian path and
whose red edges form the Hamiltonian path induced by the value order of σ. A
properly colored (PC) path is a simple path whose consecutive edges alternate
in color. Let ρ(σ) be the maximum number of vertices on such a path. The
parameter ρ_min(n) := min over σ of ρ(σ) controls a matching-edge invariant
appearing in work of Norin, Steiner, Thomassé, and Wollan (arXiv:2505.08634,
2025) on Lovász's vertex-transitive path conjecture: any linear lower bound
ρ_min(n) ≥ cn with c > 0 would push the best bound on that problem from
Ω(n^{9/14}) to Ω(n^{2/3}).

## Key results

- **Conversion Lemma.** A PC path on k vertices in H_σ yields a simple path
  in the associated graph G(n, σ) using exactly k matching edges. The
  construction is explicit, so any uniform lower bound on ρ transfers to the
  matching-edge parameter f(n).

- **Tight family.** The block-reversal permutation σ_m ∈ S_{4m+3} (composition
  (2, 3, 1, 3, 1, ..., 3, 1, 3, 2) with each block reversed) achieves
  ρ(σ_m) = 2m+3 = (n+3)/2. An explicit PC path attaining this value is given,
  and a matching upper bound follows from a local block calculus.

- **Refutation of the half-bound conjecture.** The conjecture
  ρ(σ) ≥ ⌈(n+3)/2⌉, supported by exhaustive enumeration through n ≤ 14, is
  false. At n = 18, the permutation
  σ* = (1, 0, 4, 3, 2, 5, 8, 7, 6, 11, 10, 9, 12, 15, 14, 13, 17, 16)
  has ρ = 10 < 11. It is obtained by deleting one interior singleton of σ_4.

- **Exhaustive computational verification.** ρ_min(n) is computed for every
  n ≤ 14. The n = 14 run covers all 8.72 × 10^{10} permutations of S_14 and
  finds no violation of the half-bound at that size.

- **Stacking and a conditional asymptotic bound.** A stackable variant of the
  σ* construction gives an infinite family with slope Δρ/Δn = 5/11 at every
  recorded k (verified for k ≤ 4). If the slope continues asymptotically, the
  true constant c_* := liminf_n ρ_min(n)/n satisfies c_* ≤ 5/11. A separate
  module search up to size 28 produces slope 10/28 at k = 1.

## File structure

```
.
├── paper.pdf                            full writeup
├── README.md                            this file
├── LICENSE                              MIT
└── src/
    ├── alt_verify_resumable.cpp         exhaustive S_n verifier (paper §3.4)
    ├── verify_n18_counterexample.cpp    independent check of ρ(σ*) = 10 (§4.1)
    ├── deletion_test.cpp                singleton-deletion experiment (§4.2)
    ├── module_search.cpp                module search behind the stacking family (§4.3)
    ├── block_deficit_census.cpp         block-by-block deficit statistics
    └── landscape_analysis.cpp           adjacent-transposition local-minima census (§6)
```

## Build and run

All sources are C++17 with OpenMP. Tested with `g++` 12+ on Linux and the
MinGW-w64 toolchain on Windows. Compile any source with:

```
g++ -O3 -march=native -std=c++17 -fopenmp src/<file>.cpp -o <file>
```

Sample invocations:

```
# Verify the half-bound conjecture exhaustively at n=10 (seconds)
./alt_verify_resumable 10

# Same at n=14 (about 5 days on 16 OpenMP threads; supports --resume)
./alt_verify_resumable 14 --resume

# Independent check of the n=18 counterexample (seconds)
./verify_n18_counterexample

# Adjacent-transposition local-minima census at n=11
./landscape_analysis 11

# Block-composition deficit statistics at n=12
./block_deficit_census 12
```

`alt_verify_resumable` writes a per-chunk log (`n<n>_chunks.log`) and can
restart from where it left off via `--resume`; running it with no arguments
prints the full flag list.

## Open questions

- Is c_* > 0? Equivalently, is there an absolute constant c > 0 such that
  ρ(σ) ≥ cn − O(1) for every σ?
- Does ρ(σ) ≥ ρ(τ_σ) hold for every σ, where τ_σ is the block quotient
  permutation? (Verified empirically through n ≤ 11.)
- Stacking conjecture: for M = (3, 1, 3, 3, 1), does the family σ^{(k)}
  satisfy ρ(σ^{(k)}) = 5k + 5 for every k ≥ 0? A positive answer gives
  c_* ≤ 5/11.

## Notes

The verification scripts in `src/` were written with assistance from Claude
Code.

## Contact

Alejandro Lizardi, alejlizardi05@gmail.com
