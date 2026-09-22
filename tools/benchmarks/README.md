# How to drive the other codes

One page per code MEQ is compared against, written so that picking one up again
after six months costs an hour rather than an afternoon. Each says the command
that works, the traps that cost real time, and — the part that matters most —
**what that code is not doing that MEQ is**.

| | |
|---|---|
| [`HOW_TO_DRIVE_FREEGS4E.md`](HOW_TO_DRIVE_FREEGS4E.md) | free boundary, Python, von Hagenow + Picard. The source of every `ref-n*/` case, and the one whose `psi_bndry` is a grid-node artefact |
| [`HOW_TO_DRIVE_DESC.md`](HOW_TO_DRIVE_DESC.md) | fixed and free boundary, inverse spectral, JAX. Works on the fixed ladder; its free-boundary solve mostly does not converge |
| [`HOW_TO_DRIVE_NICE.md`](HOW_TO_DRIVE_NICE.md) | free boundary, C++ P1, direct Newton. The fastest code here |
| [`HOW_TO_DRIVE_CHEASE.md`](HOW_TO_DRIVE_CHEASE.md) | fixed boundary, inverse, Fortran. The cheapest conversion and the closest agreement with MEQ |
| [`HOW_TO_DRIVE_TSC.md`](HOW_TO_DRIVE_TSC.md) | free boundary, vintage Fortran, card input. Carries a sizing bug that fires on any machine with no passive structure |

The harnesses themselves stay where they are — `tools/freegs4e-benchmark/`,
`tools/desc-benchmark/`, `tools/chease-benchmark/`, and `meq-race/` inside the
NICE and TSC trees. These are notes about driving them, not a fifth harness.

## The one lesson that is common to all five

**ASK WHAT INPUT THE TWO ARMS SHARE, AND WHETHER THAT SHARED INPUT COULD PRODUCE
THE AGREEMENT BY ITSELF.** Three separate times in this campaign a comparison
agreed because one side had been handed the other's answer, or because two paths
believed distinct called the same code:

* MEQ reproduces `freegs4e`'s `psi_bndry` on the limited case because MEQ is
  **given** that number — it is the value at the one grid node that set it;
* DESC's `--boundary reference` and `--blend 1.0` agree to every digit because
  they are one code path reached two ways;
* a DESC free-boundary run seeded from the reference reported a tiny error
  because the optimiser **had not moved**, and the harness never read
  `result["success"]`.

The cheap general defence is a refinement control: change the shared input and
see whether the agreement follows it. A real agreement is insensitive; a
manufactured one tracks. On the limited case `freegs4e`'s boundary node moves
0.6 m between the 129² and 513² grids, and MEQ's agreement follows it — which
is the test firing.

M-162 and M-164 in `MEASUREMENTS.md` are where this is recorded with numbers.
