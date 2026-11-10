# Fortran precision refactor

Updated: 2026-11-10. Branch: `issue-24-arbitrary-precision`.
Status: planning complete; implementation has not started.
Planning reference: `9367bad97db388f467e3d986d23bc5777ddb0dbd`.
Related work: [Issue #24](https://github.com/brandynlucca/SpheroidalWaves.jl/issues/24).

## Objective and constraints

Replace the CMake translator with four shared Fortran kernels: real/complex prolate and real/complex oblate. Compile each for double, quad, and MPFR/MPC using explicit precision definitions and compiler preprocessing.

- Preserve numerical algorithms, native defaults, public interfaces, normalization, branches, result shapes, and thread safety.
- Retain three libraries; MPFR precision remains selectable at runtime without a native-precision bottleneck.
- Reuse `deps/mpfr_real.f90`, `mpfr_complex.f90`, `mpfr_batch.f90`, and the arithmetic C bindings. No Python or separate C++ numerical implementation.
- Keep numerical operations in Fortran; restrict macros to declarations, constants, and names. No replacement source parser.
- Use literal operation names and at most one physical line per new code comment. Preserve upstream attribution.
- Track remaining Issue #24 features separately. Julia runtime memory reduction is outside this refactor.

## Current implementation

`deps/*_swf.f90` and `deps/*_swf_quad.f90` contain native kernels. `deps/generate_mpfr_fortran.cmake` translates quad kernels into MPFR sources. Keep the generator until its replacement passes validation.

Integration: `CMakeLists.txt`, `scripts/build_from_source.jl`, `contrib/yggdrasil/build_tarballs.jl`, and CI/Coverage/Yggdrasil workflows. Preserve library products, entrypoints, and separate backend module directories.

## Ordered work

### 1. Capture a reproducible baseline

- [ ] Record starting commit/diff, toolchain, build flags, and library hashes. Keep baseline and migration builds separate.
- [ ] Inventory double/quad differences, including `knd1` and precision-dependent branches; preserve intentional behavior.
- [ ] Assign every generator transformation an explicit replacement: declarations, literals/conversions, tolerances, I/O/labels, cache equality/bytes/cleanup, and workspace checks.
- [ ] Require intended library paths and ABI symbols before every validation run; missing backends must fail.
- [ ] Capture native/Julia reference results and memory/time measurements with fixed inputs and commands. Recapture missing temporary baselines; historical counts are not current validation.

Gate: reproducible baseline and complete transformation inventory recorded.

### 2. Establish debugging and precision definitions

- [ ] Remove global `-O3`; preserve Release optimization and enable genuine Debug symbols, bounds checks, backtraces, and unoptimized stepping for Fortran/C.
- [ ] Prove a small shared `.F90` example for all three precisions. Record source/header layout and literal syntax here.
- [ ] Test decimal/complex constants beyond 113 bits under multiple rounding modes, conversions, integer division, and initialization. MPFR literals must avoid intermediate native rounding.
- [ ] Account for GNU traditional preprocessing; avoid assumed C stringification/token-pasting semantics. Separate runtime MPFR initialization from native constant expressions.
- [ ] Apply definitions per target; verify `#include` dependencies, module isolation, MPFR with quad disabled, and direct compiler builds without the translator.

Gate: arithmetic/ownership tests pass; debugger transcripts show maintained source lines, arguments, backtraces, stepping, and MPFR value inspection for each backend.

### 3. Migrate kernels individually

Complete each row across double, quad, and MPFR before moving to the next.

| Kernel | Status | Commit and validation evidence |
| --- | --- | --- |
| Real prolate | Pending | — |
| Real oblate | Pending | — |
| Complex prolate | Pending | — |
| Complex oblate | Pending | — |

For each kernel:

- Make declarations, constants, conversions, and initialization explicit; eliminate ambiguous implicit numeric declarations.
- Move tolerances, cache policies/cleanup, and overflow-safe workspace accounting into Fortran. Count lazy arrays and reject oversized requests before allocation.
- Preserve lazy fallback arrays, cache budgets, initialized unused dimensions, skipped-loop initialization, and oblate `maxn >= maxj`.
- Preserve wrapper/module/symbol names. Compile migrated sources for all three backends; retain existing build paths for remaining kernels.
- Compare raw native results and public APIs against the baseline and independent references. Exercise skipped/executed fallbacks, cold/reused/replaced caches, and mode switches.

Gate per row: numerical, cache/rounding, debugging, memory, and performance comparisons pass; explain differences before acceptance. Julia fallbacks cannot establish native correctness.

### 4. Validate the complete replacement

- [ ] Run native cache/mode checks against double and quad, plus MPFR arithmetic/ABI checks. Add missing quad coverage; `native_cache` currently links only double.
- [ ] Run all four Julia test groups: `api-and-precision`, `parameter-derivatives`, `coordinates-and-endpoints`, and `functions-and-coefficients`.
- [ ] Cover eigenvalues, coefficients, waves/derivatives, combined calls, sweeps, inverse calculations, and normalization quantities across supported inputs, kinds, endpoints, and static limits.
- [ ] Require improvement against independent nonzero-parameter references at 160/256 bits. Compare raw MPFR113 with native quad separately from public guarded `precision_bits=113` calls.
- [ ] Check ownership, cache replacement, concurrent calls, and precision/rounding restoration after success/failure. Test oversized-request rejection and subsequent recovery.
- [ ] Compare memory/time with the baseline, separating native RSS, precompilation, consumer RSS, warm allocations, and repeated-call retention.

Gate: all required checks pass using the intended backends and existing accuracy tolerances.

### 5. Finish build and distribution integration

- [ ] Update source/install lists, source-build verification, and workflows. Preserve double-only and MPFR-without-quad configurations.
- [ ] Build/audit the BinaryBuilder recipe for existing targets using the migration commit, not its stale source pin. Do not execute target-built generators during cross-compilation.
- [ ] Validate JLL products, symbols, and dependencies on supported CI platforms with Julia 1.10 and stable.
- [ ] After replacement paths pass, remove the translator, `.gitignore` exception, and redundant kernels. Check references and build from a clean checkout without generated files.
- [ ] Add concise build/debug instructions and verify ordinary Julia precompilation and backend overrides.

Gate: clean source and BinaryBuilder builds/tests pass without the custom translator.

## Validation and execution notes

Reuse both native tests and existing Julia suites/reference data: quad/oblate accuracy, angular conventions, radial Wronskians/Hankel kinds, endpoints, complex branch continuation, and threading. No reference-generation dependency.

Serialize heavy builds/tests with time/address-space limits; earlier runs used a 3 GiB cap and 1 GiB Julia heap hint. The hint is not a limit. Start single-threaded; run dedicated threading checks separately. Keep temporary harnesses/binaries/logs ignored; record portable summaries here. Establish fresh performance budgets rather than hard-coding historical RSS.

## Session handoff

Start by reading this plan and inspecting the branch, commits, and working tree. Preserve user changes. Update checkboxes/status only after gates pass.

Before ending, record commit/diff, completed phase/kernel, exact test commands/results, failures, and next action below. Do not depend on `/tmp` surviving. Parallelize independent reviews/preparation; keep kernel migration and heavy validation sequential.

- Checkout verified on 2026-11-10: `9367bad97db388f467e3d986d23bc5777ddb0dbd`; this plan is the only untracked file. Implementation/build/CI changes are committed.
- Existing work: MPFR/MPC-backed Fortran, broad Julia `precision_bits` support, cache/workspace limits, lazy fallback allocation, and CI/Yggdrasil integration. The CMake translator remains required.
- Historical validation: 678 wave checks passed; the native comparison had 84 successes and 12 expected workspace rejections, with unchanged statuses/checksums after memory fixes. These are prior-session results, not checks rerun against this checkout.
- Historical measurements: matched native MPFR113 peaks fell from about 20.7 to 11.5 MiB; Julia angular/radial temporary allocations fell about 49%/14%. These are workload-specific baselines to recapture, not universal limits.
- Old `/tmp` libraries and comparison logs are absent. Rebuild locally; Yggdrasil is not a prerequisite. Cross-platform packaging/release still needs verification; direct recipe runs must replace the stale source pin.
- Remaining API limits: `accuracy` and coordinate-zero searches lack `precision_bits`; wave Jacobians require `order=1`; arbitrary-precision Jacobians require `h=nothing`.
- Refactor work completed: implementation survey and this plan only; no kernel migration or new test run.
- Next phase: 1, baseline capture.
- Decisions pending in phase 2: precise constant syntax and shared-source/header layout.
- Next action: establish the baseline from the actual checkout, then fix Debug build options and prove the shared precision definitions on a small example.
