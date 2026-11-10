# Fortran precision refactor

Updated: 2026-11-10. Branch: `issue-24-arbitrary-precision`.
Status: baseline capture and shared precision proof in progress; kernel migration has not started.
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

### Native source survey

All four double/quad kernel pairs differ numerically only in `selected_real_kind(8)` versus `selected_real_kind(33)`; remaining differences are comments, spacing, and continuation layout. Both complex kernels currently set `knd1 == knd`. Preserve their separate declarations and conversion sites, including Bouwkamp branches. Keep the existing wrappers and ABI names.

Shared text does not imply identical precision policies: retain `kindd`/`kindq` branches, native accuracy targets, and runtime `precision`/`range` behavior. MPFR currently selects quad branches while deriving precision and tolerances from runtime bits.

### Translator replacement inventory

Line references below refer to `deps/generate_mpfr_fortran.cmake` at the planning commit.

| Current transformation | Maintained replacement |
| --- | --- |
| Parse statements, comments, continuations, arguments (10–106) | Compiler-preprocessed `.F90`; preserve attribution and OpenMP directives. No replacement parser. |
| Rewrite `real`/`dble`, `cmplx`/`dcmplx`, and kind arguments (108–160) | Explicit backend conversion names; preserve integer division and `knd1` conversions. |
| Rewrite floating and complex literals outside strings (162–238) | Native constants or direct decimal-string MPFR construction; explicit complex construction. |
| Remove `format`; replace I/O with labeled `continue` (246–257) | Conditional diagnostic I/O with preserved branch labels. |
| Replace real/complex declarations; remove scalar defaults (258–267) | Explicit declaration macros and runtime MPFR initialization; preserve native constant initialization. |
| Replace kinds and module names (268–271) | Target-specific kind/name definitions; preserve native branches and MPFR quad-branch selection. |
| Replace integer `minacc` assignments (276–278) | Maintained MPFR policy `max(8,int(mpfr_bits*0.3010299956639812d0)-8)` under the original conditions. |
| Replace cache equality and element sizes (279–284) | Exact MPFR equality; native epsilon-scaled equality; MPFR descriptor bytes + 88 + `8*ceil(bits/64)` storage accounting. |
| Infer and guard workspace (288–389) | Explicit per-kernel dimensions and overflow-safe accounting before `main`: two copies per real, four per complex, lazy arrays, 1 MiB overhead, and cache allowance; preserve status 4 rejection. |
| Release scalar qleg cache storage (399–406) | Deallocate MPFR scalar data during cleanup; assigning zero would retain allocations. |
| Inject imports and implicit MPFR typing (413–419) | Explicit imports and `implicit none`; declare implicit locals and bare `dimension` arrays. |
| Insert workspace call and complete cache cleanup (420–427) | Maintained guard and `mpfr_clear_caches`, preserving batch entry/exit cleanup and precision/rounding/status restoration. |
| Translate four quad kernels into generated files (435–440) | Compile each shared source per backend with isolated module directories; retain translation for unmigrated kernels. |

Prove literal syntax without C stringification or token pasting, for example `SW_LITERAL(0.1_knd,"0.1")`. Signed decimal text matters under directed rounding: `-mpreal("0.1")` differs from `mpreal("-0.1")`; the translator currently uses the former. Explain any correction against the baseline. Initialize runtime constants at active precision/rounding; derived-type copies retain their stored precision. Explicit initialization must not rely on unallocated MPFR values behaving as zero.

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

- Initial checkout verified on 2026-11-10: `9367bad97db388f467e3d986d23bc5777ddb0dbd`; this plan was the only untracked file and was committed as `1f71346`.
- Existing work: MPFR/MPC-backed Fortran, broad Julia `precision_bits` support, cache/workspace limits, lazy fallback allocation, and CI/Yggdrasil integration. The CMake translator remains required.
- Historical validation: 678 wave checks passed; the native comparison had 84 successes and 12 expected workspace rejections, with unchanged statuses/checksums after memory fixes. These are prior-session results, not checks rerun against this checkout.
- Historical measurements: matched native MPFR113 peaks fell from about 20.7 to 11.5 MiB; Julia angular/radial temporary allocations fell about 49%/14%. These are workload-specific baselines to recapture, not universal limits.
- Old `/tmp` libraries and comparison logs are absent. Rebuild locally; Yggdrasil is not a prerequisite. Cross-platform packaging/release still needs verification; direct recipe runs must replace the stale source pin.
- Remaining API limits: `accuracy` and coordinate-zero searches lack `precision_bits`; wave Jacobians require `order=1`; arbitrary-precision Jacobians require `h=nothing`.
- Refactor preparation: native source survey and translator replacement inventory recorded; no kernel migration or validation gate completed.
- Current phases: baseline capture and shared precision proof in progress.
- Decisions pending in phase 2: precise constant syntax and shared-source/header layout.
- Next action: establish the baseline from the actual checkout, then fix Debug build options and prove the shared precision definitions on a small example.

### Commit checkpoint: baseline and precision scaffolding

Paused at the user's request on 2026-11-10. HEAD is `1f71346`; the changes below are uncommitted. Suggestion: `build: prepare precision refactor validation and Debug builds`.

- Frozen baseline: `git archive 9367bad97db388f467e3d986d23bc5777ddb0dbd`, extracted to `.scratch/precision-refactor/baseline-source`; builds/logs are in sibling `baseline-build` and `baseline-logs` directories. Keep these separate from `debug-build` and `debug-logs`.
- Toolchain: Ubuntu 26.04 x86_64, GCC/GFortran 15.2.0, CMake 4.2.3, Julia 1.12.6, MPFR 4.2.2, MPC 1.3.1, GMP 6.3.0. Missing headers and GDB 17.1 were extracted from Ubuntu packages into `.scratch/precision-refactor/sysroot`; no system installation.
- Fresh baseline configured with `-DCMAKE_BUILD_TYPE=Release -DSWF_BUILD_TESTS=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON`. Local dependency options were `-DCMAKE_PREFIX_PATH="$PWD/.scratch/precision-refactor/sysroot/usr" -DMPFR_LIBRARY=/usr/lib/x86_64-linux-gnu/libmpfr.so.6 -DMPC_LIBRARY=/usr/lib/x86_64-linux-gnu/libmpc.so.3`.
- `ulimit -v 3145728; /usr/bin/time -v timeout 900 cmake --build .scratch/precision-refactor/baseline-build --parallel 1` passed: 158.87 s elapsed, 902304 KiB peak RSS. Original flags include both global and Release `-O3`.
- `ctest --test-dir .scratch/precision-refactor/baseline-build -j1 --output-on-failure` passed both tests in 1.76 s under the same cap; tests explicitly use four OpenMP threads.
- Native quad coverage passed independently: compile `test/native_cache.f90` with `gfortran -O3 -fopenmp -ffree-line-length-none`, include `baseline-build/mod_quad`, link `baseline-build/lib/libspheroidal_batch_quad.so` with its absolute rpath, then run with `OMP_NUM_THREADS=4 timeout 60`. New CMake target reproduces this coverage.
- Debug build passed with `-DCMAKE_BUILD_TYPE=Debug` and otherwise the same options, cap, timeout, and serialized build command. Compile commands show `-O0 -g3 -fcheck=all -fbacktrace` for Fortran and `-O0 -g3` for C; all three shared precision proof executables compiled. Debug tests have **not** run.
- Added opt-in required-backend path/ABI checks, enabled them in CI/Coverage, and prepared `scripts/precision_baseline.jl` for 48 raw cases per backend, each repeated twice, at MPFR113. These new Julia checks have **not** executed successfully yet.
- Julia setup initially failed because the local registry lacks the Manifest-pinned JLL and Preferences package. Exact pinned sources/artifact were fetched from Julia's package server into an isolated depot. The first comparison attempt failed while importing missing Preferences, before any numerical execution; Preferences is now present, and retry is pending. The Manifest was not changed.
- Native library SHA256 values: double `dae62294642ae332ce9dce0d2cc7c5eb35c3204f8bb99d0af1a60364934856c4`; quad `086121161646d78143e3eb22c65324b0956f45d3c44e5516ec1b8e51e6a84ecd`; MPFR `3b4b7b79e7f28863b793eb42d2bc97655d8ed7f4ea16f5b9b87cbcb1a6337f1a`.
- Prepared `deps/precision.inc` and `test/shared_precision.F90`: explicit types, dual native/string literal spelling, named conversions, native parameter versus runtime MPFR initialization, signed rounding, ownership, and precision changes. The long exact literal currently relies on GNU's unlimited line-length option; check portability before accepting phase 2.
- No kernel migration; translator and redundant kernels remain. No phase gate is marked complete. No debugger, Julia reference/performance, full suite, or packaging validation is claimed.

Resume with Debug CTest, then the raw baseline harness and Julia wave references. For Julia, set `JULIA_DEPOT_PATH="$PWD/.scratch/precision-refactor/depot:$HOME/.julia"`, all three `SPHEROIDALWAVES_LIBRARY_*` variables to the absolute baseline library paths, `SPHEROIDALWAVES_REQUIRE_BACKENDS=true`, and Julia/OpenBLAS/OpenMP thread counts to one. Use the 3 GiB cap, `/usr/bin/time -v`, `timeout 240`, and `julia --startup-file=no --heap-size-hint=1G --compiled-modules=existing --project=. scripts/precision_baseline.jl .scratch/precision-refactor/baseline-logs/native.tsv`. Complete debugger and direct-compiler proof checks before migration.
