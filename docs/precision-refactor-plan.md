# Fortran precision refactor

Local handoff only: do not commit or push this file, and do not add it to `.gitignore`. The user makes numerical-work commits; pause with suggested commit text and the completed checkboxes from Issue #31.

Updated: 2026-11-10. Branch: `issue-24-arbitrary-precision`.
Status: phases 1–2 and all four kernel migrations validated; real prolate is committed as `f35c719`, real oblate as `48d67d3`, complex prolate as `cfddc84`, and complex oblate as `adf75d1`. Build integration is committed as `3be2f1a`. CI fixes and the mixed-geometry MPFR regression are uncommitted; platform and artifact validation remain open.
Planning reference: `9367bad97db388f467e3d986d23bc5777ddb0dbd`.
Related work: [Issue #24](https://github.com/brandynlucca/SpheroidalWaves.jl/issues/24).

## Objective and constraints

Replace the CMake translator with four shared Fortran kernels: real/complex prolate and real/complex oblate. Compile each for double, quad, and MPFR/MPC using explicit precision definitions and compiler preprocessing.

- Preserve numerical algorithms, native defaults, public interfaces, normalization, branches, result shapes, and thread safety.
- Retain three libraries; MPFR precision remains selectable at runtime without a native-precision bottleneck.
- Reuse `deps/mpfr_real.f90`, `mpfr_complex.f90`, `mpfr_batch.f90`, and the arithmetic C bindings. No Python or separate C++ numerical implementation.
- Keep numerical operations in Fortran; restrict macros to declarations, constants, and names. No replacement source parser.
- Use literal operation names and at most one physical line per new code comment. Preserve upstream attribution.
- Keep the main CMake build readable and compact: consolidate repeated backend setup and organize native test registration separately. Do this during final integration, after the individual kernel migrations; the user explicitly requested no CMake restructuring during the complex-prolate batch.
- Track remaining Issue #24 features separately. Julia runtime memory reduction is outside this refactor.

## Current implementation

The four `deps/*_kernel.F90` sources now compile directly for double, quad, and MPFR. The old `*_swf.f90`/`*_swf_quad.f90` files and `generate_mpfr_fortran.cmake` remain unused but retained until final integration and artifact validation pass.

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

The proof uses `deps/precision.inc` and `test/shared_precision.F90`, with `SW_REAL`, `SW_COMPLEX`, and `SW_LITERAL(0.1_sw_kind,"0.1")`. Named Fortran conversions avoid C stringification/token pasting. Signed decimal text matters under directed rounding: `-mpreal("0.1")` differs from `mpreal("-0.1")`; the translator currently uses the former. Explain any correction against the baseline. Initialize runtime constants at active precision/rounding; derived-type copies retain their stored precision. Explicit initialization must not rely on unallocated MPFR values behaving as zero.

## Ordered work

### 1. Capture a reproducible baseline

- [x] Record starting commit/diff, toolchain, build flags, and library hashes. Keep baseline and migration builds separate.
- [x] Inventory double/quad differences, including `knd1` and precision-dependent branches; preserve intentional behavior.
- [x] Assign every generator transformation an explicit replacement: declarations, literals/conversions, tolerances, I/O/labels, cache equality/bytes/cleanup, and workspace checks.
- [x] Require intended library paths and ABI symbols before every validation run; missing backends must fail.
- [x] Capture native/Julia reference results and memory/time measurements with fixed inputs and commands. Recapture missing temporary baselines; historical counts are not current validation.

Gate: reproducible baseline and complete transformation inventory recorded.

### 2. Establish debugging and precision definitions

- [x] Remove global `-O3`; preserve Release optimization and enable genuine Debug symbols, bounds checks, backtraces, and unoptimized stepping for Fortran/C.
- [x] Prove a small shared `.F90` example for all three precisions. Record source/header layout and literal syntax here.
- [x] Test decimal/complex constants beyond 113 bits under multiple rounding modes, conversions, integer division, and initialization. MPFR literals must avoid intermediate native rounding.
- [x] Account for GNU traditional preprocessing; avoid assumed C stringification/token-pasting semantics. Separate runtime MPFR initialization from native constant expressions.
- [x] Apply definitions per target; verify `#include` dependencies, module isolation, MPFR with quad disabled, and direct compiler builds without the translator.

Gate: arithmetic/ownership tests pass; debugger transcripts show maintained source lines, arguments, backtraces, stepping, and MPFR value inspection for each backend.

### 3. Migrate kernels individually

Complete each row across double, quad, and MPFR before moving to the next.

| Kernel | Status | Commit and validation evidence |
| --- | --- | --- |
| Real prolate | Complete | `f35c719`; shared source for all three backends, native tests, 678 wave checks, raw/fallback comparisons and debugger checks passed. |
| Real oblate | Complete | `48d67d3`; shared source for all three backends; native, wave, cache/threading, fallback, debugger and optional-build checks passed. |
| Complex prolate | Complete (`cfddc84`) | Shared source for all three backends; Debug/Release, 678 wave checks, cache/threading, fallback, debugger, ABI, performance and optional-build checks pass. |
| Complex oblate | Complete (`adf75d1`) | Shared source for all three backends; native, Julia, cache/threading, fallback, debugger, ABI, performance and optional-build checks pass. |

For each kernel:

- Make declarations, constants, conversions, and initialization explicit; eliminate ambiguous implicit numeric declarations.
- Move tolerances, cache policies/cleanup, and overflow-safe workspace accounting into Fortran. Count lazy arrays and reject oversized requests before allocation.
- Preserve lazy fallback arrays, cache budgets, initialized unused dimensions, skipped-loop initialization, and oblate `maxn >= maxj`.
- Preserve wrapper/module/symbol names. Compile migrated sources for all three backends; retain existing build paths for remaining kernels.
- Compare raw native results and public APIs against the baseline and independent references. Exercise skipped/executed fallbacks, cold/reused/replaced caches, and mode switches.

Gate per row: numerical, cache/rounding, debugging, memory, and performance comparisons pass; explain differences before acceptance. Julia fallbacks cannot establish native correctness.

### 4. Validate the complete replacement

- [x] Simplify `CMakeLists.txt` after all four migrations: remove transitional source lists and generation plumbing, consolidate repeated target setup, and separate native-test registration. Review the resulting main build flow for readability; moving duplicated boilerplate to another file alone is insufficient.
- [x] Run native cache/mode checks against double and quad, plus MPFR arithmetic/ABI checks, after the complete replacement. Quad coverage now exists.
- [ ] Run all four Julia test groups: `api-and-precision`, `parameter-derivatives`, `coordinates-and-endpoints`, and `functions-and-coefficients`.
- [ ] Cover eigenvalues, coefficients, waves/derivatives, combined calls, sweeps, inverse calculations, and normalization quantities across supported inputs, kinds, endpoints, and static limits.
- [ ] Require improvement against independent nonzero-parameter references at 160/256 bits. Compare raw MPFR113 with native quad separately from public guarded `precision_bits=113` calls.
- [ ] Check ownership, cache replacement, concurrent calls, and precision/rounding restoration after success/failure. Test oversized-request rejection and subsequent recovery.
- [ ] Compare memory/time with the baseline, separating native RSS, precompilation, consumer RSS, warm allocations, and repeated-call retention.

Gate: all required checks pass using the intended backends and existing accuracy tolerances.

### 5. Finish build and distribution integration

- [x] Update source/install lists, source-build verification, and workflows. Preserve double-only and MPFR-without-quad configurations.
- [ ] Build/audit the BinaryBuilder recipe for existing targets using the migration commit, not its stale source pin. Do not execute target-built generators during cross-compilation.
- [ ] Validate JLL products, symbols, and dependencies on supported CI platforms with Julia 1.10 and stable.
- [ ] After replacement paths pass, remove the translator, `.gitignore` exception, and redundant kernels. Check references and build from a clean checkout without generated files.
- [ ] Verify ordinary Julia precompilation and backend overrides. Contributor build/debug documentation is deferred at the user's request.

Gate: clean source and BinaryBuilder builds/tests pass without the custom translator.

## Validation and execution notes

Reuse both native tests and existing Julia suites/reference data: quad/oblate accuracy, angular conventions, radial Wronskians/Hankel kinds, endpoints, complex branch continuation, and threading. No reference-generation dependency.

Serialize heavy builds/tests with time/address-space limits; earlier runs used a 3 GiB cap and 1 GiB Julia heap hint. The hint is not a limit. Start single-threaded; run dedicated threading checks separately. Keep temporary harnesses/binaries/logs ignored; record portable summaries here. Establish fresh performance budgets rather than hard-coding historical RSS.

## Session handoff

HEAD is `adf75d1` (the user's complex-oblate refactor commit); complex prolate is `cfddc84`, real prolate `f35c719`, and real oblate `48d67d3`. Baseline/Debug scaffolding is committed as `4bae30e`; the oblate bounds fix and proof cleanup are `5982979`. The local plan remains untracked after its file-only GitHub deletion in `9ba3460`; `.gitignore` is unchanged. Keep the translator and redundant kernels until artifact validation. The user declined contributor documentation: `deps/README.md` was removed; do not recreate it unless requested. Quote Issue #31 checkbox text literally when reporting completion.

### Foundation evidence before kernel migration

The baseline source is `git archive 9367bad97db388f467e3d986d23bc5777ddb0dbd`, extracted to `.scratch/precision-refactor/baseline-source`. Separate sibling directories hold `baseline-build`, `baseline-logs`, `debug-build`, and `debug-logs`; future work must retain that separation.

Toolchain: Ubuntu 26.04 x86_64, GCC/GFortran 15.2.0, CMake 4.2.3, Julia 1.12.6, MPFR 4.2.2, MPC 1.3.1, GMP 6.3.0, GDB 17.1. Headers/GDB were extracted into the scratch `sysroot`; pinned Julia dependencies use an isolated `depot`, with no Manifest change.

| Validation | Result | Elapsed / peak RSS |
| --- | --- | --- |
| Frozen Release build and native tests | Build passed; CTest 2/2; independent quad cache check passed | Build 158.87 s / 902304 KiB; CTest 1.76 s |
| Frozen raw backend comparison | 48 double, 48 quad, 46 MPFR successes and 2 expected status-4 rejections; every case twice | 15.32 s / 385796 KiB |
| Julia wave suite against frozen libraries | 678/678 checks | 3:48.68 / 1004036 KiB |
| Julia wave suite against current Release libraries | 678/678 checks | Test runtime 3:49.7; see `release-logs/precision-waves.time` |
| Current Debug CTest | 6/6, including double/quad cache, MPFR arithmetic/ABI and all precision proofs | CTest 6.82 s / timed-process peak 89876 KiB |
| Direct compiler precision proof | Double, quad and MPFR passed without the translator | 1.70 s / 73856 KiB |
| Current clean Release build and CTest | Build passed; 6/6 tests | Build 152.75 s / 905128 KiB; CTest 7.60 s |
| Current raw backend comparison | Byte-identical TSV to frozen baseline, including statuses and accuracy estimates | 16.32 s / 385652 KiB |
| MPFR without quad / double only | Clean builds and CTest 4/4 and 2/2; double-only default is Release | Builds 24.23 s / 492028 KiB and 21.64 s / 427636 KiB |

Raw MPFR113 differs from native quad by at most `3.38e-31` relatively. Expected workspace rejections are complex-prolate fixtures 2/3, mode 3. `baseline-logs/native.tsv` SHA256 is `8e3296d1fd64ffc84355030194d20996633066464b139301605ca4446b99a863`.

Frozen library SHA256 values:

- Double: `dae62294642ae332ce9dce0d2cc7c5eb35c3204f8bb99d0af1a60364934856c4`
- Quad: `086121161646d78143e3eb22c65324b0956f45d3c44e5516ec1b8e51e6a84ecd`
- MPFR: `3b4b7b79e7f28863b793eb42d2bc97655d8ed7f4ea16f5b9b87cbcb1a6337f1a`

Debug compile commands confirm Fortran `-O0 -g3 -fcheck=all -fbacktrace` and C `-O0 -g3`. Bounds checks exposed the oblate singleton `match(2)` access; both oblate kernels now allocate an extra match slot and zero the array. The proof checks `1+1e-40` against an exact rational at 160/256/384 bits under all five rounding modes, signed literals, conversions, initialization and ownership; its longest line is 121 columns. Header dependencies and isolated module directories were verified.

GDB transcripts in `debug-logs/gdb-{double,quad,mpfr}.log` show maintained source lines, arguments, stepping and backtraces; MPFR inspection includes the 160-bit value as hexadecimal text. GDB 17.1 crashes on `finish` from a derived-type return; setting a temporary breakpoint on the caller line and continuing works (`proof-mpfr.gdb`). No Python was used.

### Reproduction and next action

Run serially from the repository root with `ulimit -v 3145728`. Configure baseline/Debug builds with `-DSWF_BUILD_TESTS=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON` and the appropriate `-DCMAKE_BUILD_TYPE`; local dependency options are `-DCMAKE_PREFIX_PATH="$PWD/.scratch/precision-refactor/sysroot/usr" -DMPFR_LIBRARY=/usr/lib/x86_64-linux-gnu/libmpfr.so.6 -DMPC_LIBRARY=/usr/lib/x86_64-linux-gnu/libmpc.so.3`. Build using `timeout 900 cmake --build <build-directory> --parallel 1`.

For Julia, source `.scratch/precision-refactor/environment.sh`: it selects the frozen libraries, requires backend paths/symbols, sets Julia/OpenBLAS/OpenMP threads to one, and selects the isolated depot. The equivalent depot setting on another checkout is `JULIA_DEPOT_PATH="$PWD/.scratch/precision-refactor/depot:$HOME/.julia"`; regenerate absolute library paths there. Native threading tests explicitly use four OpenMP threads. Timed commands used:

```sh
/usr/bin/time -v timeout 240 julia --startup-file=no --heap-size-hint=1G --compiled-modules=existing --project=. scripts/precision_baseline.jl .scratch/precision-refactor/baseline-logs/native.tsv
/usr/bin/time -v timeout 300 julia --startup-file=no --heap-size-hint=1G --compiled-modules=existing --project=. -e 'include("test/require_backends.jl"); include("test/precision_waves.jl")'
/usr/bin/time -v timeout 180 ctest --test-dir .scratch/precision-refactor/debug-build -j1 --output-on-failure
/usr/bin/time -v timeout 180 bash .scratch/precision-refactor/direct-proof.sh
```

The direct-proof script records complete compiler/link commands. Debugger commands are in `debug-logs/proof.gdb` and `proof-mpfr.gdb`; run them with the scratch `sysroot/usr/bin/gdb -batch -x <command-file> <proof-executable>`. Logs and `/usr/bin/time -v` reports live beside each build's results.

Next: pause for the user's complex-oblate commit, then simplify CMake and complete cross-kernel/build/distribution validation. Do not stage this local plan. Full four-group Julia validation and BinaryBuilder/JLL/platform gates remain open; local testing needs no release, and direct recipe builds must replace the stale source pin.

Earlier Issue #31 checkpoints: `4bae30e` completes Debug/optimization; `5982979` adds no checklist completion. Real prolate completes its consolidation and production preprocessor selection; real oblate completes its consolidation. See the complex-oblate checkpoint below for current guidance. Contributor documentation is deferred at the user's request. Do not edit the issue automatically.

During explicit-initialization migration, audit odd-tail paired paths involving `eigorder`/`eigt(li+2)` and terminal complex `nar` assignment; these are review concerns, not established failures. Remaining API limits: `accuracy` and coordinate-zero searches lack `precision_bits`; wave Jacobians require `order=1`, and arbitrary-precision Jacobians require `h=nothing`.

### Real-prolate checkpoint

Seven code/test files: `CMakeLists.txt`, `deps/precision.inc`, `deps/precision.F90`, `deps/prolate_kernel.F90`, `deps/generate_mpfr_fortran.cmake`, `test/prolate_cache.F90`, and `test/mpfr_prolate_threads.f90`. The shared kernel preserves the ABI, native policies, runtime MPFR precision, lazy arrays, and workspace accounting. It declares previously implicit numeric locals and arguments explicitly. Four `dnorm` intermediates intentionally retain native default-real precision; MPFR uses `mp_real`. Nested first-iteration comparison and initial endpoint `jmax` avoid two undefined reads.

CMake now translates only the three unmigrated geometries. Keep literal module names in nested `#ifdef` branches: CMake's Makefile scanner did not expand macro module names or handle `#elif defined` module selection reliably. GFortran 15 also crashed when a generic was imported under both original and alias names; explicit conversion wrappers avoid that bug. Kernel native conversions still expand to the intrinsic, avoiding wrapper overhead.

Evidence lives in `.scratch/precision-refactor/prolate-migration/`:

- Clean Release build: 139.04 s, 904828 KiB peak. Final Debug/Release native tests passed 10/10 in 9.76/8.85 s, peaks 104116/102960 KiB; see each `*-logs/final-ctest.log`.
- MPFR-without-quad and double-only clean builds passed; final tests including permanent fallback fixtures passed 7/7 and 3/3 in 4.15/0.06 s, peaks 104120/28032 KiB.
- Raw TSV is byte-identical to the frozen baseline, SHA256 `8e3296d1fd64ffc84355030194d20996633066464b139301605ca4446b99a863`: 48 double, 48 quad, 46 MPFR successes plus two expected status-4 rejections, each twice. Raw run: 16.42 s / 385896 KiB.
- Current-library wave checks: 678/678, test time 4m25.4s, wall 4:31.53, peak 1049956 KiB. `validate-release.sh` records the exact backend overrides and commands.
- Wrapper compatibility: `scaled_batch.f90`, `psms_batch_fortran.f90`, and `oblate_batch_fortran.f90` remain unchanged and linked. The raw baseline directly executes scaled angular/radial and native eigenvalue entrypoints. `wrapper-abi.jl` additionally calls double prolate/oblate angular wrappers in both normalization modes and radial wrappers for kinds 1–4 at multiple coordinates: all 12 batch calls succeed with finite outputs and byte-identical snapshots against the frozen library (`frozen/wrapper-abi.log`, `migrated/wrapper-abi.log`). MPFR continues through `mpfr_batch.f90`.
- `prolate_driver.F90`, `run-driver.sh`, and `native-gates.sh` reproduce direct cold/warm outputs for eigenvalues, radial/angular values and derivatives, coefficients, endpoint identities, and fallback cases. Release snapshots match frozen Release exactly; Debug snapshots match legacy Debug exactly. Existing Debug/Release differences in double fallback results remain unchanged.
- Natural `r2eta` hits: double `(c,m,l,x1)=(200,50,108,0.1)`; quad `(500,100,264,0.1)`; MPFR113 `(200,50,107,0.1)`. Tests request 80/300/80 degrees. MPFR uses a test-only 512 MiB workspace allowance and restores it; default remains 256 MiB. Release fallback time/RSS: double 0.10 s/29120 KiB, quad 1.20 s/43200 KiB, MPFR 2.60 s/103736 KiB, matching baseline closely.
- `debug-logs/gdb-{double,quad,mpfr}.log` records maintained source, arguments, recurrence stepping, backtraces and MPFR value formatting. `fixture-{double,quad,mpfr}.log` records actual fallback execution and an allocated saved table. Step past procedure entry before inspecting dynamic array bounds in GDB.
- Warm quad benchmark, 100 calls per case over 18 cases, three alternating pairs: frozen 9.835/9.189/9.572 s; migrated 9.158/9.195/10.653 s. No consistent slowdown. Initial 20-call MPFR totals: frozen 4.012 s, migrated 3.896 s; RSS about 24 MiB in both.

Suggested commit title: `Share the real-prolate kernel across double, quad and MPFR`.

### Real-oblate checkpoint

Five code/test files: `CMakeLists.txt`, `deps/precision.inc`, `deps/oblate_kernel.F90`, `test/oblate_cache.F90`, and `test/mpfr_oblate_threads.f90`. CMake translates only the two complex geometries. Native defaults, exact MPFR constants, runtime precision, all 25 workspace groups, cache budgets, both lazy table families and exported C interfaces are preserved. The six historical default-real temporaries in `r1eta`/`r2eta` remain native default-real; MPFR uses `mp_real`.

Initialization fixes skip angular work for eigenvalue-only calls, define omitted derivatives and accuracy fields, initialize zero-term/empty-tail recurrence state, and suppress an orphan diagnostic write. Cache fixes preserve MPFR coefficient ownership, reconstruct two lookahead entries after cache hits, and recompute normalization when coefficients are reused. Audit equivalent cache contracts in the earlier prolate and remaining complex kernels during final cross-kernel validation.

Evidence is under `.scratch/precision-refactor/oblate-migration/`, compared against the committed `f35c719` libraries in `prolate-migration/`:

- Clean Release build: 125.91 s, 902144 KiB peak. Final Debug/Release tests: 14/14 each, 23.78/20.42 s. Clean MPFR-without-quad: 10/10; double-only: 4/4. Heavy runs were serialized under a 3 GiB address-space cap.
- Julia wave checks: 678/678, wall 3:47.78, peak 1027584 KiB. Raw ABI cases: 48 double, 48 quad, 46 MPFR successes plus two expected status-4 rejections, each twice. Raw numerical values are unchanged; only one MPFR radial accuracy field changes from 0 to 31 across its four output rows. New TSV SHA256: `849a50447a3d0d5ac61dec13c8e151961a6d042cfda05150ab7a0a5fa62deb4f`.
- Expanded matrix: 90 cases/backend, 1949 rows/backend, both Debug and Release. Every numerical value matches the prior library. Changes are limited to formerly undefined value-only angular accuracy and origin radial accuracy. Explicit wrapper checks execute 12 prolate/oblate angular/radial batches, both normalizations and radial kinds 1–4, with identical old/new outputs; exported C symbols are unchanged.
- Natural R1 fallback: `(c,m,lnum,x)=(500,50,2,0.2)` enters `r1eta` at degree 50 in all three backends. GDB inspects the allocated saved table. Natural R2 GDB fixtures: double `(500,0,350,0.2)`; MPFR113 `(500,0,310,0.25)` enters at degree 281. Permanent second-kind fixtures request 310 degrees. The MPFR tests temporarily allow 512 MiB and restore caller settings. Quad's bounded natural R2 search did not select that method; two analytic direct recurrence cases execute and check `r2eta` for every backend, without claiming natural quad selection.
- Cold/warm physical snapshots: native double/quad are byte-identical to baseline. Old MPFR R1 fixture returns NaNs in all four radial channels for both degrees; repaired results are finite, cold/warm identical, and agree with quad to about `1.94e-33` relative error. The larger MPFR R2 fixture remains byte-identical to baseline and agrees with quad within `2.80e-16`, consistent with native quad's 15-digit target. `physical-fallbacks/validation.log` passes, explicitly distinguishing the old NaN failure from all strict comparisons.
- `debug-logs/gdb-*.log` records maintained source locations, stepping, arguments/backtraces and MPFR value formatting. Native poison builds with uninitialized-real markers and floating-point traps pass zero-term and empty-tail fixtures. Use whole exact-size derived-type output arrays in tests: an optimized GFortran test temporary caused an invalid free before that fixture was corrected.
- Alternating 50-call benchmarks, three pairs: quad paired medians within 2%; MPFR changes about +0.95%, +4.31%, +0.07% across three cases, with RSS below 17 MiB. A small cost in the middle MPFR case remains possible; the initial roughly 10% slowdown did not persist.

Suggested commit title: `Share the real-oblate kernel across double, quad and MPFR`.

### Complex-prolate checkpoint

Six code/test files: `CMakeLists.txt`, `deps/precision.inc`, `deps/complex_prolate_kernel.F90`, `test/shared_precision.F90`, `test/complex_prolate_cache.F90`, and `test/mpfr_complex_prolate_threads.f90`. Keep this local plan out of the commit. CMake now translates only complex oblate; the generator and redundant original kernels remain available. CMake readability cleanup is explicitly deferred to final integration at the user's request.

The shared source retains `knd1` declarations, exact MPFR complex construction, native accuracy branches, runtime precision, lazy tables, all 19 workspace groups and existing interfaces. It fixes Legendre coefficient ownership/lookahead/normalization, initializes omitted derivatives and endpoint Bessel state, and bounds the eigenvalue recurrence against its allocated coefficient arrays. Initialize saved Bessel term counts once per order, not per degree: later degrees use the preceding count to choose truncation.

Evidence under `.scratch/precision-refactor/complex-prolate-migration/`, baseline commit `48d67d3aceb228158b117425c6e01214dbb55874`:

- Fresh raw baseline and migrated Release TSV are byte-identical: SHA256 `849a50447a3d0d5ac61dec13c8e151961a6d042cfda05150ab7a0a5fa62deb4f`; 48 double, 48 quad, 46 MPFR successes plus two expected workspace rejections, each twice at 113 bits.
- Expanded Debug matrix: 80 combinations/backend, each cold and warm; all 5,425 comparable output rows match exactly. Old double rejects all 15 combinations of `(c,m,lnum,x1)=(2+25i,0,4,2)` with a coefficient-array bounds error; the new implementation completes them, with eigenvalues agreeing with quad to about `1.96e-16` scaled error. Other cases cover endpoints, tiny imaginary parameters, both angular normalizations and all output modes.
- Final Debug tests: 18/18, 100.95 s, peak 166240 KiB. Release: 18/18, 95.95 s, peak 164656 KiB. Clean Release build: 147.95 s, peak 901944 KiB. All heavy work is serialized under a 3 GiB address-space cap.
- Permanent fallback fixtures enter `r2eta` naturally: double `(200+0.25i,50,108,0.1)`; quad `(500+0.25i,100,264,0.1)`; MPFR113 `(200+0.25i,50,110,0.1)`. They request 80/300/80 degrees and check finite results and Wronskians. MPFR temporarily allows 512 MiB, restoring the original budget afterward. GDB also confirms integral fallback execution; quad exercises the Legendre fallback.
- GDB logs `debug-logs/gdb-*.log` record maintained source, arguments, stepping, backtraces and MPFR value inspection. `branches-*.log` and `extra-quad-2.log` record fallback execution. The four-thread MPFR regression checks 159/160/192/256 bits, five rounding modes, caller-context restoration, rejection isolation and recovery.
- Julia wave checks: 678/678, test time 3m58.6s, wall 4:04.23, peak 1029544 KiB. Clean MPFR-without-quad: 13/13 tests; double-only: 5/5. `release-gates.sh` and `remaining-gates.sh` record the commands and explicit backend paths.
- Double/quad fallback snapshots exactly match baseline, cold and warm; all migrated cold/warm outputs match. MPFR eigenvalues and R1 are unchanged, while R2 changes at degrees 110–129 following the cache repairs. All 20 changed degrees improve against a 256-bit reference: worst relative error drops from `7.17e-22` to `3.92e-27`. The reference run uses a test-only 1 GiB workspace budget because 512 MiB correctly rejects the 256-bit request; actual peak RSS is 190180 KiB, wall 6.91 s, under the same 3 GiB process cap. See `fallback-precision.log` and `reference_driver.F90`.
- Existing real wrappers execute 12 angular/radial batches with identical outputs; exported C symbols remain unchanged (24 double, 45 quad, 21 MPFR). The raw ABI snapshot also executes `scaled_batch` and both complex geometries.
- Three alternating benchmark pairs, 20 warm combined calls at `(c,m,lnum,x1)=(5+2i,1,3,2)`: median new/old ratios 0.993 double, 0.959 quad, 1.004 MPFR; no consistent slowdown. Warm RSS stays below 16 MiB. Release fallback RSS is 30592/48256/164516 KiB for double/quad/MPFR. `validation-summary.log` contains the comparisons.

Suggested commit title: `Share the complex-prolate kernel across double, quad and MPFR`.

At the complex-prolate checkpoint, Issue #31's consolidation checkbox was already checked; general migration checkboxes still awaited complex oblate. No commits, staging, pushes or issue edits were performed by the assistant.

### Complex-oblate checkpoint

Five code/test files: `CMakeLists.txt`, `deps/precision.inc`, `deps/complex_oblate_kernel.F90`, `test/complex_oblate_cache.F90`, and `test/mpfr_complex_oblate_threads.f90`. All four shared kernels now build directly for double, quad, and MPFR, without CMake translation. Retain the translator and redundant kernels until final integration. The broader CMake cleanup remains deferred to the next batch; no contributor documentation is requested.

The migration preserves native defaults, separate `knd1` declarations/conversions, exact MPFR constants, runtime precision, 25 workspace groups, lazy radial tables, and existing wrapper entrypoints. Cache repairs preserve coefficient ownership, reconstruct the two lookahead entries, and recompute normalization with reused coefficients. Initialization fixes define omitted angular derivatives, zero-tail R1 convergence state, and terminal radial accuracy for odd paired batches.

Evidence is in `.scratch/precision-refactor/complex-oblate-migration/`, against the user's `cfddc84` libraries in `complex-prolate-migration/`:

- Fresh raw baseline SHA256: `849a50447a3d0d5ac61dec13c8e151961a6d042cfda05150ab7a0a5fa62deb4f`; 48 double, 48 quad, 46 MPFR successes and two expected workspace rejections, each twice.
- Full Debug CTest passed 22/22 in 277.35 s, peak 166452 KiB. After the terminal-accuracy fix, all three affected cache tests passed again in 74.62 s. Final Release CTest passed 22/22 in 246.47 s, peak 164488 KiB. The initial fresh Release build reached a 300 s timeout without compiler errors; its bounded continuation passed with a 600 s allowance. All heavy work remains serialized under 3 GiB.
- Expanded matrix: 90 mode combinations and 2230 recorded rows per backend. Double/quad rows are unchanged. MPFR repairs 80 formerly nonfinite rows at `c=2+25i`, including 60 radial-value rows and 20 Wronskian diagnostics. Other rows are unchanged.
- Against MPFR256, the repaired MPFR113 high-imaginary values have worst relative error `1.36e-28`, within their 26-digit target; native double/quad errors are `5.33e-13`/`1.67e-29`. The MPFR256 Wronskian error is `5.77e-26`. At 113 bits this Wronskian loses over 44 digits to cancellation, so the native tests check finite radial results and reuse there instead.
- Odd paired fixture `c=10+0.25i, m=0, lnum=3, x=2`: frozen Debug leaves the last accuracy sentinel untouched in all three backends. Repaired values are 13 digits in double and 31 in quad/MPFR; every numerical value is unchanged. The exploratory raw `c=20+0.25i` odd batch falls below the kernel's documented even-batch threshold; its last double Wronskian cannot pass the tolerance. Its numerical outputs match the frozen implementation, and native double/quad wrappers already pad this case to an even batch. No numerical tolerance was weakened.
- GDB records maintained source, arguments, stepping, backtraces and MPFR value formatting for all three backends. R1 fallback is selected naturally at `c=500+0.25i, m=l=50, x=0.2`; MPFR R2 fallback at `c=500+0.25i, m=0, l=282, x=0.25`. Two analytic zero/one-term fixtures execute both radial recurrences for every backend. Complex fallback tests temporarily allow 1 GiB workspace, then restore caller settings; production defaults are unchanged.
- Julia waves passed 678/678 (8m45.2s test time, 8:58.84 wall, 1051276 KiB peak). Final raw output is byte-identical to the fresh baseline; all 12 real-wrapper calls and exported C symbols (24 double, 45 quad, 21 MPFR) are unchanged. Native fallback snapshots match the baseline, and all repaired MPFR snapshots are finite and cold/warm identical; four previously nonfinite fallback rows are repaired.
- Clean MPFR-without-quad passed 16/16 tests; clean double-only passed 6/6 and retains the default Release configuration.
- The repaired `c=500+0.25i` MPFR113 fallback agrees with native quad to `1.03e-29` relative error. Separate-process timings remained noisy even with CPU pinning. An interleaved C ABI harness loads distinct old/new libraries with local deep binding, verifies identical outputs, and alternates calls on CPU 0: double ratios 0.999/0.970/0.999 and MPFR 1.001/1.018/0.971 for eigen/angular/radial modes. Quad ratios were 1.026/1.139/0.900; a focused angular repeat with reversed library load order gave 0.988. No consistent regression was established; retain the individual measurements rather than claiming exact timing equality.

Suggested commit title: `Share the complex-oblate kernel across double, quad and MPFR`.

Issue #31: complex-oblate consolidation is already checked. This checkpoint completes seven implementation items: explicit declarations/initialization; precision-aware constants/constructors/conversions; maintained accuracy/truncation policies; cache policies/cleanup; overflow-safe workspace checks; conditional I/O; and CMake compilation of all shared kernels with isolated module directories and unchanged interfaces. Leave source-build/CI/Yggdrasil integration, translator removal, and CMake simplification open. Contributor documentation remains deferred. No staging, commits, pushes, issue edits, or `.gitignore` changes were made.

### Build-integration checkpoint

Uncommitted batch after `adf75d1`: root CMake is 100 lines (previously 265), with 33 lines in `test/CMakeLists.txt`. Backend setup and geometry test registration are shared. Library names, module isolation, wrapper sources, optimization modes, and install products are preserved; native tests run serially with explicit OpenMP thread counts and 180-second limits.

The optional Julia source-build script verifies maintained kernels and headers, supports backend/configuration/build-directory selection, honors `FC` and `CMAKE_GENERATOR`, locates actual libraries across lib/bin/configuration directories, verifies loading and ABI symbols, and returns nonzero on failure. Ordinary package installation still uses JLL artifacts without this script. CI uses the script and adds Debug/double-only/MPFR-without-quad jobs. Coverage uses the same source-build path. The direct recipe now pins `adf75d1e31886132cb9a37bb27feed6e6a590d33` and builds serially; the workflow still substitutes the selected revision and requires all three backends before quad/regression/wave artifact tests. No contributor documentation was added.

Evidence: `.scratch/precision-refactor/build-integration/`; commands are in `validate-release.sh`, `prepare-source.sh`, `validate-configurations.sh`, and `check-integration.sh`. All heavy commands run serially under a 3 GiB virtual-memory limit, with a 1 GiB Julia heap hint.

- Fresh Release source build, 22/22 CTests (249.95 s), and installation passed; combined wall time 9:58.91, peak 901836 KiB.
- Debug build from a copy without legacy kernels or translator passed 22/22 CTests (315.44 s) and installation; combined wall 6:42.91, peak 492696 KiB. Fortran/C flags retain `-O0` and debug information; Fortran bounds/backtrace flags remain active. DWARF names all four maintained kernels for all backends.
- Removing quad wrappers from that copy still permits MPFR-without-quad: 16/16 CTests (151.77 s) and installation passed; combined wall 7:06.25, peak 901888 KiB. Removing MPFR Fortran/C files also permits double-only: 6/6 CTests (1.08 s), installation, wall 1:07.19, peak 427560 KiB; no C compiler is configured.
- Fresh baseline and installed-library raw outputs are byte-identical: SHA256 `849a50447a3d0d5ac61dec13c8e151961a6d042cfda05150ab7a0a5fa62deb4f`, 142 successes and two expected workspace rejections, each twice. All exported symbols match (268 double, 308 quad, 631 MPFR, including Fortran internals), and installed dependencies resolve.
- The final source script also passes with native tests disabled. Missing maintained sources, missing `mpfr_arithmetic.h`, and invalid `FC` fail before configuration. Workflow YAML, embedded shell/Julia, and recipe syntax pass checks; the workflow's revision/version substitution was executed on a disposable copy.
- Installed-library Julia suites passed 1836/1836 checks: 820 quad, 338 numerical regressions, and 678 arbitrary-precision waves. Combined wall time 13:38.45, peak 1001368 KiB; wave test time 9m35.3s. Required backend paths and symbols were checked before execution.

BinaryBuilder target builds, generated-JLL platform runs, and the broader four-group Julia validation remain open. Retain legacy kernels, translator, and its existing `.gitignore` exception. No staging, commits, pushes, issue edits, or `.gitignore` edits were made.

Suggested commit title: `Integrate shared Fortran sources into builds, CI, and Yggdrasil`.

Eligible Issue #31 text, quoted literally:

- Update `scripts/build_from_source.jl`, CI workflows, and the Yggdrasil recipe for the shared sources.
- The `CMakeLists.txt` file has a fair bit of repeated setups and other plumbing components that makes it difficult to debug. This should be simplified.

## CI failures and mixed-geometry audit — 2026-11-10

User steered this batch to failing Coverage/Windows CI and platform timing. PR #32, HEAD `3be2f1af045740fab9e270055566e38959319313`; CI run `42927521044`, Coverage run `42927521257`. Evidence is under `.scratch/precision-refactor/ci-fixes/`; raw job logs include Ubuntu, Windows, both Macs, and Coverage.

- Coverage's Julia 1.10.13 failure reproduced exactly: lazy recurrence factors grew the retained cache to 1,105,728 bytes against a 1,105,312-byte reservation. Production accounting now measures the existing precision value and includes tuple/array storage. Added checks at 64/160/512 bits and 16/256/2048 factors. Julia 1.10.13 and 1.12.6 each pass 1901 arbitrary-precision checks against the installed shared-source backends; wall times 35.16/52.09 s, peaks 714312/1062292 KiB. The downloaded Julia 1.10.13 archive matched its official SHA256; dependencies were resolved in a scratch environment without changing the root project/manifest.
- Windows passed all double/quad tests but all ten MPFR tests crashed with `0xc0000374`. Its compiler came from `C:\mingw64`, while arithmetic development libraries came from a separate MSYS2 installation. CI now selects absolute UCRT64 C/Fortran compiler paths and matching math libraries, uses Ninja, and places eight runtime DLLs beside the native executables. Those same DLLs travel with the backend artifact, removing compiler/MSYS2 setup from Windows Julia consumers. YAML/Bash parsing and a mocked configuration with spaces in its prefix passed. This is not yet validated on Windows; the exact crashing DLL/root cause was not observable in the old job's artifacts.
- Intel macOS subsequently failed `mpfr_prolate_threads` with two mismatches; its other 21 native tests passed. The test now prints field values, accuracy, status, precision, rounding, and context diagnostics. CI retains CMakeCache/build/CTest logs on failure. Linux passed 100 consecutive prolate runs with `MALLOC_PERTURB_=165` (8.84 s), plus Debug. The macOS mismatch remains unresolved pending a run with diagnostics; no numerical tolerance was relaxed.
- CI compilation times: Ubuntu 1:52.51, Windows 3:12.46, ARM Mac 3:36.93, Intel Mac 5:32.45. Native tests: Ubuntu 113.90 s, ARM Mac 155.13 s, Intel Mac 390.21 s. Windows's 113.35 s includes immediate MPFR crashes and is not a successful-suite timing. Intel setup was about 30 s, while its build/test step took 12:23. No timeout increases made.
- New `test/mpfr_batch_context.f90` covers all four geometries, eigenvalue/angular/radial modes, 113/160/192/256 bits, all rounding modes, two coordinates, four OpenMP workers, seeded cache cleanup across modules, caller context restoration, and eleven failure/recovery cases. Release passed in 40.64 s; Debug also passed. The existing prolate test passed in both configurations. Builds/tests remained serialized under 3 GiB virtual memory, with Julia's 1 GiB heap hint. Broader retained-memory/performance auditing is still open.

No new Issue #31 checkbox is completed by this batch. The remaining unchecked entries are generator removal and documentation; the user explicitly deferred contributor documentation. Retain legacy sources/translator pending artifact gates. Do not stage/commit this local plan or change `.gitignore`.

Suggested commit title: `Fix cache accounting and Windows CI setup; expand MPFR regression diagnostics`.
