# Backend Overrides

This page describes how backend library paths are selected and how users can override them.

## Resolution Order

At module initialization, backend libraries are resolved in this order:

1. Package artifacts from `Artifacts.toml` (default user path)
2. Local libraries in the package's `build/bin`, `build/lib`, or `build` directory (developer fallback)
3. Environment variables (explicit override)

At runtime, explicit API calls always take precedence over automatic initialization:

- `set_backend_library!("/path/to/lib"; precision=:double)`
- `set_backend_library!("/path/to/lib"; precision=:quad)`

## Environment Variables

Set one or both variables before starting Julia:

- `SPHEROIDALWAVES_LIBRARY_DOUBLE` for double precision backend
- `SPHEROIDALWAVES_LIBRARY_QUAD` for quad precision backend

Windows PowerShell example:

```powershell
$env:SPHEROIDALWAVES_LIBRARY_DOUBLE = "C:\path\to\spheroidal_batch_double.dll"
$env:SPHEROIDALWAVES_LIBRARY_QUAD = "C:\path\to\spheroidal_batch_quad.dll"
julia
```

Linux/macOS shell example:

```bash
export SPHEROIDALWAVES_LIBRARY_DOUBLE="/path/to/libspheroidal_batch_double.so"
export SPHEROIDALWAVES_LIBRARY_QUAD="/path/to/libspheroidal_batch_quad.so"
julia
```

## Programmatic Override

Use explicit runtime configuration when paths are known inside application code:

```julia
using SpheroidalWaves

set_backend_library!("/path/to/libspheroidal_batch_double.so"; precision=:double)
set_backend_library!("/path/to/libspheroidal_batch_quad.so"; precision=:quad)
```

## Inspect Active Paths

```julia
backend_library(precision=:double)
backend_library(precision=:quad)
```

## Notes

- Released artifact entries provide out-of-box backend loading on supported platforms.
- If artifacts are missing for a platform, the package build step can compile local backends which are discovered on the next package load.
- When both precision artifacts are installed, the build step exits before checking for CMake or a Fortran compiler.
- Local builds are discovered directly; no generated runtime configuration file is needed.
- `Artifacts.toml` is the standard Julia package artifact manifest, not executable runtime configuration.
- Shared Fortran memoization is retained but made thread-private, so each Julia thread has independent caches.
- A small Julia lock protects only backend path and dynamic-library handle management; it is not held during numerical calls.
- If a configured path does not exist, it is ignored and a warning is emitted.
- If no backend is configured for a requested precision, calls fail with a clear error message.

