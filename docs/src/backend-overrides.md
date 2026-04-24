# Backend Overrides

This page describes how backend library paths are selected and how users can override them.

## Resolution Order

At module initialization, backend libraries are resolved in this order:

1. Environment variables (preferred override path)
2. Local generated config file `deps/library_config.jl` (fallback)

At runtime, explicit API calls always take precedence over automatic initialization:

- `set_backend_library!("/path/to/lib"; precision=:double)`
- `set_backend_library!("/path/to/lib"; precision=:quad)`

## Environment Variables

Set one or both variables before starting Julia:

- `SPHEROIDALWAVEFUNCTIONS_LIBRARY_R8` for double precision backend
- `SPHEROIDALWAVEFUNCTIONS_LIBRARY_R16` for quad precision backend

Windows PowerShell example:

```powershell
$env:SPHEROIDALWAVEFUNCTIONS_LIBRARY_R8 = "C:\path\to\spheroidal_batch_r8.dll"
$env:SPHEROIDALWAVEFUNCTIONS_LIBRARY_R16 = "C:\path\to\spheroidal_batch_r16.dll"
julia
```

Linux/macOS shell example:

```bash
export SPHEROIDALWAVEFUNCTIONS_LIBRARY_R8="/path/to/libspheroidal_batch_r8.so"
export SPHEROIDALWAVEFUNCTIONS_LIBRARY_R16="/path/to/libspheroidal_batch_r16.so"
julia
```

## Programmatic Override

Use explicit runtime configuration when paths are known inside application code:

```julia
using SpheroidalWaveFunctions

set_backend_library!("/path/to/libspheroidal_batch_r8.so"; precision=:double)
set_backend_library!("/path/to/libspheroidal_batch_r16.so"; precision=:quad)
```

## Inspect Active Paths

```julia
backend_library(precision=:double)
backend_library(precision=:quad)
```

## Notes

- `deps/library_config.jl` is generated locally by build tooling and is intentionally gitignored.
- If a configured path does not exist, it is ignored and a warning is emitted.
- If no backend is configured for a requested precision, calls fail with a clear error message.
