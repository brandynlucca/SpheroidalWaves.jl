using SHA

build = abspath(only(ARGS))
inputs = IOBuffer()
for name in ("ImageOS", "ImageVersion", "RUNNER_ARCH")
    println(inputs, name, '=', get(ENV, name, ""))
end
cache = Dict{String,String}()
for line in eachline(joinpath(build, "CMakeCache.txt"))
    entry = match(r"^([^/#][^:]*):(BOOL|STRING|PATH|FILEPATH)=(.*)$", line)
    entry === nothing && continue
    name, _, value = entry.captures
    cache[name] = value
    println(inputs, name, '=', value)
end

function fingerprint(path)
    isfile(path) || error("Missing cache input: $path")
    println(inputs, path, '=', bytes2hex(open(sha256, path)))
end

fingerprint(@__FILE__)
sources = read(`git ls-files -z -- CMakeLists.txt deps test/CMakeLists.txt 'test/*.f90' 'test/*.F90' 'test/*.c' 'test/*.h'`, String)
for path in sort!(split(sources, '\0'; keepempty=false))
    fingerprint(path)
end
for (name, value) in sort!(collect(cache))
    if endswith(name, "_LIBRARY") || occursin(r"^CMAKE_(C|Fortran)_COMPILER$", name)
        fingerprint(value)
    end
end
for (name, header) in (("MPFR", "mpfr.h"), ("MPC", "mpc.h"), ("GMP", "gmp.h"))
    fingerprint(joinpath(cache[name * "_INCLUDE_DIR"], header))
end
for language in ("C", "Fortran")
    compiler = cache["CMAKE_$(language)_COMPILER"]
    print(inputs, read(`$compiler --version`, String))
    for library in ("libgcc.a", "libgfortran.a", "libquadmath.a")
        path = strip(read(`$compiler -print-file-name=$library`, String))
        isfile(path) && fingerprint(path)
    end
end
# Windows arithmetic DLLs are installed before restoring cached build products.
if Sys.iswindows()
    for path in sort!(readdir(joinpath(build, "bin"); join=true))
        endswith(path, ".dll") && fingerprint(path)
    end
end
key = bytes2hex(sha256(take!(inputs)))
println("Native build inputs: $key")
open(io -> println(io, "key=$key"), ENV["GITHUB_OUTPUT"], "a")
