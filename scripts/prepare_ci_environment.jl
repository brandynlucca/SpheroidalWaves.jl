using Pkg, SHA

@assert Sys.ARCH === (ENV["EXPECTED_JULIA_ARCH"] == "x64" ? :x86_64 : :aarch64)
Pkg.instantiate(;allow_autoprecomp=false)

inputs = IOBuffer()
println(inputs, VERSION, ' ', Base.GIT_VERSION_INFO.commit, ' ', Sys.MACHINE)
println(inputs, get(ENV, "ImageVersion", ""))
println(inputs, "bounds=yes inline=yes optimize=1 cpu=", ENV["JULIA_CPU_TARGET"])
println(inputs, abspath(pwd()), ' ', DEPOT_PATH)
for path in ("Project.toml", "Manifest.toml", ".github/workflows/CI.yml",
        "scripts/prepare_ci_environment.jl", "scripts/precompile_ci_environment.jl")
    print(inputs, path, '=', bytes2hex(open(sha256, path)), '\n')
end
for (directory, _, files) in walkdir("src"), file in sort!(files)
    path = joinpath(directory, file)
    print(inputs, path, '=', bytes2hex(open(sha256, path)), '\n')
end
key = bytes2hex(sha256(take!(inputs)))
println("Julia compilation inputs: $key")
open(ENV["GITHUB_OUTPUT"], "a") do io
    println(io, "key=$key")
    println(io, "depot=", first(DEPOT_PATH))
    println(io, "version=$VERSION")
end
