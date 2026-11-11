using SpheroidalWaves, Test
@assert isempty(SpheroidalWaves._backend_handles)

# Preserve both Julia 1.10 preference variants used by the existing test matrix.
if VERSION < v"1.11"
    preferences = "LocalPreferences.toml"
    @assert !isfile(preferences)
    try
        write(preferences, "[SpheroidalWaves]\nprecompile_workload = false\n")
        run(`$(Base.julia_cmd()) --startup-file=no --project=. --check-bounds=yes --inline=yes -O1 -e 'using SpheroidalWaves, Test'`)
    finally
        rm(preferences; force=true)
    end
end
