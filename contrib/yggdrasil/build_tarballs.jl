using BinaryBuilder

name = "SpheroidalWaves"
version = v"0.4.1"

# For a local BinaryBuilder trial, point this at the source checkout. For the
# Yggdrasil submission, replace this selection with a GitSource pinned to the
# published commit containing the backend changes. Never build a moving branch.
local_source = get(ENV,"SPHEROIDALWAVES_SOURCE_DIR","")
revision = get(ENV,"SPHEROIDALWAVES_SOURCE_REVISION","")
sources = if isempty(local_source)
    occursin(r"^[0-9a-f]{40}$",revision) ||
        error("Set SPHEROIDALWAVES_SOURCE_REVISION to the published backend commit (40 hex digits)")
    [GitSource("https://github.com/brandynlucca/SpheroidalWaves.jl.git",revision)]
else
    [DirectorySource(abspath(local_source);target="SpheroidalWaves")]
end

script = raw"""
cd ${WORKSPACE}/srcdir/SpheroidalWaves*
cmake -S . -B build-binarybuilder \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_BUILD_TYPE=Release
cmake --build build-binarybuilder --parallel ${nproc}
cmake --install build-binarybuilder
install_license LICENSE
"""

# Start with the platforms currently distributed by the package. Quad precision
# and the Fortran/OpenMP runtime must be validated on each before expanding it.
platforms = expand_gfortran_versions([
    Platform("x86_64","linux";libc="glibc"),
    Platform("x86_64","windows"),
    Platform("x86_64","macos"),
    Platform("aarch64","macos"),
])

products = [
    LibraryProduct(["libspheroidal_batch_double","spheroidal_batch_double"],:libspheroidal_batch_double),
    LibraryProduct(["libspheroidal_batch_quad","spheroidal_batch_quad"],:libspheroidal_batch_quad),
]

dependencies = [Dependency("CompilerSupportLibraries_jll")]

build_tarballs(ARGS,name,version,sources,script,platforms,products,dependencies;
               preferred_gcc_version=v"12",julia_compat="1.10")
