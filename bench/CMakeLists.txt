# bench/ — space.irbem's Google Benchmark suite: the one place every performance claim comes from.
#
# OFF by default and never part of a normal build. A benchmark binary is not a test: it takes
# minutes, it wants an idle machine, and building it into every configure would put a several-second
# FetchContent of google/benchmark in the path of anyone who just wanted the headers.
#
#   cmake -S bench -B build/bench -DCMAKE_BUILD_TYPE=Release -DCHEATAH_SPACE_BENCH=ON
#   cmake --build build/bench -j
#
# `scripts/bench_run.sh` is the supported entry point and does not need CMake at all — it compiles
# the same one translation unit with the same flags directly, which is why the flag list below and
# the one in that script must agree. They are checked against each other by `bench_run.sh check`.
cmake_minimum_required(VERSION 3.20)
project(cheatah_space_bench CXX)

option(CHEATAH_SPACE_BENCH "Build space.irbem's benchmark suite (fetches google/benchmark)" OFF)
if(NOT CHEATAH_SPACE_BENCH)
    message(STATUS "cheatah-space bench: OFF (set -DCHEATAH_SPACE_BENCH=ON to build it)")
    return()
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

get_filename_component(_repo_root ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)
set(CHEATAH_DIR "${_repo_root}/../cheatah" CACHE PATH "cheatah source checkout")
set(CHEATAH_GPU_LINALG_DIR "${_repo_root}/../cheatah-gpu-linalg" CACHE PATH "cheatah-gpu-linalg checkout")

# Same pin as cheatah-gpu-linalg/bench and cheatah's tests/benchmarks. A benchmark harness that
# drifts between repos produces numbers that cannot be compared between repos.
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
include(FetchContent)
FetchContent_Declare(benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        192ef10025eb2c4cdd392bc502f0c852196baa48)
FetchContent_MakeAvailable(benchmark)

add_executable(cheatah_space_bench ${CMAKE_CURRENT_SOURCE_DIR}/irbem_bench.cpp)

# -O3 -march=native, and -ffp-contract=off so the host lane's arithmetic is the arithmetic the
# golden tests pinned. NEVER -ffast-math: it would let the compiler reassociate the spherical
# harmonic sums and the RK4 accumulation, and the resulting binary would be measuring different
# mathematics from the one the correctness suite verified.
target_compile_options(cheatah_space_bench PRIVATE -O3 -march=native -ffp-contract=off)
target_include_directories(cheatah_space_bench PRIVATE
    ${_repo_root}
    ${CHEATAH_DIR}/stdlib/ndarray
    ${CHEATAH_DIR}/stdlib/builtins
    ${CHEATAH_DIR}/stdlib/fixarray)
target_link_libraries(cheatah_space_bench PRIVATE benchmark::benchmark ${CMAKE_DL_LIBS})

# The device lane. Absent, `__has_include` in dispatch.hpp compiles the whole suite host-only and
# the GPU rows report themselves skipped — which is a truthful table on a machine with no GPU, not
# a broken build.
if(EXISTS ${CHEATAH_GPU_LINALG_DIR}/consumer.cmake)
    include(${CHEATAH_GPU_LINALG_DIR}/consumer.cmake)
    target_link_libraries(cheatah_space_bench PRIVATE cheatah_gpu_linalg::consumer)
    target_compile_definitions(cheatah_space_bench PRIVATE
        CHEATAH_SPACE_IRBEM_SPV_DIR="${_repo_root}/build/gpu/shaders")
else()
    message(STATUS "cheatah-space bench: no cheatah-gpu-linalg at ${CHEATAH_GPU_LINALG_DIR} — "
                   "host lanes only")
endif()
