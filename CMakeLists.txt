# cheatah-space — the first cheatah standard-library extension.
#
# Provides the `space` package (astronomy / space-physics) as hand-authored, header-only C++
# (templated with concepts), in the `cheatah::space::*` namespaces. A downstream project opts
# in with `biome add cheatah-space`, then `import space.time` (etc.).
#
# ZERO DEPENDENCIES — like cheatah itself, this extension pulls in nothing external. The only
# build input is the cheatah toolchain (for purrc + the stdlib headers the dev-time C++ tests
# compile against). NASA CDF and PRBEM/IRBEM are NEVER required: they are optional, dev-only
# references used to verify and benchmark space.cdf / space.irbem, gated behind a flag.
#
# Modules (see each subdir README for status & roadmap):
#   space.time   — Julian dates & epoch conversions (incl. the CDF_EPOCH bridge)   [WORKING]
#   space.cdf    — NASA Common Data Format I/O, from scratch in C++                 [roadmap]
#   space.irbem  — radiation-belt / magnetic-field models, reimplemented clean-room [roadmap]
cmake_minimum_required(VERSION 3.24)
project(cheatah-space LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- The `space` package, surfaced as a header-only INTERFACE target -----------------------
# Consumers resolve `import space.time` from these headers on purrc's import path; nothing to
# compile or link. The include dir is the repo root so `space/space.hpp` & `space/<mod>/...`
# resolve exactly as purrc expects.
add_library(cheatah_space INTERFACE)
add_library(cheatah::space ALIAS cheatah_space)
target_include_directories(cheatah_space INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_compile_features(cheatah_space INTERFACE cxx_std_20)

# Everything below is for building/testing THIS repo standalone; when cheatah-space is fetched
# as a sub-project (the biome path) the parent only wants the headers above.
if(NOT PROJECT_IS_TOP_LEVEL)
    return()
endif()

# A single, predictable output tree (so scripts/coverage.sh finds ./bin/cheatah_space_tests).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Floating point: the contract, not a preference -------------------------------------------
# space.irbem's reference lane must produce the SAME BITS at every optimization level, because its
# self-goldens are byte comparisons and its differential tests against the IRBEM oracle are read as
# physics, not as noise. FMA contraction silently changes the last bit of `a*b+c`, so it is off —
# and it is set GLOBALLY rather than per target, because a flag that must be remembered at each
# add_executable() is a flag that will be forgotten at one of them, and one target compiled
# differently is the whole divergence back again.
# -ffast-math is never used anywhere in this repo.
add_compile_options(-ffp-contract=off)

# --- Standalone dev: wire git hooks so a push runs the qa_gate ----------------------------
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git" OR EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../.git")
    execute_process(
        COMMAND git config core.hooksPath .githooks
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_QUIET ERROR_QUIET)
    message(STATUS "cheatah-space: git pre-push hook -> scripts/qa_gate.sh")
endif()

# The sibling checkout the qa_gate + dev-time tests need: the cheatah toolchain (purrc + the
# cheatah runtime + the stdlib sources the C++ tests compile). Override to relocate it.
set(CHEATAH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../cheatah" CACHE PATH "cheatah toolchain checkout")

# --- Options -----------------------------------------------------------------------------------
option(CHEATAH_SPACE_BUILD_TESTS "Build the C++ unit tests (fetches GoogleTest)" OFF)
option(CHEATAH_SPACE_SANITIZE    "Build the unit tests with ASan + UBSan" OFF)
option(CHEATAH_SPACE_GPU         "Build space.irbem's GPU lane (needs slangc + cheatah-gpu-linalg)" ON)

# --- space.irbem's GPU lane ---------------------------------------------------------------------
# space.irbem's field-line and drift-shell integrals run on the GPU. The kernels are written once in
# Slang and compiled to SPIR-V by slangc at build time — the same arrangement cheatah-gpu-linalg
# uses for kernels/linalg.slang, and the build rule below is deliberately its shape.
#
# Every kernel binds set 0 from slot 0 (the ABI the contexts implement), so two kernels visible in
# one compilation would overlap. Each is therefore guarded by -DKERNEL_<NAME> and compiled alone —
# which is why this is a table of (entry, guard) rather than one invocation over the file.
#
# OFF by default: the CPU lane is complete on its own, and a researcher without a Vulkan SDK must
# still be able to build and use space.irbem.
if(CHEATAH_SPACE_GPU)
    set(CHEATAH_GPU_LINALG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../cheatah-gpu-linalg"
        CACHE PATH "cheatah-gpu-linalg checkout (provides the device Context and the consumer seam)")
    # The GPU stack is CONSUMED, never reimplemented: this include brings in the include dirs, the
    # volk archive and the shader dir, and defines cheatah_gpu_linalg::consumer. Nothing here names
    # a backend or links -lvulkan.
    # ON by default, so a missing toolchain must DEGRADE rather than abort — a researcher with no
    # Vulkan SDK still gets the entire CPU library, which is complete on its own. Only an EXPLICIT
    # -DCHEATAH_SPACE_GPU=ON with a missing dependency is an error, because then the user asked for
    # something we cannot deliver and silence would be the wrong answer.
    find_program(CHEATAH_SPACE_SLANGC slangc
        HINTS "$ENV{VULKAN_SDK}/bin" "$ENV{HOME}/Tools/vulkan-sdk/1.4.341.1/x86_64/bin")
    set(_gpu_missing "")
    if(NOT EXISTS "${CHEATAH_GPU_LINALG_DIR}/build/consumer.cmake")
        set(_gpu_missing "${_gpu_missing} cheatah-gpu-linalg (${CHEATAH_GPU_LINALG_DIR}/build/consumer.cmake)")
    endif()
    if(NOT CHEATAH_SPACE_SLANGC)
        set(_gpu_missing "${_gpu_missing} slangc (ships in the Vulkan SDK)")
    endif()
    if(_gpu_missing)
        if(DEFINED CACHE{CHEATAH_SPACE_GPU} AND CHEATAH_SPACE_GPU_EXPLICIT)
            message(FATAL_ERROR "CHEATAH_SPACE_GPU was requested but is missing:${_gpu_missing}")
        endif()
        message(WARNING
            "cheatah-space: GPU lane DISABLED — missing:${_gpu_missing}. The CPU lane is complete; "
            "space.irbem will report itself GPU-unavailable at runtime rather than failing.")
        set(CHEATAH_SPACE_GPU OFF)
    endif()
endif()

if(CHEATAH_SPACE_GPU)
    include("${CHEATAH_GPU_LINALG_DIR}/build/consumer.cmake")

    set(_irbem_slang "${CMAKE_CURRENT_SOURCE_DIR}/space/irbem/gpu/irbem.slang")
    set(_irbem_spv_dir "${CMAKE_BINARY_DIR}/shaders")
    file(MAKE_DIRECTORY "${_irbem_spv_dir}")

    # "entry|GUARD" — one of the THREE places a kernel must be registered, with dispatch.hpp's
    # registered_kernels and irbem.slang's entry point. RegistryMatchesTheSlangEntryPointsAndThe
    # CMakeShaderTable asserts all three lists agree. Omitting a row HERE is the silent one: slangc
    # is never asked to compile the kernel, so shader_path() never exists, every existence guard
    # takes the host lane, and the batch reports Result<bool>::value == false forever.
    set(_irbem_kernels
        "irbem_dipole_f32|KERNEL_DIPOLE"
        "irbem_igrf_f32|KERNEL_IGRF"
                "irbem_trace_i_f32|KERNEL_TRACE"
        "irbem_shell_foot_f32|KERNEL_SHELL"
        "irbem_t89_f32|KERNEL_T89"
        "irbem_trace_path_f32|KERNEL_TRACE_PATH"
        "irbem_trace_total_f32|KERNEL_TRACE_TOTAL"
        "irbem_mead_f32|KERNEL_MEAD"
        "irbem_opq_f32|KERNEL_OPQ"
        "irbem_opd_f32|KERNEL_OPD"
        "irbem_om97_f32|KERNEL_OSTAPENKO")

    set(_irbem_spv_outputs "")
    foreach(_spec IN LISTS _irbem_kernels)
        string(REPLACE "|" ";" _parts "${_spec}")
        list(GET _parts 0 _entry)
        list(GET _parts 1 _guard)
        set(_spv "${_irbem_spv_dir}/${_entry}.spv")
        add_custom_command(OUTPUT "${_spv}"
            COMMAND "${CHEATAH_SPACE_SLANGC}" "${_irbem_slang}"
                    -D${_guard} -target spirv -entry ${_entry} -stage compute -o "${_spv}"
            DEPENDS "${_irbem_slang}"
            COMMENT "slangc ${_entry}.spv"
            VERBATIM)
        list(APPEND _irbem_spv_outputs "${_spv}")
    endforeach()
    add_custom_target(cheatah_space_shaders ALL DEPENDS ${_irbem_spv_outputs})

    # Consumers of the GPU lane get the seam plus the shader location. The runtime env var
    # CHEATAH_SPACE_IRBEM_SPV_DIR still overrides this, which is how a test relocates the blobs.
    add_library(cheatah_space_gpu INTERFACE)
    add_library(cheatah::space_gpu ALIAS cheatah_space_gpu)
    target_link_libraries(cheatah_space_gpu INTERFACE cheatah::space cheatah_gpu_linalg::consumer)
    target_compile_definitions(cheatah_space_gpu INTERFACE
        CHEATAH_SPACE_IRBEM_SPV_DIR="${_irbem_spv_dir}")
    message(STATUS "cheatah-space: GPU lane ON — slangc=${CHEATAH_SPACE_SLANGC}, shaders=${_irbem_spv_dir}")
endif()

# --- C++ unit tests: GoogleTest, fetched ONLY under this option --------------------------------
include(FetchContent)
if(CHEATAH_SPACE_BUILD_TESTS)
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        52eb8108c5bdec04579160ae17225d66034bd723)  # pinned (matches ../cheatah)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    enable_testing()
    add_subdirectory(tests)
endif()

# --- The cheatah (.purr) system suite, run through the qa_gate as one ctest ----------------
# The gate compiles each systests/test_*.purr against the `space` package with purrc and runs
# it, which also gates that the C++ templates instantiate cleanly (scalar + ndarray paths).
include(CTest)
if(BUILD_TESTING)
    add_test(NAME qa_gate
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/qa_gate.sh
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(qa_gate PROPERTIES ENVIRONMENT "CHEATAH_DIR=${CHEATAH_DIR}")
endif()
