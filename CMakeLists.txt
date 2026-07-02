# cheatah-space — the first cheatah standard-library extension.
#
# Provides the `space` package (astronomy / space-physics) as hand-authored, header-only C++
# (templated with concepts), in the `cheatah::space::*` namespaces. A downstream project opts
# in with `biome add cheatah-space`, then `import space.time` (etc.).
#
# ZERO DEPENDENCIES — like cheatah itself, this extension pulls in nothing external. The only
# build input is the cheatah toolchain (for purrc + the cmake helpers when building/testing
# standalone). NASA CDF and PRBEM/IRBEM are NEVER required: they are optional, dev-only
# references used to verify and benchmark space.cdf / space.irbem, gated behind a flag.
#
# Modules (see each subdir README for status & roadmap):
#   space.time   — Julian dates & epoch conversions (incl. the CDF_EPOCH bridge)   [WORKING]
#   space.cdf    — NASA Common Data Format I/O, from scratch in C++                 [roadmap]
#   space.irbem  — radiation-belt / magnetic-field models, reimplementing IRBEM     [roadmap]
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

# --- Standalone dev: wire git hooks so a push runs the qa_gate ----------------------------
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git" OR EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../.git")
    execute_process(
        COMMAND git config core.hooksPath .githooks
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_QUIET ERROR_QUIET)
    message(STATUS "cheatah-space: git pre-push hook -> scripts/qa_gate.sh")
endif()

# --- Tests: the cheatah (.purr) suite, run through the qa_gate as one ctest ----------------
# The gate compiles each tests/*.purr against the `space` package with purrc and runs it,
# which also gates that the C++ templates instantiate cleanly (scalar + ndarray paths).
include(CTest)
if(BUILD_TESTING)
    add_test(NAME qa_gate
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/qa_gate.sh
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    # Let the gate find the sibling toolchain; override -DCHEATAH_DIR to relocate it.
    set(CHEATAH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../cheatah" CACHE PATH "cheatah toolchain checkout")
    set_tests_properties(qa_gate PROPERTIES ENVIRONMENT "CHEATAH_DIR=${CHEATAH_DIR}")
endif()
