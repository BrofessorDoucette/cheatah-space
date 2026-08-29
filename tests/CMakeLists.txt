# C++ unit tests for the space package's hand-authored headers.
#
# The single binary `cheatah_space_tests` is what scripts/coverage.sh instruments (clang source-based
# coverage) and what security/run-valgrind.sh memchecks. Linking the cheatah::space INTERFACE target
# puts the repo root on the include path so `#include "space/space.hpp"` resolves.
#
# space.time is pure arithmetic over the cheatah prelude (`cheatah.hpp`) and — on the vectorized
# path — the `ndarray` module. The prelude and ndarray's arithmetic operators are header-only, but
# ndarray's shape machinery (broadcast_shapes) is compiled, so we COMPILE the module's sources into
# our own STATIC target from the sibling ../cheatah checkout (CHEATAH_DIR) — the standard way a C++
# program consumes a cheatah stdlib module. Building from source means the stdlib is compiled with
# THIS build's compiler and flags — so the ASan preset sanitizes it and the coverage build
# instruments it (no archive/config mismatch), and the tests need only a cheatah SOURCE checkout.
#
# (Consumers of cheatah-space itself never do any of this: they `biome add cheatah-space` and
# `import space.time` — biome wires CHEATAH_MODULE_PATH and the stdlib. This is only cheatah-space's
# own dev-time C++ coverage harness for the headers; see scripts/test-biome-install.sh for the real
# consumer path.)
include(GoogleTest)

set(_cheatah_stdlib ${CHEATAH_DIR}/stdlib)
if(NOT EXISTS ${_cheatah_stdlib}/ndarray/ndarray.cpp)
    message(FATAL_ERROR
        "cheatah-space tests: ${_cheatah_stdlib} is not a cheatah checkout — set CHEATAH_DIR to your "
        "cheatah source tree (default ../cheatah).")
endif()

# The one stdlib module with compiled symbols the vectorized space.time path touches. Its public
# include dirs carry both the ndarray surface and the builtins prelude (cheatah.hpp). Guarded so a
# future consumer target in this repo can reuse it without redefinition.
if(NOT TARGET cheatah_ndarray)
    add_library(cheatah_ndarray STATIC
        ${_cheatah_stdlib}/ndarray/ndarray.cpp
        ${_cheatah_stdlib}/ndarray/ufunc_simd.cpp
    )
    target_include_directories(cheatah_ndarray PUBLIC
        ${_cheatah_stdlib}/ndarray
        ${_cheatah_stdlib}/builtins
        ${_cheatah_stdlib}/fixarray
    )
    target_compile_features(cheatah_ndarray PUBLIC cxx_std_20)
    set_target_properties(cheatah_ndarray PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

add_executable(cheatah_space_tests
    alloc_counter.cpp
    space_time_test.cpp
    space_civil_test.cpp
    space_calendar_reform_test.cpp
    space_cdf_leapseconds_test.cpp
    space_cdf_format_test.cpp
    space_cdf_records_test.cpp
    space_cdf_reader_test.cpp
    space_cdf_mapping_test.cpp
    irbem_frames_test.cpp
    irbem_policy_test.cpp
    irbem_context_test.cpp
    irbem_datetime_test.cpp
    irbem_coords_rotations_test.cpp
    irbem_coords_geodetic_test.cpp
    irbem_coords_helio_test.cpp
    irbem_igrf_test.cpp
    irbem_status_test.cpp
    irbem_api_test.cpp
    irbem_field_test.cpp
    irbem_field_hostpaths_test.cpp
    irbem_t89_test.cpp
    irbem_mead_test.cpp
    irbem_opq_test.cpp
    irbem_opd_test.cpp
    irbem_ostapenko_test.cpp
    irbem_batch_simd_test.cpp
    irbem_trace_api_test.cpp
    irbem_driftshell_test.cpp
    irbem_total_field_test.cpp
)
target_link_libraries(cheatah_space_tests PRIVATE
    cheatah::space
    cheatah_ndarray
    GTest::gtest_main
)

# Optional ASan + UBSan (the `asan` preset turns this on). Sanitize the stdlib we compiled too, so
# the whole binary is instrumented consistently. Valgrind runs the plain debug build.
if(CHEATAH_SPACE_SANITIZE)
    foreach(_t cheatah_space_tests cheatah_ndarray)
        target_compile_options(${_t} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    endforeach()
    target_link_options(cheatah_space_tests PRIVATE -fsanitize=address,undefined)
endif()

gtest_discover_tests(cheatah_space_tests)

# --- space.irbem's GPU lane -----------------------------------------------------------------
# Its own target, so `ctest -L gpu` selects the device tests and nothing else — labelling the whole
# suite would make the label mean nothing. Built only with CHEATAH_SPACE_GPU=ON; the seam itself is
# already covered without a device by the __has_include path in the main suite.
if(CHEATAH_SPACE_GPU)
    add_executable(cheatah_space_gpu_tests irbem_gpu_test.cpp)
    target_link_libraries(cheatah_space_gpu_tests PRIVATE
        cheatah::space_gpu
        cheatah_ndarray
        GTest::gtest_main)
    add_dependencies(cheatah_space_gpu_tests cheatah_space_shaders)
    # Vendor driver caches leak by design; leak-checking the GPU lane reports the driver, not us.
    gtest_discover_tests(cheatah_space_gpu_tests
        TEST_PREFIX "gpu:"
        PROPERTIES LABELS gpu ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
endif()
