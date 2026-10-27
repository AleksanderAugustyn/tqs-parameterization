cmake_minimum_required(VERSION 3.20)

include(FetchContent)

project(tqs-parameterization
        VERSION 0.1.0
        LANGUAGES Fortran CXX
        DESCRIPTION "Nix three-quadratic-surface (3QS) nuclear shape parameterization"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Fetch GCC compiler options package ---
FetchContent_Declare(
        gcc-compiler-options
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/gcc-compiler-options.git
        GIT_TAG 2.0.0
)
FetchContent_MakeAvailable(gcc-compiler-options)

# --- Fetch foundational Fortran modules ---
FetchContent_Declare(
        fortran-foundations
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/fortran-foundations.git
        GIT_TAG 3.0.0
)
FetchContent_MakeAvailable(fortran-foundations)

# --- Preprocessor & build type ---
set(CMAKE_Fortran_PREPROCESS ON)
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "RelWithDebInfo" "Release")
endif ()

# --- Compiler flags via gcc-compiler-options ---
include(GCCCompilerOptions/FortranCompilerOptions)
create_fortran_library_interface(TARGET tqs_parameterization_flags)

# --- Source files ---
set(TQS_PARAM_SOURCES
        src/tqs_workers_mod.f08
        src/tqs_resolve_mod.f08
        src/tqs_parameterization_mod.f08
        src/tqs_parameterization_c_api_mod.f08
)

# --- No LTO for the library's own objects ---
# One-shot == cached and repeat == first call are BITWISE claims. The
# arithmetic lives in the workers and resolve modules and the tiers in the
# front module; without LTO no kernel can be inlined into a tier routine, so
# every tier path executes the same machine code. With -flto=auto -ffast-math
# a kernel is inlined into each caller and optimized per call site, and the
# same call then returns different last bits from different call sites. A
# source-file property lands after the gcc-opts interface flags, so it wins.
# Consumers keep LTO for their own code.
set_source_files_properties(${TQS_PARAM_SOURCES} PROPERTIES COMPILE_OPTIONS "-fno-lto")

# --- Static library target ---
add_library(tqs_parameterization STATIC ${TQS_PARAM_SOURCES})
target_link_libraries(tqs_parameterization
        PUBLIC fortran_foundations
        PRIVATE tqs_parameterization_flags
)
target_compile_definitions(tqs_parameterization PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO>
        $<$<CONFIG:Release>:RELEASE>
)
set_target_properties(tqs_parameterization PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)
target_include_directories(tqs_parameterization
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
add_library(TqsParameterization::tqs_parameterization ALIAS tqs_parameterization)

# --- Shared library target (for FFI consumers: Python ctypes, etc.) ---
add_library(tqs_parameterization_shared SHARED ${TQS_PARAM_SOURCES})
target_link_libraries(tqs_parameterization_shared
        PUBLIC fortran_foundations
        PRIVATE tqs_parameterization_flags
)
target_compile_definitions(tqs_parameterization_shared PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO>
        $<$<CONFIG:Release>:RELEASE>
)
set_target_properties(tqs_parameterization_shared PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_shared
        OUTPUT_NAME tqs_parameterization
)
target_include_directories(tqs_parameterization_shared
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files_shared>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
add_library(TqsParameterization::tqs_parameterization_shared
        ALIAS tqs_parameterization_shared)

# --- C++ interface target (consumes the shared library via the C header) ---
add_library(tqs_parameterization_cxx INTERFACE)
target_include_directories(tqs_parameterization_cxx INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(tqs_parameterization_cxx INTERFACE
        tqs_parameterization_shared)
target_compile_features(tqs_parameterization_cxx INTERFACE cxx_std_20)
add_library(TqsParameterization::tqs_parameterization_cxx
        ALIAS tqs_parameterization_cxx)

# --- Install the shared library into the Python package (wheel builds) ---
install(TARGETS tqs_parameterization_shared
        LIBRARY DESTINATION tqs_parameterization)

# --- Tests (skipped for wheel builds via -DTQS_PARAM_BUILD_TESTS=OFF) ---
option(TQS_PARAM_BUILD_TESTS "Build the Fortran/C++ test suite" ON)
if (TQS_PARAM_BUILD_TESTS)
    enable_testing()

add_library(tqs_param_test_utils STATIC
        tests/test_utils_mod.f08
        tests/tqs_test_fixtures_mod.f08
        tests/tqs_test_tiers_mod.f08
)
target_link_libraries(tqs_param_test_utils
        PUBLIC tqs_parameterization fortran_foundations
        PRIVATE tqs_parameterization_flags)
set_target_properties(tqs_param_test_utils PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests)
target_include_directories(tqs_param_test_utils PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests>)

set(TQS_PARAM_TEST_SUITES status workers resolve core neck moments rtheta error boundary equivalence statelessness breakdown)
foreach (suite IN LISTS TQS_PARAM_TEST_SUITES)
    add_executable(tqs_param_${suite}_test tests/tqs_param_${suite}_test.f08)
    target_link_libraries(tqs_param_${suite}_test PRIVATE
            tqs_parameterization tqs_param_test_utils tqs_parameterization_flags)
    set_target_properties(tqs_param_${suite}_test PROPERTIES
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests/${suite})
    add_test(NAME ${suite} COMMAND tqs_param_${suite}_test)
endforeach ()

add_executable(tqs_param_c_api_smoke_test tests/c_api_smoke_test.cpp)
target_link_libraries(tqs_param_c_api_smoke_test PRIVATE tqs_parameterization_cxx)
add_test(NAME c_api_smoke COMMAND tqs_param_c_api_smoke_test)

# Contract family 3, concurrency: many threads computing on ONE shared cache.
find_package(Threads REQUIRED)
add_executable(tqs_param_thread_stress_test tests/thread_stress_test.cpp)
target_link_libraries(tqs_param_thread_stress_test PRIVATE
        tqs_parameterization_cxx Threads::Threads)
add_test(NAME thread_stress COMMAND tqs_param_thread_stress_test)
endif ()
