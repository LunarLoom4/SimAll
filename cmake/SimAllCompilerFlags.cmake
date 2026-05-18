# =============================================================================
# SimAllCompilerFlags.cmake
# Centralized compiler/linker flag policy. All subsystem libraries link the
# INTERFACE target `SimAll::CompilerFlags` to inherit these flags consistently.
# =============================================================================

if(TARGET SimAll_CompilerFlags)
    return()
endif()

add_library(SimAll_CompilerFlags INTERFACE)
add_library(SimAll::CompilerFlags ALIAS SimAll_CompilerFlags)

# ---- Warnings -----------------------------------------------------------------
if(MSVC)
    target_compile_options(SimAll_CompilerFlags INTERFACE
        /W4 /permissive- /Zc:__cplusplus /Zc:preprocessor /MP /utf-8
        /wd4251   # DLL-interface warning for STL members
        /wd4275)
    target_compile_definitions(SimAll_CompilerFlags INTERFACE
        _CRT_SECURE_NO_WARNINGS
        NOMINMAX
        WIN32_LEAN_AND_MEAN)
else()
    target_compile_options(SimAll_CompilerFlags INTERFACE
        -Wall -Wextra -Wpedantic
        -Wno-unused-parameter
        -fvisibility=hidden -fvisibility-inlines-hidden)
endif()

# ---- Vectorization ------------------------------------------------------------
if(SIMALL_USE_AVX512)
    if(MSVC)
        target_compile_options(SimAll_CompilerFlags INTERFACE /arch:AVX512)
    else()
        target_compile_options(SimAll_CompilerFlags INTERFACE -mavx512f -mavx512dq)
    endif()
elseif(SIMALL_USE_AVX2)
    if(MSVC)
        target_compile_options(SimAll_CompilerFlags INTERFACE /arch:AVX2)
    else()
        target_compile_options(SimAll_CompilerFlags INTERFACE -mavx2 -mfma)
    endif()
endif()

# ---- Release tuning -----------------------------------------------------------
if(NOT MSVC)
    target_compile_options(SimAll_CompilerFlags INTERFACE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:RelWithDebInfo>:-O2 -g>
        $<$<CONFIG:Debug>:-O0 -g>)
endif()

# ---- Sanitizers (Debug only) --------------------------------------------------
if(SIMALL_ENABLE_SANITIZERS AND NOT MSVC)
    target_compile_options(SimAll_CompilerFlags INTERFACE
        $<$<CONFIG:Debug>:-fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer>)
    target_link_options(SimAll_CompilerFlags INTERFACE
        $<$<CONFIG:Debug>:-fsanitize=address -fsanitize=undefined>)
endif()

target_compile_features(SimAll_CompilerFlags INTERFACE cxx_std_20)
