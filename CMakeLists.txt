cmake_minimum_required(VERSION 3.10)
project(EntanglementCpp LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Try to find OpenMP
find_package(OpenMP)

if(NOT OpenMP_CXX_FOUND)
    # Try Homebrew libomp locations (Apple Silicon and Intel)
    find_library(_LIBOMP_LIB NAMES omp libomp
        PATHS /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib
        NO_DEFAULT_PATH)
    find_path(_LIBOMP_INCLUDE NAMES omp.h
        PATHS /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include
        NO_DEFAULT_PATH)

    if(_LIBOMP_LIB AND _LIBOMP_INCLUDE)
        message(STATUS "Found libomp: ${_LIBOMP_LIB}, includes: ${_LIBOMP_INCLUDE}")

        # Create an imported INTERFACE target similar to OpenMP::OpenMP_CXX
        add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)

        set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${_LIBOMP_INCLUDE}"
            INTERFACE_LINK_LIBRARIES "${_LIBOMP_LIB}"
        )

        # Compiler and linker flags differ on Apple Clang
        if(APPLE)
            # Apple Clang needs the preprocessor flag to enable OpenMP support
            set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                INTERFACE_COMPILE_OPTIONS "-Xpreprocessor;-fopenmp"
            )
        else()
            set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                INTERFACE_COMPILE_OPTIONS "-fopenmp"
            )
        endif()

        # Mark OpenMP as found for the rest of the script
        set(OpenMP_CXX_FOUND TRUE)
    else()
        message(STATUS "libomp not found in Homebrew prefix (/opt/homebrew/opt/libomp).")
    endif()
endif()

# Add executable
add_executable(linking_number linking_number.cpp)

# Enable debug information and warnings (avoid passing GCC/Clang flags to MSVC)
target_compile_options(linking_number PRIVATE
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wextra>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-g>
)

# Link OpenMP if found
if(OpenMP_CXX_FOUND)
    target_link_libraries(linking_number PRIVATE OpenMP::OpenMP_CXX)
else()
    message(WARNING "OpenMP not found. Building without OpenMP. To enable OpenMP on macOS, install libomp (brew install libomp) and rerun CMake.")
endif()

# Link math library on Unix systems
if(UNIX)
    target_link_libraries(linking_number PRIVATE m)
endif()
