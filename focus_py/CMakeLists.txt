cmake_minimum_required(VERSION 3.15...3.27)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find Python and pybind11
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Find Qhull - try config mode first, then module mode
find_package(Qhull CONFIG QUIET)
if(NOT Qhull_FOUND)
    find_package(Qhull MODULE QUIET)
endif()

if(NOT Qhull_FOUND)
    # Try to find it manually as a fallback
    find_path(QHULL_INCLUDE_DIR 
        NAMES libqhullcpp/Qhull.h
        PATHS 
            /usr/include
            /usr/local/include
            /opt/homebrew/include  # macOS Homebrew
            $ENV{CONDA_PREFIX}/include  # Conda environments
    )
    
    find_library(QHULLCPP_LIBRARY
        NAMES qhullcpp
        PATHS
            /usr/lib
            /usr/local/lib
            /usr/lib/x86_64-linux-gnu
            /opt/homebrew/lib
            $ENV{CONDA_PREFIX}/lib
    )
    
    find_library(QHULL_R_LIBRARY
        NAMES qhull_r
        PATHS
            /usr/lib
            /usr/local/lib
            /usr/lib/x86_64-linux-gnu
            /opt/homebrew/lib
            $ENV{CONDA_PREFIX}/lib
    )
    
    if(QHULL_INCLUDE_DIR AND QHULLCPP_LIBRARY AND QHULL_R_LIBRARY)
        set(QHULL_FOUND TRUE)
        set(QHULL_INCLUDE_DIRS ${QHULL_INCLUDE_DIR})
        set(QHULL_LIBRARIES ${QHULLCPP_LIBRARY} ${QHULL_R_LIBRARY})
        message(STATUS "Found Qhull manually:")
        message(STATUS "  Include: ${QHULL_INCLUDE_DIRS}")
        message(STATUS "  Libraries: ${QHULL_LIBRARIES}")
    else()
        message(FATAL_ERROR 
            "Qhull not found. Please install Qhull:\n"
            "  Ubuntu/Debian: sudo apt-get install libqhull-dev\n"
            "  Fedora/RHEL:   sudo dnf install qhull-devel\n"
            "  macOS:         brew install qhull\n"
            "  Conda:         conda install -c conda-forge qhull"
        )
    endif()
else()
    # Use the found Qhull package
    if(TARGET Qhull::qhullcpp)
        set(QHULL_LIBRARIES Qhull::qhullcpp Qhull::qhull_r)
    else()
        set(QHULL_LIBRARIES ${QHULL_LIBRARIES})
    endif()
    message(STATUS "Found Qhull via find_package")
endif()

# Create the Python extension module
pybind11_add_module(_focus MODULE src/pybind11_focus.cpp)

# Set include directories
target_include_directories(_focus PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${QHULL_INCLUDE_DIRS}
)

# Link Qhull libraries
target_link_libraries(_focus PRIVATE ${QHULL_LIBRARIES})

# Optimization flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

target_compile_options(_focus PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra>
)

# Install the extension module
install(TARGETS _focus LIBRARY DESTINATION focus_py)
