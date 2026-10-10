cmake_minimum_required(VERSION 3.16)
# Run on older macOS than the build machine (must be set before project()).
set(CMAKE_OSX_DEPLOYMENT_TARGET "12.0" CACHE STRING "Minimum macOS version")
project(remersher VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Eigen (header-only): use the system copy, or download it if missing.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 not found on the system; downloading Eigen 3.4.0")
  include(FetchContent)
  # SOURCE_SUBDIR points at a folder without a CMakeLists.txt, so this only downloads (Eigen's
  # own build, tests and install rules are not needed for header-only use).
  FetchContent_Declare(eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
    SOURCE_SUBDIR do-not-build)
  FetchContent_MakeAvailable(eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen3::Eigen INTERFACE ${eigen_SOURCE_DIR})
endif()

# Boost (header-only parts: graph max-flow). Too large to download here, so explain how to get it.
find_package(Boost QUIET CONFIG)
if(NOT Boost_FOUND)
  find_package(Boost QUIET MODULE)
endif()
if(NOT TARGET Boost::headers)
  message(FATAL_ERROR "Boost headers not found. Install them with one of:\n"
    "  Debian / Ubuntu:  sudo apt install libboost-dev\n"
    "  Fedora:           sudo dnf install boost-devel\n"
    "  Arch / Manjaro:   sudo pacman -S boost\n"
    "  openSUSE:         sudo zypper install boost-devel\n"
    "  macOS:            brew install boost")
endif()

# --- QuadriFlow (BSD) ------------------------------------------------------------------------
set(QF_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/quadriflow)
if(NOT EXISTS ${QF_DIR}/src/parametrizer.cpp)
  message(FATAL_ERROR "QuadriFlow is missing: run `git submodule update --init`")
endif()

set(LEMON_DIR ${QF_DIR}/3rd/lemon-1.3.1)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)  # lemon's CMakeLists predates CMake 3.5
add_subdirectory(${LEMON_DIR} ${CMAKE_BINARY_DIR}/lemon EXCLUDE_FROM_ALL)
# lemon uses `register`, which C++17 removed.
set_target_properties(lemon PROPERTIES CXX_STANDARD 14)

add_library(quadriflow STATIC
  ${QF_DIR}/src/adjacent-matrix.cpp
  ${QF_DIR}/src/dedge.cpp
  ${QF_DIR}/src/hierarchy.cpp
  ${QF_DIR}/src/loader.cpp
  ${QF_DIR}/src/localsat.cpp
  ${QF_DIR}/src/merge-vertex.cpp
  ${QF_DIR}/src/optimizer.cpp
  ${QF_DIR}/src/parametrizer.cpp
  ${QF_DIR}/src/parametrizer-flip.cpp
  ${QF_DIR}/src/parametrizer-int.cpp
  ${QF_DIR}/src/parametrizer-mesh.cpp
  ${QF_DIR}/src/parametrizer-scale.cpp
  ${QF_DIR}/src/parametrizer-sing.cpp
  ${QF_DIR}/src/subdivide.cpp
)
target_include_directories(quadriflow SYSTEM PUBLIC
  ${QF_DIR}/src ${QF_DIR}/3rd/pcg32 ${QF_DIR}/3rd/pss
  ${LEMON_DIR} ${CMAKE_BINARY_DIR}/lemon)
# EIGEN_MPL2_ONLY keeps the build free of LGPL-only Eigen modules.
target_compile_definitions(quadriflow PUBLIC EIGEN_MPL2_ONLY $<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS>)
set_target_properties(quadriflow PROPERTIES CXX_STANDARD 14)
target_compile_options(quadriflow PRIVATE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-w>)
target_compile_options(quadriflow INTERFACE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-register> $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-deprecated-register>)
target_link_libraries(quadriflow PUBLIC Eigen3::Eigen Boost::headers lemon)

# --- remersher --------------------------------------------------------------------------------
add_library(remersher_core STATIC src/mesh_io.cpp src/mesh_cleanup.cpp src/crease.cpp src/curvature.cpp src/density.cpp src/features.cpp src/quad_cleanup.cpp src/quality.cpp src/resample.cpp src/symmetry.cpp src/remesher.cpp)
target_include_directories(remersher_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(remersher_core PRIVATE quadriflow)
target_compile_options(remersher_core PRIVATE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra>)

add_executable(remersher src/cli.cpp)
target_link_libraries(remersher PRIVATE remersher_core)

install(TARGETS remersher RUNTIME DESTINATION bin)
