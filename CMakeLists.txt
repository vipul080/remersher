cmake_minimum_required(VERSION 3.16)
project(remersher VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED CONFIG)

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
target_compile_options(quadriflow INTERFACE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-register -Wno-deprecated-register>)
target_link_libraries(quadriflow PUBLIC Eigen3::Eigen Boost::headers lemon)

# --- remersher --------------------------------------------------------------------------------
add_library(remersher_core STATIC src/mesh_io.cpp src/mesh_cleanup.cpp src/quality.cpp src/resample.cpp src/remesher.cpp)
target_include_directories(remersher_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(remersher_core PRIVATE quadriflow)
target_compile_options(remersher_core PRIVATE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra>)

add_executable(remersher src/cli.cpp)
target_link_libraries(remersher PRIVATE remersher_core)

install(TARGETS remersher RUNTIME DESTINATION bin)
