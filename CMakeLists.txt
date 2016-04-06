cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)
project(pbrt-v2)

add_definitions(-DPBRT_HAS_OPENEXR -DPBRT_PROBES_NONE)

# turn on c++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Qt
set(CMAKE_AUTOMOC ON)
find_package(Qt5Core)
find_package(Qt5Network)

# OpenEXR
include_directories(/usr/include/OpenEXR)
set(OPENEXR_LIBRARIES IlmImf Half)


# Threads
set(CMAKE_THREAD_PREFER_PTHREAD)
find_package(Threads)

include_directories(src src/core)

# pbrt lib
file(GLOB CORE_FILES
    "src/core/*.cpp"
    "src/core/*.c"
    "src/core/*.h"
    "src/accelerators/*.h"
    "src/accelerators/*.cpp"
    "src/cameras/*.h"
    "src/cameras/*.cpp"
    "src/film/*.h"
    "src/film/*.cpp"
    "src/filters/*.h"
    "src/filters/*.cpp"
    "src/integrators/*.h"
    "src/integrators/*.cpp"
    "src/lights/*.h"
    "src/lights/*.cpp"
    "src/materials/*.h"
    "src/materials/*.cpp"
    "src/renderers/*.h"
    "src/renderers/*.cpp"
    "src/samplers/*.h"
    "src/samplers/*.cpp"
    "src/shapes/*.h"
    "src/shapes/*.cpp"
    "src/textures/*.h"
    "src/textures/*.cpp"
    "src/volumes/*.h"
    "src/volumes/*.cpp"
)
add_library(core ${CORE_FILES})
target_link_libraries(core RenderingServer ${CMAKE_THREAD_LIBS_INIT} ${OPENEXR_LIBRARIES})

# pbrt app
file(GLOB PBRT_FILES
    "src/main/*.h"
    "src/main/*.cpp"
)
add_executable(pbrt ${PBRT_FILES})
target_link_libraries(pbrt core)