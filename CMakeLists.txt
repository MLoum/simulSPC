cmake_minimum_required(VERSION 3.15)
project(simulSPC VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Support OpenMP
find_package(OpenMP REQUIRED)

# Bibliothèque interne regroupant la logique métier commune
add_library(simulspc_core
    src/Experiment.cpp
    src/MDF_Gaussian.cpp
    src/Particle.cpp
    src/Solvent.cpp
)

target_include_directories(simulspc_core PUBLIC include)
target_link_libraries(simulspc_core PUBLIC OpenMP::OpenMP_CXX)

# Optimisations MSVC (/MT pour linkage statique CRT si souhaité sous Windows)
if(MSVC)
    set_property(TARGET simulspc_core PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

# ---------------------------------------------------------
# Exécutable principal
# ---------------------------------------------------------
add_executable(simulSPC src/simulSPC.cpp)
target_link_libraries(simulSPC PRIVATE simulspc_core)

# ---------------------------------------------------------
# Tests
# ---------------------------------------------------------
enable_testing()

foreach(STEP 1 2 3 4)
    set(TEST_NAME test_step${STEP})
    add_executable(${TEST_NAME} tests/${TEST_NAME}.cpp)
    target_link_libraries(${TEST_NAME} PRIVATE simulspc_core)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
