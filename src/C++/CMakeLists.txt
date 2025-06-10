cmake_minimum_required(VERSION 3.14)
project(RachfordRice)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Uključi pybind11 (kao subdirektorij ili git submodule)
add_subdirectory(pybind11)
find_package(OpenMP REQUIRED)

# Svi source fileovi
set(SOURCES
    bindings.cpp
    wrapper.cpp
    Component.cpp
    ComponentData.cpp
    EOSBase.cpp
    PengRobinson.cpp
    RachfordRice.cpp
    SolveResult.cpp
    Calculations.cpp
)

# Kreiraj pybind11 modul koji Python može uvesti
pybind11_add_module(eos_cpp ${SOURCES})
target_link_libraries(eos_cpp PRIVATE OpenMP::OpenMP_CXX)

# Ako koristiš MSVC: omogući OpenMP
if (MSVC)
    target_compile_options(eos_cpp PRIVATE /openmp)
endif()