// bindings.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "SolveResult.h"
#include "Component.h"
#include "wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(eos_cpp, m) {
    py::class_<Component>(m, "Component")
        .def(py::init<std::string, std::string, double, double, double, double, double,
            double, double, double, double>())
        .def_readwrite("name", &Component::name)
        .def_readwrite("formula", &Component::formula)
        .def_readwrite("Mw", &Component::Mw)
        .def_readwrite("Tc", &Component::Tc)
        .def_readwrite("Pc", &Component::Pc)
        .def_readwrite("omega", &Component::omega)
        .def_readwrite("fraction", &Component::fraction)
        .def_readwrite("CpA", &Component::CpA)
        .def_readwrite("CpB", &Component::CpB)
        .def_readwrite("CpC", &Component::CpC)
        .def_readwrite("CpD", &Component::CpD);

    py::class_<SolveResult>(m, "SolveResult")
        .def_readwrite("T", &SolveResult::T)
        .def_readwrite("P", &SolveResult::P)
        .def_readwrite("V", &SolveResult::V)
        .def_readwrite("iterations", &SolveResult::iterations);

    m.def("mainFromPython", &mainFromPython, "Solve EOS problem from Python");
}
