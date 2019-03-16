#pragma once

#include <string>
#include <tuple>

#include <pybind11/pybind11.h>

#include <sdsl/util.hpp>
#include <sdsl/sorted_int_stack.hpp>

#include "calc.hpp"
#include "io.hpp"
#include "operations/iteration.hpp"
#include "operations/sizes.hpp"
#include "docstrings.hpp"


namespace py = pybind11;


inline auto add_sorted_int_stack(py::module& m)
{
    py::dict sorted_int_stack_dict;
    m.attr("sorted_int_stack") = sorted_int_stack_dict;

    using T = sdsl::sorted_int_stack;

    auto cls = py::class_<T>(m, "SortedIntStack")
        .def("empty", &T::empty)
        .def("size", &T::size)
        .def("top", &T::top)
        .def("pop", &T::pop)
        .def("push", &T::push)
        .def(py::init([](T::size_type x) {
            return T(x);
        }));

    add_sizes(cls);
    // add_description(cls);
    // add_serialization(cls);

    return cls;
}
