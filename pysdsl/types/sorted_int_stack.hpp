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
    using Stack = sdsl::sorted_int_stack;

    auto cls = py::class_<Stack>(m, "SortedIntStack")
        .def("empty", &Stack::empty, "Checks whether the stack is empty.")
        .def("top", [](const Stack& self) {
            if (self.size() > 0u)
                return self.top();
            throw py::index_error("top from empty stack");
        }, "Returns the topmost element of the stack.")
        .def("pop", [](Stack& self) {
            if (self.size() == 0u)
                throw py::index_error("pop from empty stack");
            auto ans = self.top(); 
            self.pop(); 
            return ans;
        }, "Removes the topmost element from the stack and returns its copy. Not thread safe.")
        .def("push", [](Stack& self, const Stack::size_type& x) {
            if (self.empty() || self.top() < x)
                self.push(x);
            else
                throw py::value_error("elements have to be pushed in strictly increasing order");
        }, "Adds new element to the top of the stack."
           "(n.b. it has to be not less than the stored ones). Not thread safe.")
        .def(py::init([](Stack::size_type x) {
            return Stack(x);
        }), py::arg("max_value"),
        "Creates a stack which can store integers not greater than max_value.");


    cls.doc() = doc_sorted_int_stack;

    add_sizes(cls);
    add_serialization(cls, 0);
    add_description(cls);

    return std::make_tuple(cls);
}
