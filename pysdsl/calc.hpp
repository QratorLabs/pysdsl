#pragma once

#include <algorithm>
#include <tuple>

#include <pybind11/pybind11.h>

#include "indexiterator.hpp"
#include "operations/sizes.hpp"


namespace py = pybind11;


using std::cbegin;
using std::cend;


template <class Sequence, typename T = typename Sequence::value_type>
inline
auto add_std_algo(py::class_<Sequence>& cls)
{
    typedef typename Sequence::value_type value_type;

    cls.def(
        "__iter__",
        [](const Sequence &s) { return py::make_iterator(cbegin(s),
                                                         cend(s)); },
        py::keep_alive<0, 1>()
    );
    cls.def(
        "__contains__",
        [](const Sequence &self, typename Sequence::value_type element) {
            return std::find(cbegin(self),
                             cend(self), element) != cend(self);
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "__getitem__",
        [](const Sequence &self, size_t position) -> T {
            if (position >= detail::size(self))
            {
                throw std::out_of_range(std::to_string(position));
            }
            return self[position];
        }
    );
    cls.def(
        "__getitem__",
        [](const Sequence &self, int64_t position) -> T {
            auto abs_position = std::abs(position);
            if (position >= 0)
            {
                throw std::exception();
            }
            if (abs_position > detail::size(self))
            {
                throw std::out_of_range(std::to_string(position));
            }
            return self[detail::size(self) - abs_position];
        }
    );
    cls.def(
        "max",
        [](const Sequence &self) -> T {
            return *std::max_element(cbegin(self), cend(self));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "argmax",
        [](const Sequence &self) {
            return std::distance(cbegin(self),
                                 std::max_element(cbegin(self),
                                                  cend(self)));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "min",
        [](const Sequence &self) -> T {
            return *std::min_element(cbegin(self), cend(self));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "argmin",
        [](const Sequence &self) {
            return std::distance(cbegin(self),
                                 std::min_element(cbegin(self),
                                                  cend(self)));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "minmax",
        [](const Sequence &self) -> std::pair<T, T> {
            auto result = std::minmax_element(cbegin(self),
                                              cend(self));
            return std::make_pair(*std::get<0>(result), *std::get<1>(result));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "sum",
        [](const Sequence &self) {
            return std::accumulate(cbegin(self), cend(self),
                                   uint64_t(0));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "all",
        [](const Sequence &self) {
            return std::all_of(cbegin(self), cend(self),
                               [] (const value_type value) -> bool
                               { return value; });
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "any",
        [](const Sequence &self) {
            return std::any_of(cbegin(self), cend(self),
                               [] (const value_type value) -> bool
                               { return value; });
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "none",
        [](const Sequence &self) {
            return std::none_of(cbegin(self), cend(self),
                                [] (const value_type value) -> bool
                                { return value; });
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "is_sorted",
        [](const Sequence &self) {
            return std::is_sorted(cbegin(self), cend(self));
        },
        py::call_guard<py::gil_scoped_release>()
    );

    return cls;
}
