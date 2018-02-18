#pragma once

#include <algorithm>

#include <pybind11/pybind11.h>

namespace py = pybind11;


template <class Sequence>
auto add_std_algo(py::class_<Sequence>& cls)
{
    typedef typename Sequence::value_type value_type;

    cls.def(
        "__iter__",
        [](const Sequence &s) { return py::make_iterator(s.begin(), s.end()); },
        py::keep_alive<0, 1>()
    );
    cls.def(
        "__contains__",
        [](const Sequence &self, typename Sequence::value_type element) {
            return std::find(self.begin(), self.end(), element) != self.end();
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "__getitem__",
        [](const Sequence &self, size_t position) {
            if (position >= self.size())
            {
                throw py::index_error(std::to_string(position));
            }
            return self[position];
        }
    );
    cls.def(
        "__getitem__",
        [](const Sequence &self, int64_t position) {
            auto abs_position = std::abs(position);
            if (position >= 0)
            {
                throw std::exception();
            }
            if (abs_position > self.size())
            {
                throw py::index_error(std::to_string(position));
            }
            return self[self.size() - abs_position];
        }
    );
    cls.def(
        "max",
        [](const Sequence &self) {
            return *std::max_element(self.begin(), self.end());
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "argmax",
        [](const Sequence &self) {
            return std::distance(self.begin(),
                                 std::max_element(self.begin(), self.end()));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "min",
        [](const Sequence &self) {
            return *std::min_element(self.begin(), self.end());
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "argmin",
        [](const Sequence &self) {
            return std::distance(self.begin(),
                                 std::min_element(self.begin(), self.end()));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "minmax",
        [](const Sequence &self) {
            auto result = std::minmax_element(self.begin(), self.end());
            return std::make_pair(*std::get<0>(result), *std::get<1>(result));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "sum",
        [](const Sequence &self) {
            return std::accumulate(self.begin(), self.end(),
                                   uint64_t(0));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "all",
        [](const Sequence &self) {
            return std::all_of(self.begin(), self.end(),
                               [] (const value_type value) -> bool
                               { return value; });
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "any",
        [](const Sequence &self) {
            return std::any_of(self.begin(), self.end(),
                               [] (const value_type value) -> bool
                               { return value; });
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "none",
        [](const Sequence &self) {
            return std::none_of(self.begin(), self.end(),
                                [] (const value_type value) -> bool
                                { return value; });
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "is_sorted",
        [](const Sequence &self) {
            return std::is_sorted(self.begin(), self.end());
        },
        py::call_guard<py::gil_scoped_release>()
    );

    return cls;
}
