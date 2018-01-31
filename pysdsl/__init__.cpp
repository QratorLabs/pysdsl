/*cppimport
<%
cfg['compiler_args'] = ['-std=c++14']
cfg['include_dirs'] = ['sdsl-lite/include']
cfg['libraries'] = ['sdsl', 'divsufsort', 'divsufsort64']
%>
*/

#include <algorithm>
#include <cstdint>
#include <string>

#include "sdsl/vectors.hpp"
#include "sdsl/enc_vector.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

using sdsl::int_vector;


template <class T, typename S=uint64_t>
auto add_class_(py::module& m, const char* name) {
    return py::class_<T>(m, name)
        .def_property_readonly("width", (uint8_t (T::*) (void) const) &T::width)
        .def_property_readonly("data", (const uint64_t* (T::*) (void) const)  &T::data)
        .def("__len__", &T::size, "The number of elements in the int_vector.")
        .def_property_readonly("size", &T::size, "The number of elements in the int_vector.")
        .def_property_readonly_static("max_size",  [](py::object /* self */) {
            return T::max_size();
        }, "Maximum size of the int_vector.")
        .def_property_readonly("size_in_mega_bytes", [](const T& self) { return sdsl::size_in_mega_bytes(self); })
        .def_property_readonly("bit_size", &T::bit_size, "The number of bits in the int_vector.")
        .def("resize", &T::resize, "Resize the int_vector in terms of elements.")
        .def("bit_resize", &T::bit_resize, "Resize the int_vector in terms of bits.")
        .def_property_readonly("capacity", &T::capacity, "Returns the size of the occupied bits of the "
                                                         "int_vector. The capacity of a int_vector is "
                                                         "greater or equal to the bit_size of the "
                                                         "vector: capacity >= bit_size).")
        .def("__getitem__", [](const T &self, size_t position) -> S {
            if (position >= self.size()) {
                throw py::index_error(std::to_string(position));
            }
            return static_cast<S>(self[position]);
        })
        .def("__setitem__", [](T &self, size_t position, S value) {
            if (position >= self.size()) {
                throw py::index_error(std::to_string(position));
            }
            self[position] = value;
        })
        .def(
            "set_to_id",
            [](T &self) {
               sdsl::util::set_to_id(self);
            },
            py::call_guard<py::gil_scoped_release>(),
            "Sets each entry of the numerical vector `v` at position `i` to value `i`"
        )
        .def(
            "set_to_value",
            [](T &self, S value) {
                sdsl::util::set_to_value(self, value);
            },
            py::call_guard<py::gil_scoped_release>(),
            py::arg("k"),
            "Set all entries of int_vector to value k. This method pre-calculates the "
            "content of at most 64 words and then repeatedly inserts these words."
        )
        .def("set_zero_bits", [](T &self) {
           sdsl::util::_set_zero_bits(self);
        }, py::call_guard<py::gil_scoped_release>(), "Sets all bits of the int_vector to 0-bits.")
        .def("set_one_bits", [](T &self) {
           sdsl::util::_set_one_bits(self);
        }, py::call_guard<py::gil_scoped_release>(), "Sets all bits of the int_vector to 1-bits.")
        .def(
            "set_random_bits",
            [](T &self, int seed) {
                 sdsl::util::set_random_bits(self, seed);
            },
            py::call_guard<py::gil_scoped_release>(),
            py::arg_v("seed", 0, "If seed = 0, the time is used to initialize the pseudo random number "
                                "generator, otherwise the seed parameter is used."),
           "Sets all bits of the int_vector to pseudo-random bits."
        )
        .def("__imod__", [](T &self, uint64_t m) {
           sdsl::util::mod(self, m);
           return self;
        }, py::is_operator())
        .def("cnt_one_bits", [](const T &self) {
           return sdsl::util::cnt_one_bits(self);
        }, "Number of set bits in vector")
        .def("cnt_onezero_bits", [](const T &self) {
           return sdsl::util::cnt_onezero_bits(self);
        }, "Number of occurrences of bit pattern `10` in vector")
        .def("cnt_zeroone_bits", [](const T &self) {
           return sdsl::util::cnt_zeroone_bits(self);
        }, "Number of occurrences of bit pattern `01` in vector")
        .def("next_bit", [](const T &self, size_t idx) {
            if (idx >= self.bit_size()) {
                throw py::index_error(std::to_string(idx));
            }
           return sdsl::util::next_bit(self, idx);
        }, py::arg("idx"), "Get the smallest position `i` >= `idx` where a bit is set")
        .def("prev_bit", [](const T &self, size_t idx) {
            if (idx >= self.bit_size()) {
                throw py::index_error(std::to_string(idx));
            }
           return sdsl::util::prev_bit(self, idx);
        }, py::arg("idx"), "Get the largest position `i` <= `idx` where a bit is set")
        .def("__str__", [](const T& self) {return sdsl::util::to_string(self); })
        .def("to_latex", [](const T& self) {return sdsl::util::to_latex_string(self); })



        .def("max", [](const T& self) {
            return std::max_element(self.begin(), self.end());
        }, py::call_guard<py::gil_scoped_release>())
        .def("min", [](const T& self) {
            return std::min_element(self.begin(), self.end());
        }, py::call_guard<py::gil_scoped_release>())
        .def("minmax", [](const T& self) {
            return std::minmax_element(self.begin(), self.end());
        }, py::call_guard<py::gil_scoped_release>())
    ;
}


PYBIND11_MODULE(pysdsl, m) {
    m.doc() = "sdsl-lite bindings for python";

    add_class_<int_vector<0>>(m, "IntVector")
        .def(py::init([](size_t size, uint64_t default_value, uint8_t bit_width) {
            return int_vector<0>(size, default_value, bit_width);
        }), py::arg("size"), py::arg("default_value") = 0, py::arg("bit_width") = 64)
        .def("expand_width", [](int_vector<0> &self, size_t width) {
           sdsl::util::expand_width(self, width);
        }, "Expands the integer width to new_width >= v.width().")
        .def("bit_compress", [](int_vector<0> &self) {
           sdsl::util::bit_compress(self);
        }, "Bit compress the int_vector. Determine the biggest value X and then set the int_width to the "
           "smallest possible so that we still can represent X.")
    ;

    add_class_<int_vector<1>, bool>(m, "BitVector")
        .def(py::init([](size_t size, bool default_value) {
            return int_vector<1>(size, default_value, 1);
        }), py::arg("size"), py::arg("default_value") = false)
        .def("flip", &int_vector<1>::flip, "Flip all bits of bit_vector")
    ;

    add_class_<int_vector<8>, uint8_t>(m, "Int8Vector")
        .def(py::init([](size_t size, uint8_t default_value) {
            return int_vector<8>(size, default_value, 8);
        }), py::arg("size"), py::arg("default_value") = 0)
    ;

    add_class_<int_vector<16>, uint16_t>(m, "Int16Vector")
        .def(py::init([](size_t size, uint16_t default_value) {
            return int_vector<16>(size, default_value, 16);
        }), py::arg("size"), py::arg("default_value") = 0)
    ;

    add_class_<int_vector<32>, uint32_t>(m, "Int32Vector")
        .def(py::init([](size_t size, uint32_t default_value) {
            return int_vector<32>(size, default_value, 32);
        }), py::arg("size"), py::arg("default_value") = 0)
    ;

    add_class_<int_vector<64>, uint64_t>(m, "Int64Vector")
        .def(py::init([](size_t size, uint64_t default_value) {
            return int_vector<64>(size, default_value, 64);
        }), py::arg("size"), py::arg("default_value") = 0)
    ;

    //add_class_<int_vector<128>, uint128_t>(m, "Int128Vector")
    //    .def(py::init([](size_t size, uint128_t default_value) {
    //        return int_vector<128>(size, default_value, 128);
    //    }), py::arg("size"), py::arg("default_value") = 0)
    //;
}
