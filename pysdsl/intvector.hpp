#pragma once

#include <string>
#include <tuple>

#include <pybind11/pybind11.h>

#include <sdsl/util.hpp>
#include <sdsl/vectors.hpp>

#include "calc.hpp"
#include "io.hpp"
#include "sizes.hpp"
#include "docstrings.hpp"


namespace py = pybind11;


template <class T, class KEY, typename S = typename T::value_type>
inline
auto add_int_class(py::module& m, py::dict& dict, KEY key,
                   const char *name, const char *doc = nullptr)
{
    auto cls = py::class_<T>(m, name)
        .def_property_readonly("width", (uint8_t(T::*)(void) const) & T::width)
        .def_property_readonly("data",
                               (const uint64_t *(T::*)(void)const) & T::data)

        .def_property_readonly("bit_size", &T::bit_size,
                               "The number of bits in the int_vector.")

        .def("resize", &T::resize,
             "Resize the int_vector in terms of elements.")
        .def("bit_resize", &T::bit_resize,
             "Resize the int_vector in terms of bits.")
        .def_property_readonly("capacity", &T::capacity, doc_capacity)

        .def(
            "__setitem__",
            [](T &self, size_t position, S value) {
                if (position >= self.size())
                {
                    throw std::out_of_range(std::to_string(position));
                }
                self[position] = value;
            }
        )

        .def("set_to_id",
             [](T &self) { sdsl::util::set_to_id(self); },
             py::call_guard<py::gil_scoped_release>(),
             "Sets each entry of the vector at position `i` to value `i`")
        .def("set_to_value",
             [](T &self, S value) { sdsl::util::set_to_value(self, value); },
             py::arg("k"),
             doc_set_to_value,
             py::call_guard<py::gil_scoped_release>()
        )
        .def("set_zero_bits",
             [](T &self) { sdsl::util::_set_zero_bits(self); },
             "Sets all bits of the int_vector to 0-bits.",
             py::call_guard<py::gil_scoped_release>())
        .def("set_one_bits",
             [](T &self) { sdsl::util::_set_one_bits(self); },
             "Sets all bits of the int_vector to 1-bits.",
             py::call_guard<py::gil_scoped_release>())
        .def(
            "set_random_bits",
            [](T &self, int seed) {
                sdsl::util::set_random_bits(self, seed);
            },
            py::arg_v(
                "seed",
                0,
                "If seed = 0, the time is used to initialize the pseudo "
                "random number generator, otherwise the seed parameter is used."
            ),
            "Sets all bits of the int_vector to pseudo-random bits.",
            py::call_guard<py::gil_scoped_release>()
        )
        .def_static(
            "rnd_positions",
            [](uint8_t log_s, uint64_t mod, uint64_t seed) {
                uint64_t mask;

                auto res = sdsl::util::rnd_positions<T>(log_s, mask, mod, seed);

                return std::make_tuple(res, mask);
            },
            py::arg("log_s"), py::arg("mod") = 0, py::arg("seed") = 0,
            "Create `2**{log_s}` random integers mod `mod` with seed `seed`",
            py::call_guard<py::gil_scoped_release>()
        )
        .def(
            "__imod__",
            [](T &self, uint64_t m) {
                sdsl::util::mod(self, m);
                return self;
            },
            py::is_operator()
        )

        .def("cnt_one_bits",
            [](const T &self) { return sdsl::util::cnt_one_bits(self); },
            "Number of set bits in vector",
            py::call_guard<py::gil_scoped_release>())
        .def("cnt_onezero_bits",
             [](const T &self) { return sdsl::util::cnt_onezero_bits(self); },
             "Number of occurrences of bit pattern `10` in vector",
             py::call_guard<py::gil_scoped_release>())
        .def("cnt_zeroone_bits",
             [](const T &self) { return sdsl::util::cnt_zeroone_bits(self); },
             "Number of occurrences of bit pattern `01` in vector",
             py::call_guard<py::gil_scoped_release>())

        .def(
            "next_bit",
            [](const T &self, size_t idx) {
                if (idx >= self.bit_size())
                {
                    throw std::out_of_range(std::to_string(idx));
                }
                return sdsl::util::next_bit(self, idx);
            },
            py::arg("idx"),
            "Get the smallest position `i` >= `idx` where a bit is set",
            py::call_guard<py::gil_scoped_release>()
        )
        .def(
            "prev_bit",
            [](const T &self, size_t idx) {
                if (idx >= self.bit_size())
                {
                    throw std::out_of_range(std::to_string(idx));
                }
                return sdsl::util::prev_bit(self, idx);
            },
            py::arg("idx"),
            "Get the largest position `i` <= `idx` where a bit is set",
            py::call_guard<py::gil_scoped_release>()
        )
    ;

    add_sizes(cls);
    add_description(cls);
    add_serialization(cls);
    add_to_string(cls);

    add_std_algo<T, S>(cls);

    if (doc) cls.doc() = doc;

    dict.attr("__setitem__")(key, cls);

    return cls;
}


inline
auto add_int_vectors(py::module& m)
{
    py::dict int_vectors_dict;

    m.attr("int_vector") = int_vectors_dict;

    return std::make_tuple(
        add_int_class<sdsl::int_vector<0>>(
            m, int_vectors_dict, "dynamic", "IntVector", doc_int_vector
        )
            .def(
                py::init([](size_t size,
                            uint64_t default_value,
                            uint8_t bit_width) {
                    return sdsl::int_vector<0>(size, default_value, bit_width);
                }),
                py::arg("size") = 0,
                py::arg("default_value") = 0,
                py::arg("bit_width") = 64,
                py::call_guard<py::gil_scoped_release>())
            .def(
                "expand_width",
                [](sdsl::int_vector<0> &self, size_t width) {
                    sdsl::util::expand_width(self, width);
                },
                "Expands the integer width to new_width >= v.width().",
                py::call_guard<py::gil_scoped_release>()
            )
            .def("bit_compress",
                [](sdsl::int_vector<0> &self) {
                    sdsl::util::bit_compress(self);
                },
                doc_bit_compress,
                py::call_guard<py::gil_scoped_release>()),

        add_int_class<sdsl::int_vector<1>, bool>(
            m, int_vectors_dict, 1ul , "BitVector"
        )
            .def(py::init([](size_t size, bool default_value) {
                    return sdsl::int_vector<1>(size, default_value, 1);
                }), py::arg("size") = 0, py::arg("default_value") = false)
            .def("flip", &sdsl::int_vector<1>::flip,
                 "Flip all bits of bit_vector",
                 py::call_guard<py::gil_scoped_release>()),

        add_int_class<sdsl::int_vector<4>, uint8_t>(
            m, int_vectors_dict, 4, "Int4Vector"
        )
            .def(py::init([](size_t size, uint8_t default_value) {
                    return sdsl::int_vector<4>(size, default_value, 4);
                }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<sdsl::int_vector<8>, uint8_t>(
            m, int_vectors_dict, 8, "Int8Vector"
        )
            .def(py::init([](size_t size, uint8_t default_value) {
                    return sdsl::int_vector<8>(size, default_value, 8);
                }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<sdsl::int_vector<16>, uint16_t>(
            m, int_vectors_dict, 16, "Int16Vector"
        )
            .def(py::init([](size_t size, uint16_t default_value) {
                     return sdsl::int_vector<16>(size, default_value, 16);
                 }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<sdsl::int_vector<24>, uint32_t>(
            m, int_vectors_dict, 24, "Int24Vector"
        )
            .def(py::init([](size_t size, uint32_t default_value) {
                     return sdsl::int_vector<24>(size, default_value, 24);
                 }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<sdsl::int_vector<32>, uint32_t>(
            m, int_vectors_dict, 32, "Int32Vector"
        )
            .def(py::init([](size_t size, uint32_t default_value) {
                    return sdsl::int_vector<32>(size, default_value, 32);
                }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<sdsl::int_vector<64>, uint64_t>(
            m, int_vectors_dict, 64, "Int64Vector"
        )
            .def(py::init([](size_t size, uint64_t default_value) {
                    return sdsl::int_vector<64>(size, default_value, 64);
                }), py::arg("size") = 0, py::arg("default_value") = 0)
    );
}