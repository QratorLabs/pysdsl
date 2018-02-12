#pragma once

#include <pybind11/pybind11.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/io.hpp>


namespace py = pybind11;


template <class T>
auto add_max_size(py::class_<T>& cls)
{
    cls.def_property_readonly_static(
        "max_size",
        [](py::object /* self */) { return T::max_size(); },
        "Maximum size of the int_vector."
    );
    return cls;
}

template <class BV, int ML>
auto add_max_size(py::class_<sdsl::dac_vector_dp<BV, ML>>& cls)
{
    return cls;
}

template <uint32_t B>
auto add_max_size(
    py::class_<sdsl::bit_vector_il<B>>& cls
)
{
    cls.def_property_readonly_static(
        "max_size",
        [](py::object /* self */) { return sdsl::bit_vector::max_size(); },
        "Maximum size of the int_vector."
    );
    return cls;
}

template <class Th, class Ts1, class Ts0>
auto add_max_size(
    py::class_<sdsl::sd_vector<Th, Ts1, Ts0>>& cls
) { return cls; }

template <uint16_t B, class R, uint16_t K>
auto add_max_size(
    py::class_<sdsl::rrr_vector<B, R, K>>& cls
) { return cls; }

template <uint32_t B>
auto add_max_size(py::class_<sdsl::hyb_vector<B>>& cls) { return cls; }


template <class T>
auto add_sizes(py::class_<T>& cls)
{
    cls.def("__len__", &T::size, "The number of elements in the int_vector.");
    cls.def_property_readonly("size", &T::size,
                              "The number of elements in the int_vector.");
    cls.def_property_readonly(
        "size_in_mega_bytes",
        [](const T &self) { return sdsl::size_in_mega_bytes(self); }
    );
    cls.def_property_readonly(
        "size_in_bytes",
        [](const T &self) { return sdsl::size_in_bytes(self); }
    );

    return add_max_size(cls);
}
