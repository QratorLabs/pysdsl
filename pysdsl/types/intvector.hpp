#pragma once

#include <string>
#include <tuple>

#include <pybind11/pybind11.h>

#include <sdsl/util.hpp>
#include <sdsl/vectors.hpp>

#include "calc.hpp"
#include "io.hpp"
#include "operations/iteration.hpp"
#include "operations/sizes.hpp"
#include "docstrings.hpp"


namespace py = pybind11;


constexpr char sym_for_width(unsigned int width) {
    switch (width) {
        case 8:
            return 'B';
        case 16:
            return 'H';
        case 32:
            return 'I';
        case 64:
            return 'Q';
        default: __builtin_unreachable();
    }
}


// checks whether width is a power of 2 (width & (width - 1) == 0)
//                              and this power is between 8 and 64
// without dummy redefinition error
template <class T,
          unsigned int width = static_cast<unsigned int>(T::fixed_int_width),
          std::enable_if_t<!(width & (width - 1)) && (width & (128u - 8u))>* dummy = nullptr>
inline auto add_int_init(py::module& m, const char* name)
{
    return py::class_<T>(m, name, py::buffer_protocol())
        .def_buffer([] (T& self) {
            return py::buffer_info(
                reinterpret_cast<void*>(self.data()),
                width / 8,
                std::string(1, sym_for_width(width)),
                1,
                { detail::size(self) },
                { width / 8 }
            ); });
}

template <class T,
          unsigned int width = static_cast<unsigned int>(T::fixed_int_width),
          std::enable_if_t<(width & (width - 1)) || !(width & (128u - 8u))>* dummy = nullptr>
inline auto add_int_init(py::module& m, const char* name)
{
    return py::class_<T>(m, name);
}
        

template <class T, typename S = typename T::value_type, typename KEY_T>
inline auto add_int_class(py::module& m, py::dict& dict, KEY_T key,
                          const char *name, const char *doc = nullptr)
{
    auto cls = add_int_init<T>(m, name)
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
                if (position >= self.size()) {
                    throw std::out_of_range(std::to_string(position)); }
                self[position] = value; })

        .def("set_to_id",
             [](T &self) { sdsl::util::set_to_id(self); },
             py::call_guard<py::gil_scoped_release>(),
             "Sets each entry of the vector at position `i` to value `i`")
        .def("set_to_value",
             [](T &self, S value) { sdsl::util::set_to_value(self, value); },
             py::arg("k"),
             doc_set_to_value,
             py::call_guard<py::gil_scoped_release>())
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
                sdsl::util::set_random_bits(self, seed); },
            py::arg_v(
                "seed",
                0,
                "If seed = 0, the time is used to initialize the pseudo "
                "random number generator, otherwise the seed parameter is used."
            ),
            "Sets all bits of the int_vector to pseudo-random bits.",
            py::call_guard<py::gil_scoped_release>())
        .def_static(
            "rnd_positions",
            [](uint8_t log_s, uint64_t mod, uint64_t seed) {
                uint64_t mask;

                auto res = sdsl::util::rnd_positions<T>(log_s, mask, mod, seed);

                return std::make_tuple(res, mask); },
            py::arg("log_s"), py::arg("mod") = 0, py::arg("seed") = 0,
            "Create `2**{log_s}` random integers mod `mod` with seed `seed`",
            py::call_guard<py::gil_scoped_release>())
        .def(
            "__imod__",
            [](T &self, uint64_t m) {
                sdsl::util::mod(self, m);
                return self; },
            py::is_operator())

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
                if (idx >= self.bit_size()) {
                    throw std::out_of_range(std::to_string(idx)); }
                return sdsl::util::next_bit(self, idx); },
            py::arg("idx"),
            "Get the smallest position `i` >= `idx` where a bit is set",
            py::call_guard<py::gil_scoped_release>())
        .def(
            "prev_bit",
            [](const T &self, size_t idx) {
                if (idx >= self.bit_size()) {
                    throw std::out_of_range(std::to_string(idx)); }
                return sdsl::util::prev_bit(self, idx); },
            py::arg("idx"),
            "Get the largest position `i` <= `idx` where a bit is set",
            py::call_guard<py::gil_scoped_release>());

    add_sizes(cls);
    add_description(cls);
    add_serialization(cls);
    add_to_string<T, S>(cls);

    add_read_access<T, S>(cls);
    add_std_algo<T, S>(cls);

    if (doc) cls.doc() = doc;

    dict.attr("__setitem__")(key, cls);

    return cls;
}


struct AddIntVectorFunctor {
    py::module& m;
    py::dict& int_vectors_dict;

    constexpr AddIntVectorFunctor(py::module& m, py::dict& int_vectors_dict) noexcept
        : m(m), int_vectors_dict(int_vectors_dict) {}

    template <size_t N>
    auto operator()(std::integral_constant<size_t, N> t) {
        using return_type = sdsl::int_vector<N>;
        std::string name = "Int" + std::to_string(N) + "Vector";
        return add_int_class<return_type, typename return_type::value_type>(
                    m, int_vectors_dict, N, name.c_str())
                .def(py::init(
                    [](size_t size, typename return_type::value_type default_value) {
                        return return_type(size, default_value, N); }),
                    py::arg("size") = 0, py::arg("default_value") = 0);
    }

    auto operator()(std::integral_constant<size_t, 0> t) {
        return add_int_class<sdsl::int_vector<0>>(
                    m, int_vectors_dict, "dynamic", "IntVector", doc_int_vector)
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
                        sdsl::util::expand_width(self, width); },
                    "Expands the integer width to new_width >= v.width().",
                    py::call_guard<py::gil_scoped_release>())
                .def("bit_compress",
                    [](sdsl::int_vector<0> &self) {
                        sdsl::util::bit_compress(self); },
                    doc_bit_compress,
                    py::call_guard<py::gil_scoped_release>());
    }

    auto operator()(std::integral_constant<size_t, 1> t) {
        return add_int_class<sdsl::int_vector<1>, bool>(
                m, int_vectors_dict, 1ul , "BitVector")
            .def(py::init(
                [](size_t size, bool default_value) {
                    return sdsl::int_vector<1>(size, default_value, 1); }),
                py::arg("size") = 0, py::arg("default_value") = false)
            .def("flip", &sdsl::int_vector<1>::flip,
                 "Flip all bits of bit_vector",
                 py::call_guard<py::gil_scoped_release>());
    }
};


template <typename... Ts>
struct SubsetFunctor {
    const std::tuple<Ts...>& tpl;

    // it doesn't work in C++14
    constexpr SubsetFunctor(std::tuple<Ts...>& tpl) noexcept
        : tpl(tpl) {}


    template <size_t N>
    const auto& operator()(std::integral_constant<size_t, N> t) const {
        return std::get<py::class_<sdsl::int_vector<N>>>(tpl);
    }
};


inline auto add_int_vectors(py::module& m)
{
    py::dict int_vectors_dict;

    m.attr("int_vector") = int_vectors_dict;

    using params = std::tuple<
        std::integral_constant<size_t, 1>,
        std::integral_constant<size_t, 4>,
        std::integral_constant<size_t, 8>,
        std::integral_constant<size_t, 16>,
        std::integral_constant<size_t, 24>,
        std::integral_constant<size_t, 32>,
        std::integral_constant<size_t, 48>,
        std::integral_constant<size_t, 64>>;

    using for_enc_vectors = std::tuple<
        std::integral_constant<size_t, 1>,
        std::integral_constant<size_t, 4>,
        std::integral_constant<size_t, 8>,
        std::integral_constant<size_t, 64>>;

    auto iv = for_each_in_tuple(params(), AddIntVectorFunctor(m, int_vectors_dict));

    // it should be a tuple of references
    // mb use of std::forward_as_tuple needed in for_each_in_tuple instead of std::make_tuple
    auto iv_as_ev_params = for_each_in_tuple(for_enc_vectors(), SubsetFunctor(iv));

    return std::make_tuple(iv, iv_as_ev_params);
}
