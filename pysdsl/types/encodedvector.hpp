#pragma once

#include <string>
#include <tuple>
#include <type_traits>

#include <pybind11/pybind11.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

#include "calc.hpp"
#include "docstrings.hpp"
#include "io.hpp"
#include "operations/iteration.hpp"
#include "operations/sizes.hpp"
#include "util/tupletricks.hpp"


namespace py = pybind11;


namespace {
template <std::size_t I> using dens = std::integral_constant<uint32_t, I>;

template <std::size_t W> using width = std::integral_constant<uint8_t, W>;
}  // namespace


auto constexpr coders = std::make_tuple(
    std::make_tuple("EliasDelta", dens<128>{}, width<0>{},
        (sdsl::coder::elias_delta*) nullptr),
    std::make_tuple("EliasGamma", dens<128>{}, width<0>{},
        (sdsl::coder::elias_gamma*) nullptr),
    std::make_tuple("Fibonacci", dens<128>{}, width<0>{},
        (sdsl::coder::fibonacci*) nullptr),
    std::make_tuple("Comma2", dens<128>{}, width<0>{},
        (sdsl::coder::comma<2>*) nullptr),
    std::make_tuple("Comma4", dens<128>{}, width<0>{},
        (sdsl::coder::comma<4>*) nullptr));


class add_enc_coders_functor
{
public:
    constexpr add_enc_coders_functor(py::module& m): m(m) {}

    template <typename Coder, uint32_t t_dens, uint8_t t_width>
    inline
    decltype(auto) operator()(const std::tuple<const char*, dens<t_dens>,
                                               width<t_width>, Coder*> &t)
    {
        using enc = sdsl::enc_vector<Coder, t_dens, t_width>;

        auto cls = py::class_<enc>(
                m,
                (std::string("EncVector") + std::get<0>(t)).c_str()
            ).def(py::init());

        add_sizes(cls);
        add_description(cls);
        add_serialization(cls);
        add_to_string(cls);

        add_read_access<enc>(cls);
        add_std_algo<enc>(cls);

        cls.doc() = "A vector `v` is stored more space-efficiently by "
                    "self-delimiting coding the deltas v[i+1]-v[i] (v[-1]:=0).";

        cls.def_property_readonly("sample_dens", &enc::get_sample_dens)
        .def(
            "sample",
            [] (const enc& self, typename enc::size_type i) {
                if (i >= self.size() / self.get_sample_dens()) {
                    throw std::out_of_range(std::to_string(i)); }
                return self.sample(i); },
            "Returns the i-th sample of the compressed vector"
            "i: The index of the sample. 0 <= i < size()/get_sample_dens()",
            py::call_guard<py::gil_scoped_release>())
        .def(
            "samples",
            [] (const enc& self) {
                const typename enc::size_type size = self.size() /
                                                        self.get_sample_dens();

                sdsl::int_vector<t_width> samples(size);
                for (std::size_t i = 0; i <= size; i++) {
                    samples[i] = self.sample(i); }
                return samples; },
            py::call_guard<py::gil_scoped_release>());

        m.attr("enc_vector").attr("__setitem__")(std::get<0>(t), cls);
        m.attr("all_compressed_integer_vectors").attr("append")(cls);

        return cls;
    }

private:
    py::module& m;
};


class add_vlc_coders_functor
{
public:
    constexpr add_vlc_coders_functor(py::module& m): m(m) {}

    template <typename Coder, uint32_t t_dens, uint8_t t_width>
    inline
    decltype(auto) operator()(const std::tuple<const char*, dens<t_dens>,
                                               width<t_width>, Coder*> &t)
    {
        using vlc = sdsl::vlc_vector<Coder, t_dens, t_width>;

        auto cls = py::class_<vlc>(m, (
                std::string("VariableLengthCodesVector") + std::get<0>(t)
        ).c_str()).def(py::init());

        add_sizes(cls);
        add_description(cls);
        add_serialization(cls);
        add_to_string(cls);

        add_read_access<vlc>(cls);
        add_std_algo<vlc>(cls);

        cls.doc() = "A vector which stores the values with "
                    "variable length codes.";

        cls.def_property_readonly("sample_dens", &vlc::get_sample_dens);

        m.attr("variable_length_codes_vector").attr(
                "__setitem__")(std::get<0>(t), cls);
        m.attr("all_compressed_integer_vectors").attr("append")(cls);

        return cls;
    }

private:
    py::module& m;
};


inline std::string key_to_string(const char* key) { return std::string(key); }

template <class KEY_T>
inline std::string key_to_string(KEY_T key) { return std::to_string(key); }


template <class Sequence, typename KEY_T>
inline
auto add_dac_vector(py::module& m, KEY_T key,
                    const char* doc = nullptr)
{
    auto name = "DirectAccessibleCodesVector" + key_to_string(key);

    auto cls = py::class_<Sequence>(m, name.c_str()).def(py::init());

    add_sizes(cls);
    add_description(cls);
    add_serialization(cls);
    add_to_string(cls);

    add_read_access<Sequence>(cls);
    add_std_algo<Sequence>(cls);

    if (doc) {
        cls.doc() = doc; }

    cls.def_property_readonly("levels", &Sequence::levels);

    m.attr("direct_accessible_codes_vector").attr("__setitem__")(key, cls);
    m.attr("all_compressed_integer_vectors").attr("append")(cls);

    return cls;
}


auto add_encoded_vectors(py::module& m)
{
    m.attr("enc_vector") = py::dict();
    m.attr("variable_length_codes_vector") = py::dict();
    m.attr("direct_accessible_codes_vector") = py::dict();
    m.attr("all_compressed_integer_vectors") = py::list();

    auto enc_classes = for_each_in_tuple(coders, add_enc_coders_functor(m));
    auto vlc_classes = for_each_in_tuple(coders, add_vlc_coders_functor(m));
    auto dac_classes = std::make_tuple(
        add_dac_vector<sdsl::dac_vector<4>>(m, 4, doc_dac_vector),
        add_dac_vector<sdsl::dac_vector<8>>(m, 8, doc_dac_vector),
        add_dac_vector<sdsl::dac_vector<16>>(m, 16, doc_dac_vector),
        add_dac_vector<sdsl::dac_vector<63>>(m, 63, doc_dac_vector),
        add_dac_vector<sdsl::dac_vector_dp<>>(m, "DP", doc_dac_vector_dp)
            .def("cost", &sdsl::dac_vector_dp<>::cost,
                 py::arg("n"), py::arg("m")),
        add_dac_vector<
            sdsl::dac_vector_dp<
                sdsl::rrr_vector<>>>(m, "DPRRR", doc_dac_vector_dp)
            .def("cost", &sdsl::dac_vector_dp<sdsl::rrr_vector<>>::cost,
                 py::arg("n"), py::arg("m"))
    );

    m.attr("DACVector") = m.attr("DirectAccessibleCodesVector4");
    m.attr("DirectAccessibleCodesVector") = m.attr(
        "DirectAccessibleCodesVector4");

    return std::tuple_cat(enc_classes, vlc_classes, dac_classes);
}
