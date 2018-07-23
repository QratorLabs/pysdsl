#pragma once

#include <string>
#include <tuple>
#include <utility>

#include <pybind11/pybind11.h>

#include <sdsl/suffix_arrays.hpp>

#include "operations/sizes.hpp"
#include "operations/iteration.hpp"
#include "docstrings.hpp"
#include "io.hpp"
#include "calc.hpp"

namespace py = pybind11;


template <class T>
inline auto add_csa_member_type(py::module& m,
                                std::string&& name,
                                const std::string& parent_name)
{
    try {
        auto cls = py::class_<T>(
            m, ("_" + name + "ofSuffixArray" + parent_name).c_str());
        add_read_access(cls);
        add_to_string(cls);
    } catch (std::runtime_error& /* ignore */) {}
}


template <class T>
inline
auto add_csa_class(py::module& m, std::string&& name, const char* doc = nullptr)
{
    auto cls = py::class_<T>(m, ("SuffixArray" + name).c_str());

    add_csa_member_type<typename T::isa_type>(m, "ISA", name);
    add_csa_member_type<typename T::bwt_type>(m, "BWT", name);
    add_csa_member_type<typename T::lf_type>(m, "LF", name);
    add_csa_member_type<typename T::psi_type>(m, "PSI", name);
    add_csa_member_type<typename T::text_type>(m, "Text", name);
    add_csa_member_type<typename T::first_row_type>(m, "FirstRow", name);

    try {
        using char2comp_type = typename T::alphabet_type::char2comp_type;
        auto cls_char2comp = py::class_<char2comp_type>(
            m, ("_Char2CompOf" + name).c_str());
        cls_char2comp.def(
            "__getitem__",
            [] (const char2comp_type& self, uint64_t c) { return self[c]; }
        );
    } catch (std::runtime_error& /* ignore */) {}
    try {
        using comp2char_type = typename T::alphabet_type::comp2char_type;
        auto cls_comp2char = py::class_<comp2char_type>(
            m, ("_Comp2CharOf" + name).c_str());
        cls_comp2char.def(
            "__getitem__",
            [] (const comp2char_type& self, uint64_t c) { return self[c]; }
        );
    } catch (std::runtime_error& /* ignore */) {}

    cls.def_property_readonly("isa", [] (const T& self ) { return &self.isa; });
    cls.def_property_readonly("bwt", [] (const T& self ) { return &self.bwt; });
    cls.def_property_readonly("lf", [] (const T& self ) { return &self.lf; });
    cls.def_property_readonly("psi", [] (const T& self ) { return &self.psi; });
    cls.def_property_readonly("text", [] (const T& self ) {
        return &self.text; });
    cls.def_property_readonly("L", [] (const T& self ) { return &self.L; });
    cls.def_property_readonly("F", [] (const T& self ) { return &self.F; });
    cls.def_property_readonly("C", [] (const T& self ) { return &self.C; });
    cls.def_property_readonly("char2comp", [] (const T& self ) {
        return &self.char2comp; });
    cls.def_property_readonly("comp2char", [] (const T& self ) {
        return &self.comp2char; });
    cls.def_property_readonly("sigma", [] (const T& self ) {
        return self.sigma; });

    cls.def(
        "extract",
        [] (const T& self, typename T::size_type begin,
            typename T::size_type end) {
            if (end >= detail::size(self)) {
                throw std::out_of_range(std::to_string(end)); }
            if (begin >= end) {
                throw std::invalid_argument("begin should be less than end"); }
            return sdsl::extract(self, begin, end); },
        py::arg("begin"),
        py::arg("end"),
        "Reconstructs the subarray T[begin:end] of the original array T\n"
        "\n\tbegin: Position of the first character which should be extracted "
        "(inclusive)"
        "\n\tend: Position of the last character which should be extracted "
        "(inclusive)\n\n"
        "Time complexity: Order{(end - begin+1) * t_{Psi} + t_{ISA} }",
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "count",
        [] (const T& self, const typename T::string_type& pattern) {
            return sdsl::count(self, pattern); },
        py::arg("pattern"),
        "Counts the number of occurrences of a pattern in a CSA",
        py::call_guard<py::gil_scoped_release>());
    cls.def(
        "locate",
        [] (const T& self, const typename T::string_type& pattern) {
            return sdsl::locate(self, pattern); },
        py::arg("pattern"),
        "Calculates all occurrences of a pattern in a CSA\n"
        "Time complexity:"
        "Order{ t_{backward_search} + z * t_{SA} },\n"
        "where `z` is the number of occurrences of pattern in the CSA",
        py::call_guard<py::gil_scoped_release>());
    cls.def(py::init(
        [] (const typename T::string_type& data)
        {
            T self;
            sdsl::construct_im(self, data,
                               sizeof(typename T::string_type::value_type));
            return self;
        }
    ));

    add_sizes(cls);
    add_description(cls);
    add_serialization(cls);
    add_to_string(cls);

    add_read_access<T>(cls);
    add_std_algo<T>(cls);

    if (doc) cls.doc() = doc;

    m.attr("suffix_array").attr("__setitem__")(name, cls);

    return cls;
}


inline auto add_csa(py::module& m)
{
    m.attr("suffix_array") = py::dict();

    auto csa_classes = std::make_tuple(
        add_csa_class<sdsl::csa_bitcompressed<>>(m, "Bitcompressed", doc_csa),
        add_csa_class<sdsl::csa_sada<>>(m, "Sadakane", doc_sada),
        add_csa_class<sdsl::csa_sada_int<>>(m, "SadakaneInt", doc_sada),
        add_csa_class<sdsl::csa_wt<>>(m, "WaveletTree", doc_csa_wt),
        add_csa_class<sdsl::csa_wt_int<>>(m, "WaveletTreeInt", doc_csa_wt));

    return csa_classes;
}
