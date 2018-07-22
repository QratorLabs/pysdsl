#pragma once

#include <string>
#include <tuple>

#include <pybind11/pybind11.h>

#include <sdsl/suffix_arrays.hpp>

#include "operations/sizes.hpp"
#include "operations/iteration.hpp"
#include "docstrings.hpp"
#include "io.hpp"
#include "calc.hpp"

namespace py = pybind11;


template <class T>
auto add_csa_class(py::module& m, std::string&& name, const char* doc = nullptr)
{
    auto cls = py::class_<T>(m, ("SuffixArray" + name).c_str());

    cls.def_property_readonly("isa", [] (const T& self ) { return self.isa; });
    cls.def_property_readonly("bwt", [] (const T& self ) { return self.bwt; });
    cls.def_property_readonly("lf", [] (const T& self ) { return self.lf; });
    cls.def_property_readonly("psi", [] (const T& self ) { return self.psi; });
    cls.def_property_readonly("text", [] (const T& self ) { return self.text; });
    cls.def_property_readonly("L", [] (const T& self ) { return self.L; });
    cls.def_property_readonly("F", [] (const T& self ) { return self.F; });
    cls.def_property_readonly("C", [] (const T& self ) { return self.C; });
    cls.def_property_readonly("char2comp", [] (const T& self ) { return self.char2comp; });
    cls.def_property_readonly("comp2char", [] (const T& self ) { return self.comp2char; });
    cls.def_property_readonly("sigma", [] (const T& self ) { return self.sigma; });

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
