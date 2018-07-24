#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

#include <pybind11/pybind11.h>

#include "calc.hpp"
#include "docstrings.hpp"
#include "io.hpp"
#include "supports.hpp"
#include "operations/sizes.hpp"
#include "operations/iteration.hpp"


namespace py = pybind11;


template <class T>
inline
auto add_bitvector_class(py::module &m, const std::string&& name,
                         const char* doc = nullptr,
                         const char* doc_rank = nullptr,
                         const char* doc_select = nullptr)
{
    auto cls = py::class_<T>(m, name.c_str()).def(py::init());

    add_sizes(cls);
    add_description(cls);
    add_serialization(cls);
    add_to_string(cls);

    add_read_access<T, bool>(cls);
    add_std_algo<T, bool>(cls);

    if (doc) cls.doc() = doc;

    cls.def(
        "get_int",
        [](const T &self, size_t idx, uint8_t len) {
            if (idx + len - 1 >= self.size()) {
                throw std::out_of_range(std::to_string(idx)); }
            if (len > 64) {
                throw std::invalid_argument("len should be <= 64"); }
            return self.get_int(idx, len); },
        py::arg("idx"),
        py::arg("len") = 64,
        "Get the integer value of the binary string of length `len` "
        "starting at position `idx`.",
        py::call_guard<py::gil_scoped_release>());

    add_rank_support(m, cls, "_" + name, "", true, "0", "1", doc_rank);
    add_select_support(m, cls, "_" + name, "", true, "0", "1", doc_select);

    m.attr("all_immutable_bitvectors").attr("append")(cls);

    return cls;
}


template <class T>
inline
auto add_bitvector_class(py::module &m, const char* name,
                         const char* doc = nullptr,
                         const char* doc_rank = nullptr,
                         const char* doc_select = nullptr)
{
    return add_bitvector_class<T>(m, std::string(name),
                                  doc, doc_rank, doc_select);
}


template <class T>
inline
auto add_bitvector_class(py::module &m, const std::string& name,
                         const char* doc = nullptr,
                         const char* doc_rank = nullptr,
                         const char* doc_select = nullptr)
{
    return add_bitvector_class<T>(m, std::string(name), // i.e. copy name
                                  doc, doc_rank, doc_select);
}


template <uint32_t t_bs>
inline auto add_bit_vector_il(py::module& m)
{
    auto cls = add_bitvector_class<sdsl::bit_vector_il<t_bs>>(
        m,
        std::string("BitVectorInterLeaved") + std::to_string(t_bs),
        doc_bit_vector_il);

    m.attr("bit_vector_interleaved").attr("__setitem__")(t_bs, cls);

    return cls;
}


template <uint16_t t_bs=63, class t_rac=sdsl::int_vector<>, uint16_t t_k=32>
inline auto add_rrr_vector(py::module& m,
                           const char* name = "RamanRamanRaoVector")
{
    auto cls = add_bitvector_class<sdsl::rrr_vector<t_bs, t_rac, t_k>>(
        m,
        std::string(name) + std::to_string(t_bs),
                    doc_rrr_vector);

    m.attr("raman_raman_rao_vectors").attr("__setitem__")(t_bs, cls);

    return cls;
}


template <uint16_t t_bs=63, class t_rac=sdsl::int_vector<>, uint16_t t_k=32>
inline auto add_rrr_vector(py::module& m, const py::class_<t_rac>& py_rac,
                           const char* name = "RamanRamanRaoVector")
{
    auto cls_name = std::string(name) + py::cast<std::string>(py_rac.attr("__name__")) + std::to_string(t_bs);

    auto cls = add_bitvector_class<sdsl::rrr_vector<t_bs, t_rac, t_k>>(
            m, cls_name, doc_rrr_vector);

    m.attr("raman_raman_rao_vectors").attr("__setitem__")(cls_name, cls);

    return cls;
}


template <class Base=sdsl::bit_vector>
inline auto add_sd_vector(py::module& m, const char* name="SDVector")
{
    auto cls = add_bitvector_class<sdsl::sd_vector<Base>>(
        m,
        std::string(name),
        doc_sd_vector);

    m.attr("sparse_bit_vectors").attr("__setitem__")(name, cls);

    return cls;
}


template <class Base=sdsl::bit_vector>
inline auto add_sd_vector(py::module& m, const py::class_<Base>& base_cls,
                          const char* name="SDVector")
{
    auto cls_name = std::string(name) + py::cast<std::string>(
        base_cls.attr("__name__"));

    auto cls = add_bitvector_class<sdsl::sd_vector<Base>>(
        m,
        cls_name,
        doc_sd_vector);

    m.attr("sparse_bit_vectors").attr("__setitem__")(cls_name, cls);

    return cls;
}


template <uint32_t k_sblock_rate>
inline auto add_hyb_vector(py::module& m)
{
    auto cls = add_bitvector_class<sdsl::hyb_vector<k_sblock_rate>>(
        m, "HybVector" + std::to_string(k_sblock_rate),
        doc_hyb_vector);

    m.attr("hybrid_bit_vectors").attr("__setitem__")(k_sblock_rate, cls);

    return cls;
}


template <class B>
auto add_bitvectors(py::module& m, py::class_<B>& bit_vector_cls)
{
    add_bitvector_supports(m, bit_vector_cls);

    m.attr("all_immutable_bitvectors") = py::list();
    m.attr("bit_vector_interleaved") = py::dict();
    m.attr("raman_raman_rao_vectors") = py::dict();
    m.attr("sparse_bit_vectors") = py::dict();
    m.attr("hybrid_bit_vectors") = py::dict();

    auto bvil_classes = std::make_tuple(add_bit_vector_il<64>(m),
                                        add_bit_vector_il<128>(m),
                                        add_bit_vector_il<256>(m),
                                        add_bit_vector_il<512>(m));

    auto hyb_classes = std::make_tuple(
        add_hyb_vector<4>(m),
        add_hyb_vector<8>(m),
        add_hyb_vector<16>(m),
        add_hyb_vector<256>(m));

    auto rrr_classes = std::make_tuple(
        add_rrr_vector<3>(m),
        add_rrr_vector<15>(m),
        add_rrr_vector<63>(m),
        add_rrr_vector<256>(m));
        //add_rrr_vector<63, sdsl::wt_int<>>(m, "RamanRamanRaoWTVector"));

    auto sd_classes = std::make_tuple(
        add_sd_vector<>(m),
        add_sd_vector<sdsl::sd_vector<>>(m, "SDVectorSD"),
        add_sd_vector(m, std::get<1>(rrr_classes)));

    return std::make_tuple(
        std::tuple_cat(bvil_classes, rrr_classes, sd_classes, hyb_classes),
        std::make_tuple(  // propagate
            std::get<2>(rrr_classes),
            std::get<0>(sd_classes),
            std::get<3>(bvil_classes)
        )
    );

}
