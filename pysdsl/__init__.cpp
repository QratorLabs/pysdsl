#include <cstdint>
#include <string>
#include <tuple>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <pybind11/pybind11.h>

#include "calc.hpp"
#include "io.hpp"
#include "operations/creation.hpp"
#include "operations/sizes.hpp"
#include "docstrings.hpp"
#include "types/encodedvector.hpp"
#include "types/intvector.hpp"
#include "supports.hpp"
#include "wavelet.hpp"


namespace py = pybind11;


template <class Sequence, typename T = typename Sequence::value_type>
inline
auto add_compressed_class(py::module &m, const std::string& name,
                          const char* doc = nullptr)
{
    auto cls = py::class_<Sequence>(m, name.c_str()).def(py::init());

    add_sizes(cls);
    add_description(cls);
    add_serialization(cls);
    add_to_string(cls);

    add_iteration(cls);
    add_std_algo<Sequence, T>(cls);

    if (doc) cls.doc() = doc;

    return cls;
}


template <class T>
inline
auto add_bitvector_class(py::module &m, const std::string&& name,
                         const char* doc = nullptr,
                         const char* doc_rank = nullptr,
                         const char* doc_select = nullptr)
{
    auto cls = add_compressed_class<T, bool>(m, name, doc);

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

    return cls;
}


template <class T>
auto add_csa(py::module& m, const char* name, const char* doc = nullptr)
{
    auto cls = py::class_<T>(m, name);

    cls.def_property_readonly("char2comp", [] (const T& self ) { return self.char2comp; });
    cls.def_property_readonly("comp2char", [] (const T& self ) { return self.comp2char; });
    cls.def_property_readonly("sigma", [] (const T& self ) { return self.sigma; });

    add_sizes(cls);
    add_description(cls);
    add_serialization(cls);
    add_to_string(cls);

    add_iteration(cls);
    add_std_algo<T>(cls);

    if (doc) cls.doc() = doc;

    return cls;
}


PYBIND11_MODULE(pysdsl, m)
{
    m.doc() = "sdsl-lite bindings for python";

    auto iv_classes = add_int_vectors(m);

    py::class_<sdsl::int_vector<1>>& bit_vector_cls = std::get<1>(iv_classes);
    add_bitvector_supports(m, bit_vector_cls);

    auto bit_vector_classes = std::make_tuple(bit_vector_cls);

    auto enc_classes = add_encoded_vectors(m);

    auto bvil_classes = std::make_tuple(
        add_bitvector_class<sdsl::bit_vector_il<64>>(m, "BitVectorIL64",
                                                     doc_bit_vector_il),
        add_bitvector_class<sdsl::bit_vector_il<128>>(m, "BitVectorIL128",
                                                      doc_bit_vector_il),
        add_bitvector_class<sdsl::bit_vector_il<256>>(m, "BitVectorIL256",
                                                      doc_bit_vector_il),
        add_bitvector_class<sdsl::bit_vector_il<512>>(m, "BitVectorIL512",
                                                      doc_bit_vector_il));

    auto rrr_classes = std::make_tuple(
        add_bitvector_class<sdsl::rrr_vector<3>>(m, "RRRVector3",
                                                  doc_rrr_vector),
        add_bitvector_class<sdsl::rrr_vector<15>>(m, "RRRVector15",
                                                  doc_rrr_vector),
        add_bitvector_class<sdsl::rrr_vector<63>>(m, "RRRVector63",
                                                  doc_rrr_vector),
        add_bitvector_class<sdsl::rrr_vector<256>>(m, "RRRVector256",
                                                   doc_rrr_vector));

    auto sd_classes = std::make_tuple(
        add_bitvector_class<sdsl::sd_vector<>>(m, std::string("SDVector"),
                                               doc_sd_vector),
        add_bitvector_class<sdsl::sd_vector<sdsl::sd_vector<>>>(m, std::string("SDVectorR"),
                                               doc_sd_vector));

    auto hyb_classes = std::make_tuple(
        add_bitvector_class<sdsl::hyb_vector<8>>(m, std::string("HybVector8"),
                                                 doc_hyb_vector),
        add_bitvector_class<sdsl::hyb_vector<16>>(m, std::string("HybVector16"),
                                                  doc_hyb_vector));

    auto wavelet_classes = std::make_tuple(
        add_wavelet_class<sdsl::wt_int<sdsl::bit_vector>>(m, "WtInt",
                                                          doc_wtint),
        add_wavelet_class<sdsl::wt_int<sdsl::bit_vector_il<>>>(m, "WtIntIL",
                                                               doc_wtint),
        add_wavelet_class<sdsl::wt_int<sdsl::rrr_vector<>>>(m, "WtIntRRR",
                                                            doc_wtint),
        add_wavelet_class<sdsl::wt_int<sdsl::sd_vector<>>>(m, "WtIntSD",
                                                           doc_wtint),
        //add_wavelet_class<sdsl::wt_int<sdsl::hyb_vector<>>>(m, "WtIntHyb",
        //                                                    doc_wtint),
        add_wavelet_class<sdsl::wt_gmr_rs<>>(m, "WtGMRrs", doc_wt_gmr_rs),
        add_wavelet_class<sdsl::wt_gmr<>>(m, "WtGMR", doc_wt_gmr),
        add_wavelet_class<sdsl::wt_ap<>>(m, "WtAP", doc_wt_ap),
        add_wavelet_class<sdsl::wt_huff<>>(m, "WtHuff", doc_wt_huff),
        add_wavelet_class<sdsl::wt_huff_int<>>(m, "WtHuffInt", doc_wt_huff),
        add_wavelet_class<sdsl::wm_int<>>(m, "WmInt", doc_wm_int),
        add_wavelet_class<sdsl::wt_blcd<>>(m, "WtBlcd", doc_wt_blcd),
        add_wavelet_class<sdsl::wt_blcd_int<>>(m, "WtBlcdInt", doc_wt_blcd),
        add_wavelet_class<sdsl::wt_hutu<>>(m, "WtHutu", doc_wt_hutu),
        add_wavelet_class<sdsl::wt_hutu_int<>>(m, "WtHutuInt", doc_wt_hutu));

    auto csa_classes = std::make_tuple(
        add_csa<sdsl::csa_bitcompressed<>>(
            m, "SuffixArrayBitcompressed", doc_csa),
        add_csa<sdsl::csa_sada<>>(m, "SuffixArraySADA", doc_sada),
        add_csa<sdsl::csa_sada_int<>>(m, "SuffixArraySADAint", doc_sada),
        add_csa<sdsl::csa_wt<>>(m, "SuffixArrayWT", doc_csa_wt),
        add_csa<sdsl::csa_wt_int<>>(m, "SuffixArrayWTint", doc_csa_wt));

    for_each_in_tuple(iv_classes, make_inits_many_functor(iv_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(enc_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(hyb_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(wavelet_classes));

    for_each_in_tuple(enc_classes, make_inits_many_functor(iv_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(enc_classes, make_inits_many_functor(enc_classes));
    //for_each_in_tuple(enc_classes, make_inits_many_functor(wavelet_classes));
#endif

    for_each_in_tuple(bvil_classes,
                      make_inits_many_functor(bit_vector_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(bvil_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(bvil_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(bvil_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(bvil_classes, make_inits_many_functor(hyb_classes));
#endif

    for_each_in_tuple(rrr_classes, make_inits_many_functor(bit_vector_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(rrr_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(rrr_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(rrr_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(rrr_classes, make_inits_many_functor(hyb_classes));
#endif

    for_each_in_tuple(sd_classes, make_inits_many_functor(bit_vector_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(sd_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(sd_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(sd_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(sd_classes, make_inits_many_functor(hyb_classes));
#endif

    for_each_in_tuple(hyb_classes, make_inits_many_functor(bit_vector_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(hyb_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(hyb_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(hyb_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(hyb_classes, make_inits_many_functor(hyb_classes));
#endif

    for_each_in_tuple(wavelet_classes, make_inits_many_functor(iv_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(wavelet_classes, make_inits_many_functor(enc_classes));
    for_each_in_tuple(wavelet_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(wavelet_classes,
                      make_inits_many_functor(wavelet_classes));
#endif

    for_each_in_tuple(iv_classes, make_pysequence_init_functor());
    for_each_in_tuple(enc_classes, make_pysequence_init_functor());

    //for_each_in_tuple(sd_classes, make_pysequence_init_functor());

    for_each_in_tuple(wavelet_classes, make_pysequence_init_functor());

    for_each_in_tuple(csa_classes, make_pysequence_init_functor());
}
