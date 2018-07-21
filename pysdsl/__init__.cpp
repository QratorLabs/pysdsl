#include <cstdint>
#include <string>
#include <tuple>

#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <pybind11/pybind11.h>

#include "docstrings.hpp"
#include "types/bitvector.hpp"
#include "types/encodedvector.hpp"
#include "types/intvector.hpp"
#include "types/suffixarray.hpp"
#include "types/wavelet.hpp"


namespace py = pybind11;


PYBIND11_MODULE(pysdsl, m)
{
    m.doc() = "sdsl-lite bindings for python";

    auto iv_classes = add_int_vectors(m);

    py::class_<sdsl::int_vector<1>>& bit_vector_cls = std::get<1>(iv_classes);

    auto bit_vector_classes = std::make_tuple(bit_vector_cls);
    auto compressed_bit_vector_classes = add_bitvectors(m, bit_vector_cls);

    auto enc_classes = add_encoded_vectors(m);

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

    auto csa_classes = add_csa(m);

    for_each_in_tuple(iv_classes, make_inits_many_functor(iv_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(enc_classes));
    for_each_in_tuple(iv_classes,
                      make_inits_many_functor(compressed_bit_vector_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(wavelet_classes));

    for_each_in_tuple(enc_classes, make_inits_many_functor(iv_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(enc_classes, make_inits_many_functor(enc_classes));
    //for_each_in_tuple(enc_classes, make_inits_many_functor(wavelet_classes));
#endif

    for_each_in_tuple(compressed_bit_vector_classes,
                      make_inits_many_functor(bit_vector_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(compressed_bit_vector_classes, make_inits_many_functor(compressed_bit_vector_classes));
#endif

    for_each_in_tuple(wavelet_classes, make_inits_many_functor(iv_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(wavelet_classes, make_inits_many_functor(enc_classes));
    for_each_in_tuple(wavelet_classes, make_inits_many_functor(compressed_bit_vector_classes));
    for_each_in_tuple(wavelet_classes,
                      make_inits_many_functor(wavelet_classes));
#endif

    for_each_in_tuple(iv_classes, make_pysequence_init_functor());
    for_each_in_tuple(enc_classes, make_pysequence_init_functor());

    //for_each_in_tuple(sd_classes, make_pysequence_init_functor());

    for_each_in_tuple(wavelet_classes, make_pysequence_init_functor());

    for_each_in_tuple(csa_classes, make_pysequence_init_functor());
}
