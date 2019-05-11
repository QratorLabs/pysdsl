#include <cstdint>
#include <string>
#include <tuple>
#include <stdexcept>

#define assert(x) if(!x) {throw std::runtime_error("assertion failed");}

#include <sdsl/vectors.hpp>

#include <pybind11/pybind11.h>

#include "docstrings.hpp"
#include "operations/creation.hpp"
#include "types/bitvector.hpp"
#include "types/encodedvector.hpp"
#include "types/intvector.hpp"
#include "types/suffixarray.hpp"
#include "types/wavelet.hpp"

namespace py = pybind11;


PYBIND11_MODULE(pysdsl, m)
{
    m.doc() = "sdsl-lite bindings for python";

    auto iv_classes_tpl = add_int_vectors(m);

    auto iv_classes = std::get<0>(iv_classes_tpl);
    auto iv_classes_as_params = std::get<1>(iv_classes_tpl);

    py::class_<sdsl::int_vector<1>>& bit_vector_cls = std::get<1>(iv_classes);

    auto bit_vector_classes = std::make_tuple(bit_vector_cls);

    auto tmp = add_bitvectors(m, bit_vector_cls);
    auto compressed_bit_vector_classes = std::get<0>(tmp);
    auto cbv_propagate = std::get<1>(tmp);

    auto enc_classes = add_encoded_vectors(m);

    auto wavelet_classes = add_wavelet(m, cbv_propagate);

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
    for_each_in_tuple(compressed_bit_vector_classes, 
                      make_inits_many_functor(compressed_bit_vector_classes));
#endif

    for_each_in_tuple(wavelet_classes, make_inits_many_functor(iv_classes));
#ifndef NOCROSSCONSTRUCTORS
    for_each_in_tuple(wavelet_classes, make_inits_many_functor(enc_classes));
    for_each_in_tuple(wavelet_classes,
                      make_inits_many_functor(wavelet_classes));
#endif

    for_each_in_tuple(iv_classes, make_pysequence_init_functor());
    for_each_in_tuple(enc_classes, make_pysequence_init_functor());
    //for_each_in_tuple(compressed_bit_vector_classes,
    //                  make_pysequence_init_functor());

    for_each_in_tuple(wavelet_classes, make_pysequence_init_functor());
    for_each_in_tuple(csa_classes, make_pysequence_init_functor());
}
