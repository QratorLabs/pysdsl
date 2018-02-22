/*cppimport
<%
cfg['compiler_args'] = ['-std=c++14', '-fvisibility=hidden']
cfg['linker_args'] = ['-fvisibility=hidden']
cfg['include_dirs'] = ['sdsl-lite/include']
cfg['libraries'] = ['sdsl', 'divsufsort', 'divsufsort64']
cfg['dependencies'] = ['converters.hpp', 'pysequence.hpp', 'io.hpp',
                       'sizes.hpp', 'calc.hpp', 'docstrings.hpp',
                       'intvector.hpp']
%>
*/

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>

#include <pybind11/pybind11.h>

#include "calc.hpp"
#include "io.hpp"
#include "converters.hpp"
#include "sizes.hpp"
#include "docstrings.hpp"
#include "intvector.hpp"


namespace py = pybind11;


using sdsl::int_vector;


template <class T>
auto add_compressed_class(py::module &m, const std::string& name,
                          const char* doc = nullptr)
{
    auto cls = py::class_<T>(m, name.c_str()).def(py::init());

    add_sizes(cls);
    add_description(cls);
    add_serialization(cls);
    add_to_string(cls);

    add_std_algo(cls);

    if (doc) cls.doc() = doc;

    return cls;
}


template <class T>
class support_helper
{
private:
    const sdsl::bit_vector& m_vec;
    const T m_support;
public:
    typedef T type;

    support_helper(const sdsl::bit_vector& vec, const T&& support):
        m_vec(vec),
        m_support(std::move(support))
    {}
    auto size() const { return m_vec.size(); }
    auto operator()(size_t idx) const { return m_support(idx); }

    operator const T&() const { return m_support; }
};


template <class Base>
auto add_support_class(py::module &m,
                       const std::string&& name,
                       const std::string&& method_name,
                       const std::string&& doc_call,
                       const char* doc = nullptr)
{
    auto cls = py::class_<Base>(m, name.c_str());

    cls.def(
        method_name.c_str(),
        [](const Base& self, size_t idx)
        {
            if (idx >= self.size())
            {
                throw std::out_of_range(std::to_string(idx));
            }
            return self(idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        py::arg("idx"),
        doc_call.c_str()
    );
    cls.attr("__call__") = cls.attr(method_name.c_str());

    add_description(cls);

    if (doc) cls.doc() = doc;

     return cls;
}


template <class S, class T>
decltype(auto) bind_support(const S*,
                            py::class_<T>& cls, const std::string& call_name,
                            const char* alt_name=nullptr)
{
    cls.def(
        call_name.c_str(),
        [](T& self)
        {
            S support;
            sdsl::util::init_support(support, &self);

            return support;
        },
        py::keep_alive<0, 1>()
    );

    if (alt_name) cls.attr(alt_name) = cls.attr(call_name.c_str());

    return cls;
}


template <class S, class T>
decltype(auto) bind_support(const support_helper<S>*,
                            py::class_<T>& cls, const std::string& call_name,
                            const char* alt_name=nullptr)
{
    cls.def(
        call_name.c_str(),
        [](T& self)
        {
            S support;
            sdsl::util::init_support(support, &self);

            return support_helper<S>(self, std::move(support));
        },
        py::keep_alive<0, 1>()
    );

    if (alt_name) cls.attr(alt_name) = cls.attr(call_name.c_str());

    return cls;
}


template <class T,
          class R0=typename T::rank_0_type, class R1=typename T::rank_1_type>
auto add_rank_support(py::module &m, py::class_<T>& cls,
                      const std::string& base_name,
                      const char* suffix = "",
                      bool defaults = true,
                      const std::string s0 = "0", const std::string s1 = "1",
                      const char* doc_rank = nullptr)
{
    add_support_class<R0>(m, base_name + "Rank" + suffix + "_" + s0, "rank",
                          "number of patterns `" + s0 + "` in the prefix "
                          "[0..idx) in supported vector", doc_rank);
    bind_support((R0 *)nullptr, cls,
                 std::string("init_rank") + suffix + "_" + s0);

    add_support_class<R1>(m, base_name + "Rank" + suffix + "_" + s1, "rank",
                          "number of patterns `" + s1 + "` in the prefix "
                          "[0..idx) in supported vector", doc_rank);
    bind_support((R1 *)nullptr, cls,
                 std::string("init_rank") + suffix + "_" + s1,
                 defaults ?
                    (std::string("init_rank") + suffix).c_str() :
                    nullptr);

    return cls;
}


template <class T,
          class S0=typename T::select_0_type,
          class S1=typename T::select_1_type>
auto add_select_support(py::module &m, py::class_<T>& cls,
                        const std::string& base_name,
                        const char* suffix = "",
                        bool defaults = true,
                        const std::string s0 = "0", const std::string s1 = "1",
                        const char* doc_select = nullptr)
{
    add_support_class<S0>(m, base_name + "Select" + suffix + "_" + s0, "select",
                          "position of the idx-th pattern `" + s0 +
                          "` in supported vector", doc_select);
    bind_support((S0 *)nullptr,
                 cls, std::string("init_select") + suffix + "_" + s0);

    add_support_class<S1>(m, base_name + "Select" + suffix + "_" + s1, "select",
                          "position of the idx-th pattern `" + s1 +
                          "` in supported vector", doc_select);
    bind_support((S1 *)nullptr, cls,
                 std::string("init_select") + suffix + "_" + s1,
                 defaults ?
                    (std::string("init_select") + suffix).c_str() :
                    nullptr);

    return cls;
}


template <class T>
auto add_bitvector_class(py::module &m, const std::string&& name,
                         const char* doc = nullptr,
                         const char* doc_rank = nullptr,
                         const char* doc_select = nullptr)
{
    auto cls = add_compressed_class<T>(m, name, doc);

    cls.def(
        "get_int",
        [](const T &self, size_t idx, uint8_t len) {
            if (idx + len - 1 >= self.size())
            {
                throw std::out_of_range(std::to_string(idx));
            }
            if (len > 64)
            {
                throw std::invalid_argument("len should be <= 64");
            }
            return self.get_int(idx, len);
        },
        py::arg("idx"),
        py::arg("len") = 64,
        "Get the integer value of the binary string of length `len` "
        "starting at position `idx`.",
        py::call_guard<py::gil_scoped_release>()
    );

    add_rank_support(m, cls, name, "", true, "0", "1", doc_rank);
    add_select_support(m, cls, name, "", true, "0", "1", doc_select);

    return cls;
}


class add_enc_coders_functor
{
public:
    add_enc_coders_functor(py::module& m): m(m) {}

    template <typename Coder, uint32_t t_dens=128, uint8_t t_width=0>
    decltype(auto) operator()(const std::pair<const char*, Coder> &t)
    {
        typedef sdsl::enc_vector<Coder, t_dens, t_width> enc;
        auto cls = add_compressed_class<enc>(
            m,
            std::string("EncVector") + std::get<0>(t),
            "A vector `v` is stored more space-efficiently by "
            "self-delimiting coding the deltas v[i+1]-v[i] (v[-1]:=0)."
        )
        .def_property_readonly("sample_dens", &enc::get_sample_dens)
        .def(
            "sample",
            [](const enc &self, typename enc::size_type i) {
                if (i >= self.size() / self.get_sample_dens())
                {
                    throw std::out_of_range(std::to_string(i));
                }
                return self.sample(i);
            },
            "Returns the i-th sample of the compressed vector"
            "i: The index of the sample. 0 <= i < size()/get_sample_dens()",
            py::call_guard<py::gil_scoped_release>()
        )
        ;

        return cls;
    }

private:
    py::module& m;
};


class add_vlc_coders_functor
{
public:
    add_vlc_coders_functor(py::module& m): m(m) {}

    template <typename Coder, uint32_t t_dens = 128, uint8_t t_width = 0>
    decltype(auto) operator()(const std::pair<const char*, Coder> &t)
    {
        typedef sdsl::vlc_vector<Coder, t_dens, t_width> vlc;

        auto cls = add_compressed_class<vlc>(
            m,
            std::string("VlcVector") + std::get<0>(t),
            "A vector which stores the values with variable length codes."
        )
        .def_property_readonly("sample_dens", &vlc::get_sample_dens)
        ;

        return cls;
    }

private:
    py::module& m;
};


PYBIND11_MODULE(pysdsl, m)
{
    m.doc() = "sdsl-lite bindings for python";

    auto iv_classes = add_int_vectors(m);

    py::class_<int_vector<1>>& bit_vector_cls = std::get<1>(iv_classes);
    auto bit_vector_classes = std::make_tuple(bit_vector_cls);

    auto const coders = std::make_tuple(
        std::make_pair("EliasDelta", sdsl::coder::elias_delta()),
        std::make_pair("EliasGamma", sdsl::coder::elias_gamma()),
        std::make_pair("Fibonacci", sdsl::coder::fibonacci()),
        std::make_pair("Comma2", sdsl::coder::comma<2>()),
        std::make_pair("Comma4", sdsl::coder::comma<4>())
    );

    auto enc_classes = for_each_in_tuple(coders, add_enc_coders_functor(m));
    auto vlc_classes = for_each_in_tuple(coders, add_vlc_coders_functor(m));
    auto dac_classes = std::make_tuple(
        add_compressed_class<sdsl::dac_vector<>>(m, "DacVector",
                                                 doc_dac_vector)
            .def_property_readonly("levels", &sdsl::dac_vector<>::levels),
        add_compressed_class<sdsl::dac_vector_dp<>>(m, "DacVectorDP",
                                                    doc_dac_vector_dp)
            .def("cost", &sdsl::dac_vector_dp<>::cost, py::arg("n"),
                 py::arg("m"))
            .def_property_readonly("levels", &sdsl::dac_vector_dp<>::levels)
    );

    auto bvil_classes = std::make_tuple(
        add_bitvector_class<sdsl::bit_vector_il<64>>(m, "BitVectorIL64",
                                                     doc_bit_vector_il),
        add_bitvector_class<sdsl::bit_vector_il<128>>(m, "BitVectorIL128",
                                                      doc_bit_vector_il),
        add_bitvector_class<sdsl::bit_vector_il<256>>(m, "BitVectorIL256",
                                                      doc_bit_vector_il),
        add_bitvector_class<sdsl::bit_vector_il<512>>(m, "BitVectorIL512",
                                                      doc_bit_vector_il)
    );

    auto rrr_classes = std::make_tuple(
        add_bitvector_class<sdsl::rrr_vector<3>>(m, "RRRVector3",
                                                  doc_rrr_vector),
        add_bitvector_class<sdsl::rrr_vector<15>>(m, "RRRVector15",
                                                  doc_rrr_vector),
        add_bitvector_class<sdsl::rrr_vector<63>>(m, "RRRVector63",
                                                  doc_rrr_vector),
        add_bitvector_class<sdsl::rrr_vector<256>>(m, "RRRVector256",
                                                   doc_rrr_vector)
    );

    auto sd_classes = std::make_tuple(
        add_bitvector_class<sdsl::sd_vector<>>(m, std::string("SDVector"),
                                               doc_sd_vector)
    );

    auto hyb_classes = std::make_tuple(
        add_bitvector_class<sdsl::hyb_vector<8>>(m, std::string("HybVector8"),
                                                 doc_hyb_vector),
        add_bitvector_class<sdsl::hyb_vector<16>>(m, std::string("HybVector16"),
                                                  doc_hyb_vector)
    );

    add_rank_support<sdsl::bit_vector,
                    sdsl::rank_support_v5<0, 1>,
                    sdsl::rank_support_v5<1, 1>>(
        m, bit_vector_cls, "BitVector", "V5", true, "0", "1", doc_rank_v5
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_v5<00, 2>,
                     sdsl::rank_support_v5<01, 2>>(
        m, bit_vector_cls, "BitVector", "V5", false, "00", "01", doc_rank_v5
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_v5<10, 2>,
                     sdsl::rank_support_v5<11, 2>>(
        m, bit_vector_cls, "BitVector", "V5", false, "10", "11", doc_rank_v5
    );
    add_rank_support<sdsl::bit_vector,
                    sdsl::rank_support_v<0, 1>,
                    sdsl::rank_support_v<1, 1>>(
        m, bit_vector_cls, "BitVector", "V", false, "0", "1", doc_rank_v
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_v<00, 2>,
                     sdsl::rank_support_v<01, 2>>(
        m, bit_vector_cls, "BitVector", "V", false, "00", "01", doc_rank_v
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_v<10, 2>,
                     sdsl::rank_support_v<11, 2>>(
        m, bit_vector_cls, "BitVector", "V", false, "10", "11", doc_rank_v
    );
    bit_vector_cls.attr("init_rank") = bit_vector_cls.attr("init_rankV5");
    bit_vector_cls.attr("init_rank_0") = bit_vector_cls.attr("init_rankV5_0");
    bit_vector_cls.attr("init_rank_1") = bit_vector_cls.attr("init_rankV5_1");

    add_select_support<sdsl::bit_vector,
                       support_helper<sdsl::select_support_mcl<0, 1>>,
                       support_helper<sdsl::select_support_mcl<1, 1>>>(
        m, bit_vector_cls, "BitVector", "MCL", true, "0", "1", doc_select_mcl
    );
    add_select_support<sdsl::bit_vector,
                       support_helper<sdsl::select_support_mcl<10, 2>>,
                       support_helper<sdsl::select_support_mcl<11, 2>>>(
        m, bit_vector_cls, "BitVector", "MCL", false, "10", "11", doc_select_mcl
    );
    bit_vector_cls.attr("init_select") = bit_vector_cls.attr("init_selectMCL");
    bit_vector_cls.attr("init_select_0") = bit_vector_cls.attr(
        "init_selectMCL_0"
    );
    bit_vector_cls.attr("init_select_1") = bit_vector_cls.attr(
        "init_selectMCL_1"
    );

    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_scan<0, 1>,
                     sdsl::rank_support_scan<1, 1>>(
        m, bit_vector_cls, "BitVector", "Scan", false, "0", "1", doc_rank_scan
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_scan<00, 2>,
                     sdsl::rank_support_scan<01, 2>>(
        m, bit_vector_cls, "BitVector", "Scan", false, "00", "01",
        doc_rank_scan
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_scan<10, 2>,
                     sdsl::rank_support_scan<11, 2>>(
        m, bit_vector_cls, "BitVector", "Scan", false, "10", "11",
        doc_rank_scan
    );

    add_select_support<sdsl::bit_vector,
                       support_helper<sdsl::select_support_scan<0, 1>>,
                       support_helper<sdsl::select_support_scan<1, 1>>>(
        m, bit_vector_cls, "BitVector", "Scan", false, "0", "1",
        doc_select_scan
    );
    add_select_support<sdsl::bit_vector,
                       support_helper<sdsl::select_support_scan<10, 2>>,
                       support_helper<sdsl::select_support_scan<01, 2>>>(
        m, bit_vector_cls, "BitVector", "Scan", false, "10", "01",
        doc_select_scan
    );

    for_each_in_tuple(iv_classes, make_inits_many_functor(iv_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(enc_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(vlc_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(dac_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(iv_classes, make_inits_many_functor(hyb_classes));

    for_each_in_tuple(enc_classes, make_inits_many_functor(iv_classes));
    for_each_in_tuple(enc_classes, make_inits_many_functor(enc_classes));
    for_each_in_tuple(enc_classes, make_inits_many_functor(vlc_classes));
    for_each_in_tuple(enc_classes, make_inits_many_functor(dac_classes));

    for_each_in_tuple(vlc_classes, make_inits_many_functor(iv_classes));
    for_each_in_tuple(vlc_classes, make_inits_many_functor(enc_classes));
    for_each_in_tuple(vlc_classes, make_inits_many_functor(vlc_classes));
    for_each_in_tuple(vlc_classes, make_inits_many_functor(dac_classes));

    for_each_in_tuple(dac_classes, make_inits_many_functor(iv_classes));
    for_each_in_tuple(dac_classes, make_inits_many_functor(enc_classes));
    for_each_in_tuple(dac_classes, make_inits_many_functor(vlc_classes));
    for_each_in_tuple(dac_classes, make_inits_many_functor(dac_classes));

    for_each_in_tuple(bvil_classes,
                      make_inits_many_functor(bit_vector_classes));
    for_each_in_tuple(bvil_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(bvil_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(bvil_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(bvil_classes, make_inits_many_functor(hyb_classes));

    for_each_in_tuple(rrr_classes, make_inits_many_functor(bit_vector_classes));
    for_each_in_tuple(rrr_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(rrr_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(rrr_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(rrr_classes, make_inits_many_functor(hyb_classes));

    for_each_in_tuple(sd_classes, make_inits_many_functor(bit_vector_classes));
    for_each_in_tuple(sd_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(sd_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(sd_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(sd_classes, make_inits_many_functor(hyb_classes));

    for_each_in_tuple(hyb_classes, make_inits_many_functor(bit_vector_classes));
    for_each_in_tuple(hyb_classes, make_inits_many_functor(bvil_classes));
    for_each_in_tuple(hyb_classes, make_inits_many_functor(rrr_classes));
    for_each_in_tuple(hyb_classes, make_inits_many_functor(sd_classes));
    for_each_in_tuple(hyb_classes, make_inits_many_functor(hyb_classes));

    for_each_in_tuple(iv_classes, make_pysequence_init_functor());
    for_each_in_tuple(enc_classes, make_pysequence_init_functor());
    for_each_in_tuple(vlc_classes, make_pysequence_init_functor());
    for_each_in_tuple(dac_classes, make_pysequence_init_functor());

    //for_each_in_tuple(sd_classes, make_pysequence_init_functor());
}
