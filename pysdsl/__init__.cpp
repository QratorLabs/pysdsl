/*cppimport
<%
cfg['compiler_args'] = ['-std=c++14', '-fvisibility=hidden']
cfg['linker_args'] = ['-fvisibility=hidden']
cfg['include_dirs'] = ['sdsl-lite/include']
cfg['libraries'] = ['sdsl', 'divsufsort', 'divsufsort64']
cfg['dependencies'] = ['converters.hpp', 'pysequence.hpp',
                       'sizes.hpp', 'calc.hpp', 'docstrings.hpp']
%>
*/

#include <cstdint>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>

#include <pybind11/pybind11.h>

#include "calc.hpp"
#include "converters.hpp"
#include "sizes.hpp"
#include "docstrings.hpp"


namespace py = pybind11;


using sdsl::int_vector;


template <class T>
auto add_io(py::class_<T>& cls)
{
    cls.def(
        "__str__",
        [](const T &self) {
            const size_t max_output = 100;

            std::stringstream fout;
            fout << '[';
            size_t count = 0;
            for (auto i: self)
            {
                if (count) fout << ", ";

                fout << i;

                if (count >= max_output)
                {
                    fout << ", ...(" << self.size() - count - 1 << " more)";
                    break;
                }
                count++;
            }
            fout << ']';
            return fout.str();
        }
    );

    cls.def(
        "store_to_file",
        [](const T &self, const std::string& file_name) {
            return sdsl::store_to_file(self, file_name);
        },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>()
    );

    cls.def_static(
        "load_from_file",
        [](const std::string& file_name) {
            T self;
            if (sdsl::load_from_file(self, file_name))
            {
                return self;
            }
            throw std::exception();
        },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>()
    );

    cls.def(
        "store_to_checked_file",
        [](const T &self, const std::string& file_name) {
            return sdsl::store_to_checked_file(self, file_name);
        },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>()
    );

    cls.def_static(
        "load_from_checkded_file",
        [](const std::string& file_name) {
            T self;
            if (sdsl::load_from_checked_file(self, file_name))
            {
                return self;
            }
            throw std::exception();
        },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>()
    );

    cls.def(
        "write_structure_json",
        [](const T& self, const std::string& file_name) {
            std::ofstream fout;
            fout.open(file_name, std::ios::out | std::ios::binary);
            if (!fout.good()) throw std::runtime_error("Can't write to file");
            sdsl::write_structure<sdsl::JSON_FORMAT>(self, fout);
            if (!fout.good()) throw std::runtime_error("Error during write");
            fout.close();
        },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>()
    );

    cls.def_property_readonly(
        "structure_json",
        [](const T& self) {
            std::stringstream fout;
            sdsl::write_structure<sdsl::JSON_FORMAT>(self, fout);
            return fout.str();
        },
        py::call_guard<py::gil_scoped_release>()
    );

    cls.def_property_readonly(
        "structure",
        [](const T& self) {
            std::stringstream fout;
            sdsl::write_structure<sdsl::JSON_FORMAT>(self, fout);
            auto json = py::module::import("json");
            return json.attr("loads")(fout.str());
        }
    );

    return cls;
}


template <class T, typename S = uint64_t>
auto add_int_class(py::module &m, const char *name, const char *doc = nullptr)
{
    auto cls = py::class_<T>(m, name)
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
                if (position >= self.size())
                {
                    throw std::out_of_range(std::to_string(position));
                }
                self[position] = value;
            }
        )

        .def("set_to_id",
             [](T &self) { sdsl::util::set_to_id(self); },
             py::call_guard<py::gil_scoped_release>(),
             "Sets each entry of the vector at position `i` to value `i`")
        .def("set_to_value",
             [](T &self, S value) { sdsl::util::set_to_value(self, value); },
             py::arg("k"),
             doc_set_to_value,
             py::call_guard<py::gil_scoped_release>()
        )
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
                sdsl::util::set_random_bits(self, seed);
            },
            py::arg_v(
                "seed",
                0,
                "If seed = 0, the time is used to initialize the pseudo "
                "random number generator, otherwise the seed parameter is used."
            ),
            "Sets all bits of the int_vector to pseudo-random bits.",
            py::call_guard<py::gil_scoped_release>()
        )
        .def_static(
            "rnd_positions",
            [](uint8_t log_s, uint64_t mod, uint64_t seed) {
                uint64_t mask;

                auto res = sdsl::util::rnd_positions<T>(log_s, mask, mod, seed);

                return std::make_tuple(res, mask);
            },
            py::arg("log_s"), py::arg("mod") = 0, py::arg("seed") = 0,
            "Create `2**{log_s}` random integers mod `mod` with seed `seed`",
            py::call_guard<py::gil_scoped_release>()
        )
        .def(
            "__imod__",
            [](T &self, uint64_t m) {
                sdsl::util::mod(self, m);
                return self;
            },
            py::is_operator()
        )

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
                if (idx >= self.bit_size())
                {
                    throw std::out_of_range(std::to_string(idx));
                }
                return sdsl::util::next_bit(self, idx);
            },
            py::arg("idx"),
            "Get the smallest position `i` >= `idx` where a bit is set",
            py::call_guard<py::gil_scoped_release>()
        )
        .def(
            "prev_bit",
            [](const T &self, size_t idx) {
                if (idx >= self.bit_size())
                {
                    throw std::out_of_range(std::to_string(idx));
                }
                return sdsl::util::prev_bit(self, idx);
            },
            py::arg("idx"),
            "Get the largest position `i` <= `idx` where a bit is set",
            py::call_guard<py::gil_scoped_release>()
        )
    ;

    add_sizes(cls);
    add_io(cls);

    add_std_algo(cls);

    if (doc) cls.doc() = doc;

    return cls;
}


template <class T>
auto add_compressed_class(py::module &m, const std::string& name,
                          const char* doc = nullptr)
{
    auto cls = py::class_<T>(m, name.c_str()).def(py::init());

    add_sizes(cls);
    add_io(cls);

    add_std_algo(cls);

    if (doc) cls.doc() = doc;

    return cls;
}


template <class T>
auto add_bitvector_class(py::module &m, const std::string&& name,
                         const char* doc = nullptr)
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

    return cls;
}


template <class T>
class support_helper
{
private:
    const sdsl::bit_vector& m_vec;
    const T& m_support;
public:
    typedef T type;

    support_helper(const sdsl::bit_vector& vec, const T& support):
        m_vec(vec),
        m_support(support)
    {}
    auto size() const { return m_vec.size(); }
    auto operator()(size_t idx) const { return m_support(idx); }
};


template <class S>
class bind_support_functor
{
public:
    bind_support_functor(const char* call_name): m_call_name(call_name)
    {}

    template <class T>
    decltype(auto) operator() (py::class_<T>& cls) const
    {
        cls.def(
            m_call_name,
            [](T& self)
            {
                S support;
                sdsl::util::init_support(support, &self);

                return support;
            },
            py::keep_alive<0, 1>()
        );
        return cls;
    }

private:
    const char* m_call_name;
};

template <class S>
class bind_support_functor<support_helper<S>>
{
public:
    bind_support_functor(const char* call_name): m_call_name(call_name)
    {}

    template <class T>
    decltype(auto) operator() (py::class_<T>& cls) const
    {
        cls.def(
            m_call_name,
            [](T& self)
            {
                S support;
                sdsl::util::init_support(support, &self);

                return support_helper<S>(self, support);
            },
            py::keep_alive<0, 1>()
        );
        return cls;
    }

private:
    const char* m_call_name;
};


template <class Base, class...T>
auto add_support_class(py::module &m, const std::string&& name,
                       std::tuple<T...>& bind,
                       const char* bind_name,
                       const char* doc = nullptr)
{
    auto cls = py::class_<Base>(m, name.c_str());

    cls.def(
        "__call__",
        [](const Base& self, size_t idx)
        {
            if (idx >= self.size())
            {
                throw std::out_of_range(std::to_string(idx));
            }
            return self(idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        py::arg("idx")
    );

    if (doc) cls.doc() = doc;

    for_each_in_tuple(bind, bind_support_functor<Base>(bind_name));

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

    auto iv_classes = std::make_tuple(
        add_int_class<int_vector<0>>(m, "IntVector", doc_int_vector)
            .def(
                py::init([](size_t size,
                            uint64_t default_value,
                            uint8_t bit_width) {
                    return int_vector<0>(size, default_value, bit_width);
                }),
                py::arg("size") = 0,
                py::arg("default_value") = 0,
                py::arg("bit_width") = 64,
                py::call_guard<py::gil_scoped_release>())
            .def(
                "expand_width",
                [](int_vector<0> &self, size_t width) {
                    sdsl::util::expand_width(self, width);
                },
                "Expands the integer width to new_width >= v.width().",
                py::call_guard<py::gil_scoped_release>()
            )
            .def("bit_compress",
                [](int_vector<0> &self) { sdsl::util::bit_compress(self); },
                doc_bit_compress,
                py::call_guard<py::gil_scoped_release>()),

        add_int_class<int_vector<1>, bool>(m, "BitVector")
            .def(py::init([](size_t size, bool default_value) {
                    return int_vector<1>(size, default_value, 1);
                }), py::arg("size") = 0, py::arg("default_value") = false)
            .def("flip", &int_vector<1>::flip,
                 "Flip all bits of bit_vector",
                 py::call_guard<py::gil_scoped_release>()),

        add_int_class<int_vector<4>, uint8_t>(m, "Int4Vector")
            .def(py::init([](size_t size, uint8_t default_value) {
                    return int_vector<4>(size, default_value, 4);
                }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<int_vector<8>, uint8_t>(m, "Int8Vector")
            .def(py::init([](size_t size, uint8_t default_value) {
                    return int_vector<8>(size, default_value, 8);
                }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<int_vector<16>, uint16_t>(m, "Int16Vector")
            .def(py::init([](size_t size, uint16_t default_value) {
                     return int_vector<16>(size, default_value, 16);
                 }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<int_vector<24>, uint32_t>(m, "Int24Vector")
            .def(py::init([](size_t size, uint32_t default_value) {
                     return int_vector<24>(size, default_value, 24);
                 }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<int_vector<32>, uint32_t>(m, "Int32Vector")
            .def(py::init([](size_t size, uint32_t default_value) {
                    return int_vector<32>(size, default_value, 32);
                }), py::arg("size") = 0, py::arg("default_value") = 0),

        add_int_class<int_vector<64>, uint64_t>(m, "Int64Vector")
            .def(py::init([](size_t size, uint64_t default_value) {
                    return int_vector<64>(size, default_value, 64);
                }), py::arg("size") = 0, py::arg("default_value") = 0)
    );

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
        //add_bitvector_class<sdsl::bit_vector_il<64>>(m, "BitVectorIL64",
        //                                              doc_bit_vector_il),
        add_bitvector_class<sdsl::bit_vector_il<512>>(m, "BitVectorIL512",
                                                      doc_bit_vector_il)
    );

    auto rrr_classes = std::make_tuple(
        add_bitvector_class<sdsl::rrr_vector<63>>(m, "RRRVector63",
                                                  doc_rrr_vector)
    );

    auto sd_classes = std::make_tuple(
        add_bitvector_class<sdsl::sd_vector<>>(m, std::string("SDVector"),
                                               doc_sd_vector)
    );

    auto hyb_classes = std::make_tuple(
        add_bitvector_class<sdsl::hyb_vector<16>>(m, std::string("HybVector16"),
                                                  doc_hyb_vector)
    );

    add_support_class<sdsl::rank_support_v<0, 1>>(
        m, "RankSupportV_0", bit_vector_classes, "init_rank_0", doc_rank_v
    );
    add_support_class<sdsl::rank_support_v<1, 1>>(
        m, "RankSupportV_1", bit_vector_classes, "init_rank", doc_rank_v
    );
    add_support_class<sdsl::rank_support_v<10, 2>>(
        m, "RankSupportV_10",  bit_vector_classes, "init_rank_10", doc_rank_v
    );
    add_support_class<sdsl::rank_support_v<11, 2>>(
        m, "RankSupportV_11", bit_vector_classes, "init_rank_11", doc_rank_v
    );

    add_support_class<sdsl::rank_support_v5<0, 1>>(
        m, "RankSupportV5_0", bit_vector_classes, "init_rankV5_0", doc_rank_v5
    );
    add_support_class<sdsl::rank_support_v5<1, 1>>(
        m, "RankSupportV5_1", bit_vector_classes, "init_rankV5", doc_rank_v5
    );
    add_support_class<sdsl::rank_support_v5<10, 2>>(
        m, "RankSupportV5_10", bit_vector_classes, "init_rankV5_10", doc_rank_v5
    );
    add_support_class<sdsl::rank_support_v5<11, 2>>(
        m, "RankSupportV5_11", bit_vector_classes, "init_rankV5_11", doc_rank_v5
    );

    add_support_class<sdsl::rank_support_scan<0, 1>>(
        m, "RankSupportScan_0",  bit_vector_classes, "init_rankScan_0",
        doc_rank_scan
    );
    add_support_class<sdsl::rank_support_scan<1, 1>>(
        m, "RankSupportScan_1", bit_vector_classes, "init_rankScan",
        doc_rank_scan);
    add_support_class<sdsl::rank_support_scan<10, 2>>(
        m, "RankSupportScan_10", bit_vector_classes, "init_rankScan_10",
        doc_rank_scan);
    add_support_class<sdsl::rank_support_scan<11, 2>>(
        m, "RankSupportScan_11", bit_vector_classes, "init_rankScan_11",
        doc_rank_scan);

    add_support_class<sdsl::rank_support_il<0>>(m, "RankSupportIL_0",
                                                bvil_classes, "init_rank_0");
    add_support_class<sdsl::rank_support_il<1>>(m, "RankSupportIL_1",
                                                bvil_classes, "init_rank");

    add_support_class<sdsl::rank_support_rrr<0, 63>>(
        m, "RankSupportRRR63_0", rrr_classes, "init_rank_0"
    );
    add_support_class<sdsl::rank_support_rrr<1, 63>>(
        m, "RankSupportRRR63_1", rrr_classes, "init_rank"
    );

    add_support_class<sdsl::rank_support_sd<0>>(m, "RankSupportSD_0",
                                                sd_classes, "init_rank_0");
    add_support_class<sdsl::rank_support_sd<1>>(m, "RankSupportSD_1",
                                                sd_classes, "init_rank");

    add_support_class<sdsl::rank_support_hyb<0>>(m, "RankSupportHyb_0",
                                                 hyb_classes, "init_rank_0");
    add_support_class<sdsl::rank_support_hyb<1>>(m, "RankSupportHyb_1",
                                                 hyb_classes, "init_rank");

    add_support_class<support_helper<sdsl::select_support_mcl<0, 1>>>(
        m, "SelectSupportMCL_0", bit_vector_classes, "init_select_0",
        doc_select_mcl
    );
    add_support_class<support_helper<sdsl::select_support_mcl<1, 1>>>(
        m, "SelectSupportMCL_1", bit_vector_classes, "init_select",
        doc_select_mcl
    );
    add_support_class<support_helper<sdsl::select_support_mcl<10, 2>>>(
        m, "SelectSupportMCL_10", bit_vector_classes, "init_select_10",
        doc_select_mcl
    );
    add_support_class<support_helper<sdsl::select_support_mcl<11, 2>>>(
        m, "SelectSupportMCL_11", bit_vector_classes, "init_select_11",
        doc_select_mcl
    );

    add_support_class<support_helper<sdsl::select_support_scan<0, 1>>>(
        m, "SelectSupportScan_0", bit_vector_classes, "init_selectScan_0",
        doc_select_scan
    );
    add_support_class<support_helper<sdsl::select_support_scan<1, 1>>>(
        m, "SelectSupportScan_1", bit_vector_classes, "init_selectScan",
        doc_select_scan
    );
    add_support_class<support_helper<sdsl::select_support_scan<10, 2>>>(
        m, "SelectSupportScan_10", bit_vector_classes, "init_selectScan_10",
        doc_select_scan
    );
    add_support_class<support_helper<sdsl::select_support_scan<01, 2>>>(
        m, "SelectSupportScan_11", bit_vector_classes, "init_selectScan_11",
        doc_select_scan
    );

    add_support_class<sdsl::select_support_il<0>>(
        m, "SelectSupportIL_0", bvil_classes, "init_select_0"
    );
    add_support_class<sdsl::select_support_il<1>>(
        m, "SelectSupportIL_1", bvil_classes, "init_select"
    );

    add_support_class<sdsl::select_support_rrr<0, 63>>(
        m, "SelectSupportRRR63_0", rrr_classes, "init_select_0"
    );
    add_support_class<sdsl::select_support_rrr<1, 63>>(
        m, "SelectSupportRRR63_1", rrr_classes, "init_select"
    );

    add_support_class<sdsl::select_support_sd<0>>(
        m, "SelectSupportSD_0", sd_classes, "init_select_0"
    );
    add_support_class<sdsl::select_support_sd<1>>(
        m, "SelectSupportSD_1", sd_classes, "init_select"
    );

    //add_support_class<sdsl::select_support_hyb<>>(m, "SelectSupportHyb",
    //                                              hyb_classes);

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

    for_each_in_tuple(rrr_classes, make_inits_many_functor(bit_vector_classes));
    for_each_in_tuple(rrr_classes, make_inits_many_functor(rrr_classes));

    for_each_in_tuple(sd_classes, make_inits_many_functor(bit_vector_classes));
    for_each_in_tuple(sd_classes, make_inits_many_functor(sd_classes));

    for_each_in_tuple(hyb_classes, make_inits_many_functor(bit_vector_classes));
    for_each_in_tuple(hyb_classes, make_inits_many_functor(hyb_classes));

    for_each_in_tuple(iv_classes, make_pysequence_init_functor());
    for_each_in_tuple(enc_classes, make_pysequence_init_functor());
    for_each_in_tuple(vlc_classes, make_pysequence_init_functor());
    for_each_in_tuple(dac_classes, make_pysequence_init_functor());

    //for_each_in_tuple(sd_classes, make_pysequence_init_functor());
}
