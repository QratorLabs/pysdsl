/*cppimport
<%
cfg['compiler_args'] = ['-std=c++14', '-fvisibility=hidden']
cfg['linker_args'] = ['-fvisibility=hidden']
cfg['include_dirs'] = ['sdsl-lite/include']
cfg['libraries'] = ['sdsl', 'divsufsort', 'divsufsort64']
cfg['dependencies'] = ['converters.hpp', 'pysequence.hpp',
                       'sizes.hpp', 'calc.hpp']
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
        .def_property_readonly("capacity", &T::capacity,
                               "Returns the size of the occupied bits of the "
                               "int_vector. The capacity of a int_vector is "
                               "greater or equal to the bit_size of the "
                               "vector: capacity >= bit_size).")

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
             "Set all entries of int_vector to value k. This method "
             "pre-calculates the content of at most 64 words and then "
             "repeatedly inserts these words.",
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
    );

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
        add_int_class<int_vector<0>>(m, "IntVector",
                                     "This generic vector class could be used "
                                     "to generate a vector that contains "
                                     "integers of fixed width `w` in [1..64].")
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
                "Bit compress the int_vector. Determine the biggest value X "
                "and then set the int_width to the smallest possible so that "
                "we still can represent X.",
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
        add_compressed_class<sdsl::dac_vector<>>(
        m, "DacVector",
        "A generic immutable space-saving vector class for unsigned integers.\n"
        "The values of a dac_vector are immutable after the constructor call.\n"
        "The `escaping` technique is used to encode values.\n"
        "This is defined as follows (see [1]):\n"
        "A k-bit integer is split into `K=k/(b-1)` bits each and "
        "encoded into `K` blocks of `b` bits each. All but the last block "
        "are marked with by a 1 in the most significant bit. Escaping with "
        "`b=8` is also known as vbyte-coding (see [2]). A experimental study "
        "of using escaping for the LCP array is given in [3].\n"
        "Time complexity: Order{log n/b} worst case, where b is the number "
        "of bits in a block\nReferences:\n"
        "[1] F. Transier and P. Sanders: `Engineering Basic Search Algorithms "
        "of an In-Memory Text Search Engine`, ACM Transactions on "
        "Information Systems, Vol. 29, No.1, Article 2, 2010\n"
        "[2] H.E. Williams and J. Zobel: `Compressing integers for fast file "
        "access`, Computing Journal Vol 43, No.3, 1999\n"
        "[3] N. Brisboa, S. Ladra, G. Navarro: `Directly addressable "
        "variable-length codes'', Proceedings of SPIRE 2009."
    )
    .def_property_readonly("levels", &sdsl::dac_vector<>::levels),
    add_compressed_class<sdsl::dac_vector_dp<>>(
        m, "DacVectorDP",
        "A generic immutable space-saving vector class for unsigned integers.\n"
        "The values of a dac_vector are immutable after the constructor call.\n"
        "The \"escaping\" technique is used to encode values. Bit widths of "
        "each encoding level are chosen optimally via dynamic programming.\n"
        "References\n [1] N. Brisaboa and S. Ladra and G. Navarro: `DACs: "
        "Bringing Direct Access to Variable-Length Codes`, "
        "Information Processing and Management (IPM) 2013"
    )
    .def("cost", &sdsl::dac_vector_dp<>::cost, py::arg("n"), py::arg("m"))
    .def_property_readonly("levels", &sdsl::dac_vector_dp<>::levels)
    );

    auto bvil_classes = std::make_tuple(
        add_bitvector_class<sdsl::bit_vector_il<512>>(
            m, "BitVectorIL512",
            "A bit vector which interleaves the original bit_vector with rank "
            "information. \nThis class is a uncompressed bit vector "
            "representation. It copies the original bit_vector and interleaves "
            "the data every t_bs bits with a cumulative sum of set bits before "
            "the current position. Each cumulative sum is stored in a 64 bit "
            "word."
        )
    );

    auto rrr_classes = std::make_tuple(
        add_bitvector_class<sdsl::rrr_vector<63>>(
            m, "RRRVector63",
            "An H_0-compressed bitvector representation.\n"
            "References:\n"
            "- Rasmus Pagh, Low redundancy in dictionaries with O(1) worst "
            "case lookup time, Technical Report 1998. "
            "ftp://ftp.cs.au.dk/BRICS/Reports/RS/98/28/BRICS-RS-98-28.pdf, "
            "Section 2.\n"
            "- Rajeev Raman, V. Raman and S. Srinivasa Rao, Succinct Indexable "
            "Dictionaries with Applications to representations of k-ary trees "
            "and multi-sets. SODA 2002.\n"
            "- Francisco Claude, Gonzalo Navarro: Practical Rank/Select "
            "Queries over Arbitrary Sequences. SPIRE 2008: 176-187\n"
            "- On the fly-decoding and encoding was discovered in; Gonzalo "
            "Navarro, Eliana Providel: Fast, Small, Simple Rank/Select on "
            "Bitmaps. SEA 2012"
        )
    );

    auto sd_classes = std::make_tuple(
        add_bitvector_class<sdsl::sd_vector<>>(
            m, std::string("SDVector"),
            "A bit vector which compresses very sparse populated bit vectors "
            "by representing the positions of 1 by the Elias-Fano "
            "representation for non-decreasing sequences\n"
            "References:\n"
            "- P. Elias: ''Efficient storage and retrieval by content and "
            "address of static files'', Journal of the ACM, 1974\n"
            "- R. Fano: ''On the number of bits required to implement an "
            "associative memory'', Memorandum 61. Computer Structures Group, "
            "Project MAC, MIT, 1971\n"
            "- D. Okanohara, K. Sadakane: ''Practical Entropy-Compressed "
            "Rank/Select Dictionary'', Proceedings of ALENEX 2007."
        )
    );

    auto hyb_classes = std::make_tuple(
        add_bitvector_class<sdsl::hyb_vector<16>>(
            m, std::string("HybVector16"),
            "A hybrid-encoded compressed bitvector representation\n"
            "References:\n- Juha Karkkainen, Dominik Kempa and "
            "Simon J. Puglisi. ''Hybrid Compression of Bitvectors for the "
            "FM-Index.'' DCC 2014."
        )
    );

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
