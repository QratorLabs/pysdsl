/*cppimport
<%
cfg['compiler_args'] = ['-std=c++14', '-fvisibility=hidden']
cfg['linker_args'] = ['-fvisibility=hidden']
cfg['include_dirs'] = ['sdsl-lite/include']
cfg['libraries'] = ['sdsl', 'divsufsort', 'divsufsort64']
cfg['dependencies'] = ['converters.hpp', 'pysequence.hpp']
%>
*/

#include <algorithm>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include <sdsl/vectors.hpp>

#include <pybind11/pybind11.h>

#include "converters.hpp"


namespace py = pybind11;


using sdsl::int_vector;


template <class Sequence>
auto add_std_algo(py::class_<Sequence>& cls)
{
    cls.def(
        "__iter__",
        [](const Sequence &s) { return py::make_iterator(s.begin(), s.end()); },
        py::keep_alive<0, 1>()
    );
    cls.def(
        "__contains__",
        [](const Sequence &self, typename Sequence::value_type element) {
            return std::find(self.begin(), self.end(), element) != self.end();
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "__getitem__",
        [](const Sequence &self, size_t position) {
            if (position >= self.size())
            {
                throw py::index_error(std::to_string(position));
            }
            return self[position];
        }
    );
    cls.def(
        "__getitem__",
        [](const Sequence &self, int64_t position) {
            auto abs_position = std::abs(position);
            if (position >= 0)
            {
                throw std::exception();
            }
            if (abs_position > self.size())
            {
                throw py::index_error(std::to_string(position));
            }
            return self[self.size() - abs_position];
        }
    );
    cls.def(
        "max",
        [](const Sequence &self) {
            return *std::max_element(self.begin(), self.end());
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "min",
        [](const Sequence &self) {
            return *std::min_element(self.begin(), self.end());
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "minmax",
        [](const Sequence &self) {
            auto result = std::minmax_element(self.begin(), self.end());
            return std::make_pair(*std::get<0>(result), *std::get<1>(result));
        },
        py::call_guard<py::gil_scoped_release>()
    );
    cls.def(
        "sum",
        [](const Sequence &self) {
            return std::accumulate(self.begin(), self.end(),
                                   uint64_t(0));
        },
        py::call_guard<py::gil_scoped_release>()
    );

    return cls;
}


template <class T>
auto add_sizes(py::class_<T>& cls)
{
    cls.def("__len__", &T::size, "The number of elements in the int_vector.");
    cls.def_property_readonly("size", &T::size,
                              "The number of elements in the int_vector.");
    cls.def_property_readonly_static(
        "max_size",
        [](py::object /* self */) { return T::max_size(); },
        "Maximum size of the int_vector."
    );
    cls.def_property_readonly(
        "size_in_mega_bytes",
        [](const T &self) { return sdsl::size_in_mega_bytes(self); }
    );
    return cls;
}



template <class T>
auto add_io(py::class_<T>& cls)
{
    cls.def(
        "__str__",
        [](const T &self) { return "[" + sdsl::util::to_string(self) + "]";}
    );
    cls.def("to_latex",
            [](const T &self) { return sdsl::util::to_latex_string(self); });

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
                    throw py::index_error(std::to_string(position));
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
             py::call_guard<py::gil_scoped_release>(),
             py::arg("k"),
             "Set all entries of int_vector to value k. This method "
             "pre-calculates the content of at most 64 words and then "
             "repeatedly inserts these words."
        )
        .def("set_zero_bits",
             [](T &self) { sdsl::util::_set_zero_bits(self); },
             py::call_guard<py::gil_scoped_release>(),
             "Sets all bits of the int_vector to 0-bits.")
        .def("set_one_bits",
             [](T &self) { sdsl::util::_set_one_bits(self); },
             py::call_guard<py::gil_scoped_release>(),
             "Sets all bits of the int_vector to 1-bits.")
        .def(
            "set_random_bits",
            [](T &self, int seed) {
                sdsl::util::set_random_bits(self, seed);
            },
            py::call_guard<py::gil_scoped_release>(),
            py::arg_v(
                "seed",
                0,
                "If seed = 0, the time is used to initialize the pseudo "
                "random number generator, otherwise the seed parameter is used."
            ),
            "Sets all bits of the int_vector to pseudo-random bits."
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
            "Number of set bits in vector")
        .def("cnt_onezero_bits",
             [](const T &self) { return sdsl::util::cnt_onezero_bits(self); },
             "Number of occurrences of bit pattern `10` in vector")
        .def("cnt_zeroone_bits",
             [](const T &self) { return sdsl::util::cnt_zeroone_bits(self); },
             "Number of occurrences of bit pattern `01` in vector")

        .def(
            "next_bit",
            [](const T &self, size_t idx) {
                if (idx >= self.bit_size())
                {
                    throw py::index_error(std::to_string(idx));
                }
                return sdsl::util::next_bit(self, idx);
            },
            py::arg("idx"),
            "Get the smallest position `i` >= `idx` where a bit is set"
        )
        .def(
            "prev_bit",
            [](const T &self, size_t idx) {
                if (idx >= self.bit_size())
                {
                    throw py::index_error(std::to_string(idx));
                }
                return sdsl::util::prev_bit(self, idx);
            },
            py::arg("idx"),
            "Get the largest position `i` <= `idx` where a bit is set"
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
                    throw py::index_error(std::to_string(i));
                }
                return self.sample(i);
            },
            "Returns the i-th sample of the compressed vector"
            "i: The index of the sample. 0 <= i < size()/get_sample_dens()"
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
                py::arg("bit_width") = 64)
            .def(
                "expand_width",
                [](int_vector<0> &self, size_t width) {
                    sdsl::util::expand_width(self, width);
                },
                "Expands the integer width to new_width >= v.width()."
            )
            .def("bit_compress",
                [](int_vector<0> &self) { sdsl::util::bit_compress(self); },
                "Bit compress the int_vector. Determine the biggest value X "
                "and then set the int_width to the smallest possible so that "
                "we still can represent X."),

        add_int_class<int_vector<1>, bool>(m, "BitVector")
            .def(py::init([](size_t size, bool default_value) {
                    return int_vector<1>(size, default_value, 1);
                }), py::arg("size") = 0, py::arg("default_value") = false)
            .def("flip", &int_vector<1>::flip, "Flip all bits of bit_vector"),

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
    .def_property_readonly("levels", &sdsl::dac_vector<>::levels)
// add_compressed_class<sdsl::dac_vector_dp<>>(
//     m, "DacVectorDP", iv_classes,
//     "A generic immutable space-saving vector class for unsigned integers.\n"
//     "The values of a dac_vector are immutable after the constructor call.\n"
//     "The ,,escaping'' technique is used to encode values. Bit widths of "
//     "each encoding level are chosen optimally via dynamic programming.\n"
//     "References\n [1] N. Brisaboa and S. Ladra and G. Navarro: `DACs: "
//     "Bringing Direct Access to Variable-Length Codes`, "
//     "Information Processing and Management (IPM) 2013"
// )
// .def("cost", &sdsl::dac_vector_dp<>::cost)
// .def_property_readonly("levels", &sdsl::dac_vector_dp<>::levels)
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

    for_each_in_tuple(iv_classes, make_pysequence_init_functor());
    for_each_in_tuple(enc_classes, make_pysequence_init_functor());
    for_each_in_tuple(vlc_classes, make_pysequence_init_functor());
    for_each_in_tuple(dac_classes, make_pysequence_init_functor());

}
