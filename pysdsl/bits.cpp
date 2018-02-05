/*cppimport
<%
cfg['compiler_args'] = ['-std=c++14', '-ftemplate-depth=3600']
cfg['include_dirs'] = ['sdsl-lite/include']
cfg['libraries'] = ['sdsl', 'divsufsort', 'divsufsort64']
%>
*/

#include <string>
#include <tuple>

#include "sdsl/bits.hpp"
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace detail
{
    // Convert array into a tuple
    template<typename Array, std::size_t... I>
    decltype(auto) a2t_impl(const Array& a, std::index_sequence<I...>)
    {
        return std::make_tuple(a[I]...);
    }
}


template<typename T, std::size_t N,
         typename Indices = std::make_index_sequence<N>>
decltype(auto) as_tuple(const T (&a) [N])
{
    return detail::a2t_impl(a, Indices{});
}


auto cnt11 = [] (uint64_t x, uint64_t& c) {
    auto result = sdsl::bits::cnt11(x, c);
    return std::make_pair(result, c);
};
auto cnt10 = [] (uint64_t x, uint64_t& c) {
    auto result = sdsl::bits::cnt10(x, c);
    return std::make_pair(result, c);
};
auto cnt01 = [] (uint64_t x, uint64_t& c) {
    auto result = sdsl::bits::cnt01(x, c);
    return std::make_pair(result, c);
};
auto sel = [] (uint64_t x, uint32_t i) {
    if (i >= sizeof(sdsl::bits::ps_overflow) / sizeof(sdsl::bits::ps_overflow[0])) {
        throw py::index_error(std::to_string(i));
    }
    return sdsl::bits::sel(x, i);
};
auto sel11 = [] (uint64_t x, uint32_t i, uint32_t c=0) {
    if (i >= sizeof(sdsl::bits::ps_overflow) / sizeof(sdsl::bits::ps_overflow[0])) {
        throw py::index_error(std::to_string(i));
    }
    return sdsl::bits::sel11(x, i, c);
};


PYBIND11_MODULE(bits, m) {
    m.doc() = "bitwise tricks on 64 bit words.";


    auto bits_cls = py::class_<sdsl::bits>(m, "bits")
        .def_property_readonly_static("all_set", [](py::object /* self */) { return sdsl::bits::all_set; })
        .def_static("all_set_", []() { return sdsl::bits::all_set; },
                    "64bit mask with all bits set to 1.")

        .def_property_readonly_static("deBruijn64", [](py::object) { return sdsl::bits::deBruijn64; })
        .def_static("deBruijn64_", []() { return sdsl::bits::deBruijn64; },
                    "This constant represents a de Bruijn sequence B(k,n) for k=2 and n=6. "
                    "Details for de Bruijn sequences see "
                    "http://en.wikipedia.org/wiki/De_bruijn_sequence "
                    "deBruijn64 is used in combination with the "
                    "array lt_deBruijn_to_idx.")

        .def_property_readonly_static("lt_deBruijn_to_idx", [](py::object) { return as_tuple(sdsl::bits::lt_deBruijn_to_idx); })
        .def_static("lt_deBruijn_to_idx_", []() { return as_tuple(sdsl::bits::lt_deBruijn_to_idx); },
                    "This table maps a 6-bit subsequence S[idx...idx+5] of constant deBruijn64 to idx.")

        .def_property_readonly_static("lt_fib", [](py::object) { return as_tuple(sdsl::bits::lt_fib); })
        .def_static("lt_fib_", []() { return as_tuple(sdsl::bits::lt_fib); },
                    "Array containing Fibonacci numbers less than 2**64")

        .def_property_readonly_static("lt_cnt", [](py::object) {
            return py::bytes(reinterpret_cast<const char*>(sdsl::bits::lt_cnt), sizeof(sdsl::bits::lt_cnt));
        })
        .def_static("lt_cnt_", []() {
            return py::bytes(reinterpret_cast<const char*>(sdsl::bits::lt_cnt), sizeof(sdsl::bits::lt_cnt));
        }, "Lookup table for byte popcounts.")

        .def_property_readonly_static("lt_hi", [](py::object) { return as_tuple(sdsl::bits::lt_hi); })
        .def_static("lt_hi_", []() { return as_tuple(sdsl::bits::lt_hi); },
                    "Lookup table for most significant set bit in a byte.")

        .def_property_readonly_static("lo_set", [](py::object) { return as_tuple(sdsl::bits::lo_set); })
        .def_static("lo_set_", []() { return as_tuple(sdsl::bits::lo_set); },
                    "lo_set[i] is a 64-bit word with the i least significant bits set and the high bits not set. "
                    "lo_set[0] = 0ULL, lo_set[1]=1ULL, lo_set[2]=3ULL...")

        .def_property_readonly_static("lo_unset", [](py::object) { return as_tuple(sdsl::bits::lo_unset); })
        .def_static("lo_unset_", []() { return as_tuple(sdsl::bits::lo_unset); },
                    "lo_unset[i] is a 64-bit word with the i least significant bits not set and the high bits set. "
                    "lo_unset[0] = FFFFFFFFFFFFFFFFULL, lo_unset_set[1]=FFFFFFFFFFFFFFFEULL, ...")

        .def_property_readonly_static("lt_lo", [](py::object) {
            return py::bytes(reinterpret_cast<const char*>(sdsl::bits::lt_lo), sizeof(sdsl::bits::lt_lo)); })
        .def_static("lt_lo_", []() {
            return py::bytes(reinterpret_cast<const char*>(sdsl::bits::lt_lo), sizeof(sdsl::bits::lt_lo));
        }, "Lookup table for least significant set bit in a byte.")

        .def_property_readonly_static("lt_sel", [](py::object) {
            return py::bytes(reinterpret_cast<const char*>(sdsl::bits::lt_sel), sizeof(sdsl::bits::lt_sel));
        })
        .def_static("lt_sel_", []() {
            return py::bytes(reinterpret_cast<const char*>(sdsl::bits::lt_sel), sizeof(sdsl::bits::lt_sel));
         }, "Lookup table for select on bytes. "
            "Entry at idx = 256*j + i equals the position of the "
            "(j+1)-th set bit in byte i. Positions lie in the range [0..7].")

        .def_property_readonly_static("ps_overflow", [](py::object) { return as_tuple(sdsl::bits::ps_overflow); })
        .def_static("ps_overflow_", []() { return as_tuple(sdsl::bits::ps_overflow); },
                    "Use to help to decide if a prefix sum stored in a byte overflows.")

        .def_static("cnt", &sdsl::bits::cnt, "Counts the number of set bits in x.", py::arg("x"))
        .def_static("cnt32", &sdsl::bits::cnt32, "Counts the number of set bits in 32-bit integer x.", py::arg("x"))
        .def_static("hi", &sdsl::bits::hi, "The position (in 0..63) of the most significant set bit "
                                           "in `x` or 0 if x equals 0.", py::arg("x"))
        .def_static("lo", &sdsl::bits::lo, "The position (in 0..63) of the rightmost 1-bit in the 64bit integer x if "
                                           "x>0 and 0 if x equals 0.", py::arg("x"))

        .def_static("cnt11", (uint32_t (*) (uint64_t)) &sdsl::bits::cnt11, py::arg("x"),
                    "Count the number of consecutive and distinct 11 in the 64bit integer x.\n"
                    "x: 64bit integer to count the terminating sequence 11 of a Fibonacci code.")
        .def_static("cnt11", cnt11, "Count the number of consecutive and distinct 11 in the 64bit integer x.\n"
                                    "x: 64bit integer to count the terminating sequence 11 of a Fibonacci code.\n"
                                    "c: Carry equals msb of the previous 64bit integer.", py::arg("x"), py::arg("c"))
        .def_static("cnt10", cnt10, "Count 10 bit pairs in the word x.\n"
                                    "x: 64bit integer to count the 10 bit pairs.\n"
                                    "c: Carry equals msb of the previous 64bit integer.", py::arg("x"), py::arg("c"))
        .def_static("cnt01", cnt01, "Count 01 bit pairs in the word x.\n"
                                    "x: 64bit integer to count the 01 bit pairs.\n"
                                    "c: Carry equals msb of the previous 64bit integer.", py::arg("x"), py::arg("c"))

        .def_static("map10", &sdsl::bits::map10, py::arg("x"), py::arg("c") = 0,
                    "Map all 10 bit pairs to 01 or 1 if c=1 and the lsb=0. All other pairs are mapped to 00.")
        .def_static("map01", &sdsl::bits::map01, py::arg("x"), py::arg("c") = 1,
                    "Map all 01 bit pairs to 01 or 1 if c=1 and the lsb=0. All other pairs are mapped to 00.")

        .def_static("sel",  sel, py::arg("x"), py::arg("i"),
                    "Calculate the position of the i-th rightmost 1 bit in the 64bit integer x\n"
                    "x: 64bit integer.\ni: Argument i must be in the range [1..cnt(x)].")
        .def_static("sel11", sel11, py::arg("x"), py::arg("i"), py::arg("c") = 0,
                     "The position (in 1..63) of the i-th 11-bit-pattern which terminates a Fibonacci coded integer in x if "
                     "x contains at least i 11-bit-patterns and a undefined value otherwise.\n"
                     "x: 64 bit integer.\n"
                     "i: Index of 11-bit-pattern. i in [1..cnt11(x)]\n"
                     "c: Carry bit from word before")
        .def_static("hi11", &sdsl::bits::hi11, py::arg("x"),
                    "The position (in 1..63) of the leftmost 1 of the leftmost 11-bit-pattern which "
                    "terminates a Fibonacci coded integer in x if x contains a 11-bit-pattern and 0 otherwise.")
        .def_static("rev", &sdsl::bits::rev, "reverses a given 64 bit word")
        ;

    m.attr("all_set") = bits_cls.attr("all_set_");
    m.attr("deBruijn64") = bits_cls.attr("deBruijn64_");
    m.attr("lt_deBruijn_to_idx") = bits_cls.attr("lt_deBruijn_to_idx_");
    m.attr("lt_fib") = bits_cls.attr("lt_fib_");
    m.attr("lt_cnt") = bits_cls.attr("lt_cnt_");
    m.attr("lt_hi") = bits_cls.attr("lt_hi_");
    m.attr("lo_set") = bits_cls.attr("lo_set_");
    m.attr("lo_unset") = bits_cls.attr("lo_unset_");
    m.attr("lt_lo") = bits_cls.attr("lt_lo_");
    m.attr("lt_sel") = bits_cls.attr("lt_sel_");
    m.attr("ps_overflow") = bits_cls.attr("ps_overflow_");

    m.attr("cnt") = bits_cls.attr("cnt");
    m.attr("cnt32") = bits_cls.attr("cnt32");
    m.attr("hi") = bits_cls.attr("hi");
    m.attr("lo") = bits_cls.attr("lo");
    m.attr("cnt11") = bits_cls.attr("cnt11");
    m.attr("cnt10") = bits_cls.attr("cnt10");
    m.attr("cnt01") = bits_cls.attr("cnt01");
    m.attr("map10") = bits_cls.attr("map10");
    m.attr("map01") = bits_cls.attr("map01");
    m.attr("sel") = bits_cls.attr("sel");
    m.attr("sel11") = bits_cls.attr("sel11");
    m.attr("hi11") = bits_cls.attr("hi11");
    m.attr("rev") = bits_cls.attr("rev");
}
