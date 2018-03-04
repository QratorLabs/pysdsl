#pragma once

#include <string>
#include <vector>

#include <sdsl/wavelet_trees.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "calc.hpp"
#include "io.hpp"


namespace py = pybind11;



template <class T, bool enable = T::lex_ordered>
class add_lex_functor;

template <class T, bool enbale = sdsl::has_node_type<T>::value>
class add_traversable_functor;


template <class T>
class add_lex_functor<T, false>
{
public:
    py::class_<T>& operator() (py::class_<T>& cls)
    {
        return cls;
    }
};


template <class T>
class add_lex_functor<T, true>
{
public:
    py::class_<T>& operator() (py::class_<T>& cls)
    {
        cls.def(
            "quantile_freq",
            [] (const T& self, typename T::size_type lb,
                typename T::size_type rb, typename T::size_type q)
            {
                return sdsl::quantile_freq(self, lb, rb, q);
            },
            py::arg("lb"), py::arg("rb"), py::arg("q"),
            "Returns the q-th smallest element and its frequency in wt[lb..rb]."
            "\n\tlb: Left array bound in T"
            "\n\trb: Right array bound in T"
            "\n\tq: q-th largest element ('quantile'), 0-based indexed.",
            py::call_guard<py::gil_scoped_release>()
        );
        cls.def(
            "lex_count",
            [] (const T& self, size_t i, size_t j, typename T::value_type c)
            {
                if (j >= self.size())
                {
                    throw std::invalid_argument("j should be less than size");
                }
                if (i >= j)
                {
                    throw std::invalid_argument("i should be less than j");
                }
                return self.lex_count(i, j, c);
            },
            py::arg("i"), py::arg("j"), py::arg("c"),
            "How many values are lexicographic smaller/greater than c in "
            "[i..j-1].\n\ti: Start index (inclusive) of the interval."
            "\n\tj: End index (exclusive) of the interval."
            "\n\tc: Value c.\nreturn A triple containing:\n\trank(i, c)"
            "\n\tnumber of values smaller than c in [i..j-1]"
            "\n\tnumber of values greater than c in [i..j-1]",
            py::call_guard<py::gil_scoped_release>()
        );
        cls.def(
            "lex_smaller_count",
            [] (const T& self, size_t i, typename T::value_type c)
            {
                if (i >= self.size())
                {
                    throw std::invalid_argument("i should be less than size");
                }
                return self.lex_smaller_count(i, c);
            },
            py::arg("i"), py::arg("c"),
            "How many values are lexicographic smaller than c in [0..i-1]."
            "\n\ti: Exclusive right bound of the range."
            "\nreturn: A tuple containing:\n\trank(i, c)\n\tnumber of values "
            "smaller than c in [0..i-1]",
            py::call_guard<py::gil_scoped_release>()
        );
        cls.def(
            "symbol_lte",
            [] (const T& self, typename T::value_type c)
            {
                auto result = sdsl::symbol_lte(self, c);
                if (!std::get<0>(result))
                {
                    throw std::runtime_error("Symbol not found");
                }
                return std::get<1>(result);
            },
            py::arg("c"),
            "Returns for a symbol c the previous smaller or equal symbol in "
            "the WT"
        );
        cls.def(
            "symbol_gte",
            [] (const T& self, typename T::value_type c)
            {
                auto result = sdsl::symbol_gte(self, c);
                if (!std::get<0>(result))
                {
                    throw std::runtime_error("Symbol not found");
                }
                return std::get<1>(result);
            },
            py::arg("c"),
            "Returns for a symbol c the next larger or equal symbol in the WT"
        );
        return cls;
    }
};


template <class T>
class add_traversable_functor<T, false>
{
public:
    py::class_<T>& operator() (py::module&, py::class_<T>& cls, std::string&&)
    {
        return cls;
    }
};


template <class T>
class add_traversable_functor<T, true>
{
public:
    py::class_<T>& operator() (py::module& m, py::class_<T>& cls,
                               std::string&& name)
    {
        typedef typename T::node_type t_node;

        try
        {
            py::class_<t_node> node_cls(m, name.c_str())
    //            .def_property_readonly("sym", &t_node::sym )
            ;
        }
        catch(std::runtime_error&) {}


        cls.def("root_node", &T::root);
        cls.def("node_is_leaf", &T::is_leaf);
        cls.def(
            "node_empty",
            [] (const T& self, const t_node& node)
            { return self.empty(node); }
        );
        cls.def(
            "node_size",
            [] (const T& self, const t_node& node)
            { return self.size(node); }
        );
        cls.def("node_sym", &T::sym);
        cls.def(
            "node_expand",
            [] (const T& self, const t_node& node)
            { return self.expand(node); }
        );
        cls.def(
            "node_expand_ranges",
            [] (const T& self, const t_node& node,
                const sdsl::range_vec_type& ranges)
            {
                return self.expand(node, ranges);
            },
            py::arg("node"), py::arg("ranges")
        );
        cls.def("node_bit_vec", &T::bit_vec);
        cls.def("node_seq", &T::seq);

        cls.def(
            "intersect",
            [] (const T& self, std::vector<sdsl::range_type> ranges, size_t t)
            {
                return sdsl::intersect(self, ranges, t);
            },
            py::arg("ranges"), py::arg("t") = 0,
            "Intersection of elements in "
            "WT[s₀, e₀], WT[s₁, e₁], ...,WT[sₖ,eₖ]\n"
            "\tranges: The ranges.\n\tt: Threshold in how many distinct ranges "
            "he value has to be present. Default: t=ranges.size()\n"
            "Return a vector containing (value, frequency) - of value which "
            "are contained in t different ranges. Frequency = accumulated "
            "frequencies in all ranges. The tuples are ordered according "
            "to value, if wt is lex_ordered."
        );
        cls.def(
            "interval_symbols",
            [] (const T& self, size_t i, size_t j)
            {
                if (j > self.size())
                {
                    throw std::invalid_argument("j should be less or equal "
                                                "than size");
                }
                if (i > j)
                {
                    throw std::invalid_argument("i should be less or equal "
                                                "than j");
                }
                size_t k;
                std::vector<typename T::value_type> cs(self.sigma);
                std::vector<size_t> rank_c_i(self.sigma);
                std::vector<size_t> rank_c_j(self.sigma);

                sdsl::interval_symbols(self, i, j, k, cs, rank_c_i, rank_c_j);

                return std::make_tuple(k, cs, rank_c_i, rank_c_j);
            },
            py::arg("i"), py::arg("j"),
            "For each symbol c in wt[i..j - 1] get rank(i, c) and rank(j, c)."
        );
        return cls;
    }
};



template <class T>
inline
auto add_wavelet_class(py::module& m, const std::string&& name,
                       const char* doc= nullptr)
{
    auto cls = py::class_<T>(m, name.c_str())
        .def_property_readonly("sigma", [] (const T& self) {
            return self.sigma;
        }, "Effective alphabet size of the wavelet tree")
        .def("get_sigma", [] (const T& self) {
            return self.sigma;
        }, "Effective alphabet size of the wavelet tree")
        // .def_property_readonly("tree", [] (const T& self) {
        //     return self.tree;
        // }, "A concatenation of all bit vectors of the wavelet tree.")

        // .def_property_readonly("max_level", [] (const T& self) {
        //     return self.max_level;
        // }, "Maximal level of the wavelet tree.")
        .def_static(
            "from_bytes",
            [] (const py::bytes& bytes)
            {
                T wt;
                sdsl::construct_im(wt, std::string(bytes), 1);
                return wt;
            },
            py::arg("s"),
            "Construct from a build sequence"
        )
        .def_static(
            "parse_string",
            [] (const std::string& s)
            {
                T wt;
                sdsl::construct_im(wt, s, 'd');
                return wt;
            },
            py::arg("s"),
            "Construct from space-separated human-readable string"
        )
        .def(
            "rank",
            [] (const T& self, typename T::size_type i,
                typename T::value_type c)
            {
                if (i >= self.size())
                {
                    throw std::out_of_range(std::to_string(i));
                }
                return self.rank(i, c);
            },
            "Calculates how many values c are in the prefix [0..i-1] of the "
            "supported vector (i in [0..size]).\nTime complexity: "
            "Order(log(|Sigma|))",
            py::arg("i"), py::arg("c"),
            py::call_guard<py::gil_scoped_release>()
        )
        .def(
            "inverse_select",
            [] (const T& self, typename T::size_type i)
            {
                if (i >= self.size())
                {
                    throw std::out_of_range(std::to_string(i));
                }
                return self.inverse_select(i);
            },
            py::arg("i"),
            "Calculates how many occurrences of value wt[i] are in the prefix"
            "[0..i-1] of the original sequence, returns pair "
            "(rank(wt[i], i), wt[i])",
            py::call_guard<py::gil_scoped_release>()
        )
        .def(
            "select",
            [] (const T& self, typename T::size_type i,
                typename T::value_type c)
            {
                if (i < 1 || i >= self.size())
                {
                    throw std::out_of_range(std::to_string(i));
                }
                if (i > self.rank(self.size(), c))
                {
                    throw std::invalid_argument(
                        std::to_string(i) + " is greater than rank(" +
                        std::to_string(i) + ", " + std::to_string(c) + ")"
                    );
                }
                return self.select(i, c);
            },
            py::arg("i"), py::arg("c"),
            "Calculates the i-th occurrence of the value c in the supported "
            "vector.\nTime complexity: Order(log(|Sigma|))",
            py::call_guard<py::gil_scoped_release>()
        )
    ;

    add_lex_functor<T>()(cls);
    add_traversable_functor<T>()(m, cls, name + "Node");

    //add_sizes(cls);
    add_description(cls);
    add_serialization(cls);
    add_to_string(cls);

    add_std_algo(cls);

    if (doc) cls.doc() = doc;

     return cls;

}
