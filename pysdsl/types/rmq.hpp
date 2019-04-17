#pragma once

#include <string>
#include <tuple>
#include <utility>

#include <pybind11/pybind11.h>

#include <sdsl/rmq_support_sparse_table.hpp>
#include <sdsl/rmq_succinct_sada.hpp>

#include "operations/sizes.hpp"
#include "operations/iteration.hpp"
#include "util/tupletricks.hpp"
#include "docstrings.hpp"
#include "io.hpp"
#include "calc.hpp"



namespace detail {
// adds constructors of t_rac... containers
template <typename PybindClass, typename... t_rac> 
typename std::enable_if<sizeof...(t_rac) == 0>::type add_rac_constructor(const PybindClass&) {}

template <typename PybindClass, typename t_rac_head, typename... t_rac_tail>
void add_rac_constructor(PybindClass& cls)
{
    cls.def(py::init([](const t_rac_head* rac) {
        return typename PybindClass::type(rac);
    }));
    add_rac_constructor<PybindClass, t_rac_tail...>(cls);
}

} // namespace details


// containers names
namespace RAC_names {
    constexpr char INT_VECTOR_NAME[] = "IntVector",
                   INT16_VECTOR_NAME[] = "Int16Vector";
}


struct add_rmq_sparse_table_functor {
    py::module& m;
    const char* doc;

    constexpr add_rmq_sparse_table_functor(py::module& m, const char* doc = nullptr)
        : m(m), doc(doc) {}

    template <typename t_rac, const char* rac_name, bool t_min>
    decltype(auto) operator()(std::tuple<t_rac, 
                                         std::integral_constant<const char*, rac_name>, 
                                         std::integral_constant<bool, t_min>>)
    {
        using Table = sdsl::rmq_support_sparse_table<t_rac, t_min>;
        using size_type = typename Table::size_type;

        std::string name = 
            std::string("Range") + (t_min ? "Min" : "Max") + "QuerySparseTable_for_" + rac_name;

        auto cls = py::class_<Table>(m, name.c_str())
            .def(py::init([](const t_rac* rac) {return Table(rac);}))
            .def("set_vector", &Table::set_vector,
                 "Sets a vector rmq is processed on.")
            .def("__call__",
                (size_type (Table::*)(size_type, size_type) const)& Table::operator(),
                 (std::string("Returns an index of the ") + (t_min ? "minimal" : "maximal") +
                             " value on the segment [l,r].").c_str());
    
        add_sizes(cls);
        add_description(cls);

        // load takes two params
        // add_serialization(cls);

        cls.doc() = doc;

        std::string key = std::string(t_min ? "Min" : "Max") + "_in_" + rac_name;

        m.attr("rmq_sparse_tables").attr("__setitem__")(key, cls);
        m.attr("all_rmq_classes").attr("append")(cls);

        return cls;
    }
};


struct add_rmq_sada_functor {
    py::module& m;
    const char* doc;

    constexpr add_rmq_sada_functor(py::module& m, const char* doc = nullptr)
        : m(m), doc(doc) {}


    template <typename... t_rac, bool t_min>
    decltype(auto) operator()(std::tuple<std::tuple<t_rac...>,
                                         std::integral_constant<bool, t_min>>)
    {
        using RMQClass = typename std::conditional<t_min, sdsl::rmq_succinct_sada<>,
                                        typename sdsl::range_maximum_support_sada<>::type>::type;
        using size_type = typename RMQClass::size_type;

        std::string name =
            std::string("Range") + (t_min ? "Min" : "Max") + "QuerySuccintSada";

        auto cls = py::class_<RMQClass>(m, name.c_str())
            .def(py::init())
            .def("__call__",
                (size_type (RMQClass::*)(size_type, size_type) const)& RMQClass::operator(),
                 (std::string("Returns an index of the ") + (t_min ? "minimal" : "maximal") +
                             " value on the segment [l,r].").c_str());;

        detail::add_rac_constructor<decltype(cls), t_rac...>(cls);

        add_sizes(cls);
        add_description(cls);
        add_serialization(cls);

        cls.doc() = doc;

        m.attr("rmq_sada").attr("__setitem__")(t_min ? "Min" : "Max", cls);
        m.attr("all_rmq_classes").attr("append")(cls);

        return cls;
    }
};


struct add_rmq_sct_functor {
    py::module& m;
    const char* doc;

    constexpr add_rmq_sct_functor(py::module& m, const char* doc = nullptr)
        : m(m), doc(doc) {}


    template <typename... t_rac, bool t_min>
    decltype(auto) operator()(std::tuple<std::tuple<t_rac...>,
                                         std::integral_constant<bool, t_min>>)
    {
        using RMQClass = typename std::conditional<t_min, sdsl::rmq_succinct_sct<>,
                                                typename sdsl::range_maximum_sct<>::type>::type;
        using size_type = typename RMQClass::size_type;

        std::string name =
            std::string("Range") + (t_min ? "Min" : "Max") + "QuerySuccintSct";

        auto cls = py::class_<RMQClass>(m, name.c_str())
            .def(py::init())
            .def("__call__",
                (size_type (RMQClass::*)(size_type, size_type) const)& RMQClass::operator(),
                 (std::string("Returns an index of the ") + (t_min ? "minimal" : "maximal") +
                             " value on the segment [l,r].").c_str());


        detail::add_rac_constructor<decltype(cls), t_rac...>(cls);

        add_sizes(cls);
        add_description(cls);
        add_serialization(cls);

        cls.doc() = doc;

        m.attr("rmq_sct").attr("__setitem__")(t_min ? "Min" : "Max", cls);
        m.attr("all_rmq_classes").attr("append")(cls);

        return cls;
    }
};


// generalized (constants -> typenames) template for usage with GeneralSubsetFunctor
template <typename t_rac, typename t_min_integral_constant>
using general_rmq_sparse_table =
    py::class_<sdsl::rmq_support_sparse_table<t_rac, t_min_integral_constant::value>>;

template <typename t_min_integral_constant>
using general_rmq_sada = py::class_<
    typename std::conditional<t_min_integral_constant::value, 
                              sdsl::rmq_succinct_sada<>, 
                              typename sdsl::range_maximum_support_sada<>::type>::type>;

template <typename t_min_integral_constant>
using general_rmq_sct = py::class_<
    typename std::conditional<t_min_integral_constant::value,
                              sdsl::rmq_succinct_sct<>,
                              typename sdsl::range_maximum_sct<>::type>::type>;


inline auto add_rmq_classes(py::module& m)
{
    m.attr("rmq_sparse_tables") = py::dict();
    m.attr("rmq_sada") = py::dict();
    m.attr("rmq_sct") = py::dict();
    m.attr("all_rmq_classes") = py::list();

    using rmq_support_sparse_table_params = std::tuple<
        std::tuple<sdsl::int_vector<>,
                   std::integral_constant<const char*, RAC_names::INT_VECTOR_NAME>,
                   std::integral_constant<bool, true>>,
        std::tuple<sdsl::int_vector<>, 
                   std::integral_constant<const char*, RAC_names::INT_VECTOR_NAME>,
                   std::integral_constant<bool, false>>,
        std::tuple<sdsl::int_vector<16>, 
                   std::integral_constant<const char*, RAC_names::INT16_VECTOR_NAME>,
                   std::integral_constant<bool, true>>,
        std::tuple<sdsl::int_vector<16>, 
                   std::integral_constant<const char*, RAC_names::INT16_VECTOR_NAME>,
                   std::integral_constant<bool, false>>
    >;

    using rmq_sada_params = std::tuple<
        std::tuple<std::tuple<sdsl::int_vector<>, sdsl::int_vector<16>>,
                   std::integral_constant<bool, true>>,
        std::tuple<std::tuple<sdsl::int_vector<>, sdsl::int_vector<16>>,
                   std::integral_constant<bool, false>>
    >;

    using rmq_sct_params = std::tuple<
        std::tuple<std::tuple<sdsl::int_vector<>, sdsl::int_vector<16>>,
                   std::integral_constant<bool, true>>,
        std::tuple<std::tuple<sdsl::int_vector<>, sdsl::int_vector<16>>,
                   std::integral_constant<bool, false>>
    >;

    auto rmq_sparse_tables = for_each_in_tuple(rmq_support_sparse_table_params(), 
                                add_rmq_sparse_table_functor(m, doc_rmq_sparse_table));
    auto rmq_sada_classes = for_each_in_tuple(rmq_sada_params(),
                                add_rmq_sada_functor(m, doc_rmq_sada));
    auto rmq_sct_classes = for_each_in_tuple(rmq_sct_params(),
                                add_rmq_sct_functor(m, doc_rmq_sct));
    
    return std::make_tuple(rmq_sparse_tables, rmq_sada_classes, rmq_sct_classes);
}
