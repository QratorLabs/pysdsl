#pragma once

#include <string>
#include <tuple>
#include <utility>

#include <pybind11/pybind11.h>

#include <sdsl/rmq_support_sparse_table.hpp>

#include "operations/sizes.hpp"
#include "operations/iteration.hpp"
#include "util/tupletricks.hpp"
#include "docstrings.hpp"
#include "io.hpp"
#include "calc.hpp"



struct add_rmq_sparse_table_functor {
    py::module& m;
    const char* doc;

    constexpr add_rmq_sparse_table_functor(py::module& m, const char* doc = nullptr)
        : m(m), doc(doc) {}

    template <typename t_rac, const char* rac_name, bool t_min>
    decltype(auto) operator()(std::tuple<t_rac, 
                                std::integral_constant<const char*, rac_name>, 
                                std::integral_constant<bool, t_min>>) {
        using T = sdsl::rmq_support_sparse_table<t_rac, t_min>;

        std::string name = 
            std::string("Range") + (t_min ? "Min" : "Max") + "QuerySparseTable_for_" + rac_name;

        auto cls = py::class_<T>(m, name.c_str())
            .def_property_readonly("size", (typename T::size_type(T::*)(void) const)& T::size)
            .def(py::init([](const t_rac* rac) {return T(rac);}))
            .def("__call__", (typename T::size_type(T::*)(typename T::size_type, typename T::size_type) const)& T::operator());
    
        add_sizes(cls);
        add_description(cls);

        // load takes two params
        // add_serialization(cls);

        return cls;
    }
};


namespace RAC_names {
    const char INT_VECTOR_NAME[] = "IntVector";
}

template <typename t_rac, typename t_min_ic>
using general_rmq_sparse_table = py::class_<sdsl::rmq_support_sparse_table<t_rac, t_min_ic::value>>;


inline auto add_rmq_classes(py::module& m) {
    m.attr("rmq_sparse_table") = py::dict();


    using rmq_support_sparse_table_params = std::tuple<
        std::tuple<sdsl::int_vector<>,
                   std::integral_constant<const char*, RAC_names::INT_VECTOR_NAME>,
                   std::integral_constant<bool, true>>,
        std::tuple<sdsl::int_vector<>, 
                   std::integral_constant<const char*, RAC_names::INT_VECTOR_NAME>,
                   std::integral_constant<bool, false>>
    >;

    auto rmq_sparse_tables = for_each_in_tuple(rmq_support_sparse_table_params(), 
                                add_rmq_sparse_table_functor(m, doc_rmq_sparse_table));

    ////////////// as params //////////////////////
    using rmq_support_sparse_table_as_params = std::tuple<
        std::tuple<sdsl::int_vector<>, std::integral_constant<bool, true>>,
        std::tuple<sdsl::int_vector<>, std::integral_constant<bool, false>>
    >;

    auto rmq_sparse_tables_as_params = for_each_in_tuple(rmq_support_sparse_table_as_params(), 
                                make_general_sybset_functor<general_rmq_sparse_table>(rmq_sparse_tables));
    ///////////////////////////////////////////////
    
    return std::tuple_cat(rmq_sparse_tables);
}
