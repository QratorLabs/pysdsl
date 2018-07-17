#pragma once

#include <algorithm>
#include <tuple>

#include <pybind11/pybind11.h>

#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "types/pysequence.hpp"
#include "operations/creation.hpp"
#include "operations/iteration.hpp"
#include "util/tupletricks.hpp"


namespace py = pybind11;

namespace {
    using detail::cbegin;
    using detail::cend;
} // namespace


namespace detail
{
    template <class BindCls,
              typename py_class = py::class_<typename BindCls::type>>
    class add_init_functor
    {
    public:
        add_init_functor(BindCls &cls_to_add_def) : m_cls_to(cls_to_add_def) {}

        template <typename InputCls>
        decltype(auto) operator()(const InputCls &)
        {
            m_cls_to.def(py::init(
                [] (const typename InputCls::type& from) {
                    return construct_from<typename BindCls::type>(from); }),
                py::arg("v"),
                py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

    private:
        BindCls& m_cls_to;
    };


    template <class... From>
    class add_many_inits_to_each
    {
    public:
        add_many_inits_to_each(const std::tuple<From...>& from_each):
                               m_from_each(from_each) {}

        template <typename BindCls>
        decltype(auto) operator()(BindCls& cls) {
            return for_each(m_from_each, add_init_functor<BindCls>(cls)); }

    private:
        const std::tuple<From...>& m_from_each;
    };


    template <class T>
    class pysequence_init_functor
    {
        typedef typename T::value_type value_type;
    public:
        decltype(auto) operator()(py::class_<T>& cls) {
            return cls.def(py::init(
                [] (const py::sequence& v) {
                    return construct_from<T>(
                        sequence_wrapper<value_type>(v)); }),
                py::arg("v")); }
    };


    class add_pysequence_init_functor
    {
    public:
        add_pysequence_init_functor() {}

        template <class BindCls,
                  class py_class = py::class_<typename BindCls::type>>
        decltype(auto) operator()(BindCls& cls) {
            return pysequence_init_functor<typename BindCls::type>()(cls); }
    };
}


inline
auto make_pysequence_init_functor() {
    return detail::add_pysequence_init_functor(); }


template <class... From>
inline
auto make_inits_many_functor(const std::tuple<From...>& from_each) {
    return detail::add_many_inits_to_each<From...>(from_each); }
