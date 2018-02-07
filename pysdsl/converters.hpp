#pragma once


#include <tuple>

#include <pybind11/pybind11.h>
#include <sdsl/vectors.hpp>

#include "pysequence.hpp"


namespace py = pybind11;


namespace detail
{
    template<typename P, typename Function, std::size_t... Is>
    decltype(auto) for_each_impl(P&& t, Function f, std::index_sequence<Is...>)
    {
        return std::make_tuple(f(std::get<Is>(t))...);
    }

    template<typename P, typename Function, std::size_t... Is>
    decltype(auto) for_each_impl(P& t, Function f, std::index_sequence<Is...>)
    {
        return std::make_tuple(f(std::get<Is>(t))...);
    }

    template<typename... T, typename Function>
    decltype(auto) for_each(const std::tuple<T...>& t, Function f)
    {
        return for_each_impl(t, f, std::index_sequence_for<T...>{});
    }

    template<typename... T, typename Function>
    decltype(auto) for_each(std::tuple<T...>& t, Function f)
    {
        return for_each_impl(t, f, std::index_sequence_for<T...>{});
    }

    template <class ToCls>
    class add_inits_functor
    {
    public:
        add_inits_functor(ToCls &cls_to) : m_cls_to(cls_to) {}

        template <typename FromCls>
        decltype(auto) operator()(FromCls &)
        {
            m_cls_to.def(py::init([](const typename FromCls::type& from) {
                return typename ToCls::type(from);
            }), py::arg("v"), py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

    private:
        ToCls& m_cls_to;
    };

    template <class... From>
    class add_inits_many_functor
    {
    public:
        add_inits_many_functor(const std::tuple<From...>& from_each):
                               m_from_each(from_each) {}

        template <typename ToCls>
        decltype(auto) operator()(ToCls &cls)
        {
            return for_each(m_from_each, add_inits_functor<ToCls>(cls));
        }

    private:
        const std::tuple<From...>& m_from_each;
    };

    class add_pysequence_init_functor
    {
    public:
        add_pysequence_init_functor() {}

        template <class ToCls>
        decltype(auto) operator()(ToCls& cls)
        {
            typedef typename ToCls::type base_class;
            typedef typename base_class::value_type value_type;

            cls.def(py::init([](const py::sequence& v) {
                return base_class(sequence_wrapper<value_type>(v));
            }), py::arg("v"));
            return cls;
        }

        template <uint8_t width>
        decltype(auto) operator()(py::class_<sdsl::int_vector<width>> &cls)
        {
            typedef typename py::class_<sdsl::int_vector<width>> py_class;
            typedef typename py_class::type base_class;
            typedef typename base_class::value_type value_type;

            cls.def(py::init([](const py::sequence& v) {
                const auto vsize = v.size();
                base_class result(vsize);
                for (size_t i = 0; i < vsize; i++)
                {
                    result[i] = py::cast<value_type>(v[i]);
                }
                return result;
            }), py::arg("v"));
            return cls;
        }
    };
}


template <typename... Ts, typename F>
decltype(auto) for_each_in_tuple(const std::tuple<Ts...> &t, F f)
{
    return detail::for_each(t, f);
}


template <typename... Ts, typename F>
decltype(auto) for_each_in_tuple(std::tuple<Ts...> &t, F f)
{
    return detail::for_each(t, f);
}


template <class ToCls>
auto make_inits_functor(ToCls &cls_to)
{
    return detail::add_inits_functor<ToCls>(cls_to);
}


auto make_pysequence_init_functor()
{
    return detail::add_pysequence_init_functor();
}


template <class ToCls>
auto add_pysequence_init(ToCls &cls)
{
    return make_pysequence_init_functor()(cls);
}


template <class... From>
auto make_inits_many_functor(const std::tuple<From...>& from_each)
{
    return detail::add_inits_many_functor<From...>(from_each);
}
