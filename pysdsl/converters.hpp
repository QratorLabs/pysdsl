#pragma once


#include <tuple>

#include <pybind11/pybind11.h>

namespace py = pybind11;


namespace detail
{
    template<typename P, typename Function, std::size_t... Is>
    decltype(auto) for_each_impl(P&& t, Function f, std::index_sequence<Is...>)
    {
        return std::make_tuple(f(std::get<Is>(t))...);
    }

    template<typename... T, typename Function>
    decltype(auto) for_each(const std::tuple<T...>& t, Function f)
    {
        return for_each_impl(t, f, std::index_sequence_for<T...>{});
    }

    template <class ToCls>
    class add_inits_functor
    {
    public:
        add_inits_functor(ToCls &cls_to) : m_cls_to(cls_to) {}

        template <typename FromCls>
        decltype(auto) operator()(FromCls &t)
        {
            m_cls_to.def(py::init([](const typename FromCls::type& from) {
                //py::print("Fast!");
                return typename ToCls::type(from);
            }), py::arg("v"), py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

    private:
        ToCls& m_cls_to;
    };
}


template <typename... Ts, typename F>
decltype(auto) for_each_in_tuple(const std::tuple<Ts...> &t, F f)
{
    return detail::for_each(t, f);
}


template <class ToCls>
auto make_inits_functor(ToCls &cls_to)
{
    return detail::add_inits_functor<ToCls>(cls_to);
}
