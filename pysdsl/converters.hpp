#pragma once


#include <tuple>

#include <pybind11/pybind11.h>

namespace py = pybind11;


namespace detail
{
    template<int... Is>
    struct seq { };

    template<int N, int... Is>
    struct gen_seq : gen_seq<N - 1, N - 1, Is...> { };

    template<int... Is>
    struct gen_seq<0, Is...> : seq<Is...> { };

    template<typename T, typename F, int... Is>
    void for_each(T&& t, F f, seq<Is...>)
    {
        auto l = { (f(std::get<Is>(t)), 0)... };
    }

    template <class ToCls>
    class add_inits_functor
    {
    public:
        add_inits_functor(ToCls &cls_to) : m_cls_to(cls_to) {}

        template <typename FromCls>
        void operator()(FromCls &t)
        {
            m_cls_to.def(py::init([](const typename FromCls::type& from) {
                //py::print("Fast!");
                return typename ToCls::type(from);
            }), py::arg("v"), py::call_guard<py::gil_scoped_release>());
        }

    private:
        ToCls& m_cls_to;
    };
}


template <typename... Ts, typename F>
void for_each_in_tuple(std::tuple<Ts...> const &t, F f)
{
    detail::for_each(t, f, detail::gen_seq<sizeof...(Ts)>());
}


template <class ToCls>
auto make_inits_functor(ToCls &cls_to)
{
    return detail::add_inits_functor<ToCls>(cls_to);
}
