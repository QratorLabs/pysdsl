#pragma once

#include <algorithm>
#include <tuple>

#include <pybind11/pybind11.h>

#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "pysequence.hpp"


namespace py = pybind11;


namespace detail
{
    template<typename P, typename Function, std::size_t... Is>
    inline
    decltype(auto) for_each_impl(P&& t, Function&& f, std::index_sequence<Is...>)
    {
        return std::make_tuple(f(std::get<Is>(t))...);
    }

    template<typename P, typename Function, std::size_t... Is>
    inline
    decltype(auto) for_each_impl(P& t, Function&& f, std::index_sequence<Is...>)
    {
        return std::make_tuple(f(std::get<Is>(t))...);
    }

    template<typename... T, typename Function>
    inline
    decltype(auto) for_each(const std::tuple<T...>& t, Function&& f)
    {
        return for_each_impl(t, f, std::index_sequence_for<T...>{});
    }

    template<typename... T, typename Function>
    inline
    decltype(auto) for_each(std::tuple<T...>& t, Function&& f)
    {
        return for_each_impl(t, f, std::index_sequence_for<T...>{});
    }

    template <class BindCls>
    class add_init_functor
    {
    public:
        add_init_functor(BindCls &cls_to_add_def) : m_cls_to(cls_to_add_def) {}

        template <typename InputCls>
        decltype(auto) operator()(const InputCls &)
        {
            m_cls_to.def(py::init([](const typename InputCls::type& from) {
                return typename BindCls::type(from);
            }), py::arg("v"), py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

    private:
        BindCls& m_cls_to;
    };

    template <uint8_t B>
    class add_init_functor<py::class_<sdsl::int_vector<B>>>
    {
    public:
        typedef py::class_<sdsl::int_vector<B>> BindCls;

        add_init_functor(BindCls &cls_to_add_def) : m_cls_to(cls_to_add_def) {}

        auto operator()(const BindCls&)
        {
            m_cls_to.def(py::init([](const typename BindCls::type& from) {
                return typename BindCls::type(from); /*just a copy*/
            }), py::arg("v"), py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

        template <typename InputCls>
        decltype(auto) operator()(const InputCls &)
        {
            m_cls_to.def(
                py::init([](const typename InputCls::type& from) {
                    typename BindCls::type result(from.size());
                    std::copy(from.begin(), from.end(), result.begin());
                    return result;
                }),
            py::arg("v"), py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

    private:
        BindCls& m_cls_to;
    };

    #define CREATE_INIT_FUNCTOR(...)\
    class add_init_functor<__VA_ARGS__>\
    {\
    public:\
        typedef __VA_ARGS__ BindCls;\
        add_init_functor(BindCls &cls_to_add_def) : m_cls_to(cls_to_add_def) {}\
        auto operator()(const py::class_<sdsl::int_vector<1>>&)\
        {\
            m_cls_to.def(py::init([](const sdsl::bit_vector& from) {\
                return typename BindCls::type(from);\
            }), py::arg("v"), py::call_guard<py::gil_scoped_release>());\
            return m_cls_to;\
        }\
        auto operator()(const BindCls&)\
        {\
            m_cls_to.def(py::init([](const typename BindCls::type& from) {\
                return typename BindCls::type(from); /*just a copy*/ \
            }), py::arg("v"), py::call_guard<py::gil_scoped_release>());\
            return m_cls_to;\
        }\
        template <typename InputCls>\
        decltype(auto) operator()(const InputCls &)\
        {\
            m_cls_to.def(\
                py::init([](const typename InputCls::type& from) {\
                    sdsl::bit_vector temp(from.size());\
                    std::copy(from.begin(), from.end(), temp.begin());\
                    return typename BindCls::type(temp);\
                }),\
                py::arg("v"), py::call_guard<py::gil_scoped_release>(),\
                "\tInvolves intermediate BitVector creation"\
            );\
            return m_cls_to;\
        }\
    private:\
        BindCls& m_cls_to;\
    };

    template <uint32_t B>
    CREATE_INIT_FUNCTOR(py::class_<sdsl::bit_vector_il<B>>);

    template <uint16_t B, class R, uint16_t K>
    CREATE_INIT_FUNCTOR(py::class_<sdsl::rrr_vector<B, R, K>>);

    template <class T, class Ts0, class Ts1>
    CREATE_INIT_FUNCTOR(py::class_<sdsl::sd_vector<T, Ts0, Ts1>>);

    template <uint32_t K>
    CREATE_INIT_FUNCTOR(py::class_<sdsl::hyb_vector<K>>);

    template <class B, class R, class S0, class S1>
    class add_init_functor<py::class_<sdsl::wt_int<B, R, S0, S1>>>
    {
    public:
        typedef py::class_<sdsl::wt_int<B, R, S0, S1>> BindCls;

        add_init_functor(BindCls &cls_to_add_def) : m_cls_to(cls_to_add_def) {}

        auto operator()(const BindCls&)
        {
            m_cls_to.def(py::init([](const typename BindCls::type& from) {
                return typename BindCls::type(from); /*just a copy*/
            }), py::arg("v"), py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

        template <typename InputCls>
        decltype(auto) operator()(const InputCls &)
        {
            m_cls_to.def(
                py::init([](const typename InputCls::type& from) {
                    typename BindCls::type result;
                    sdsl::int_vector<> temp(from.size());
                    std::copy(from.begin(), from.end(), temp.begin());

                    sdsl::construct_im(result, temp);

                    return result;
                }),
                "\tInvolves intermediate IntVector creation",
                py::arg("v"), py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

        template <uint8_t width>
        decltype(auto) operator()(const py::class_<sdsl::int_vector<width>> &)
        {
            m_cls_to.def(
                py::init([](const sdsl::int_vector<width>& from) {
                    typename BindCls::type result;

                    sdsl::construct_im(result, from);

                    return result;
                }),
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
        decltype(auto) operator()(BindCls &cls)
        {
            return for_each(m_from_each, add_init_functor<BindCls>(cls));
        }

    private:
        const std::tuple<From...>& m_from_each;
    };

    class add_pysequence_init_functor
    {
    public:
        add_pysequence_init_functor() {}

        template <class BindCls>
        decltype(auto) operator()(BindCls& cls)
        {
            typedef typename BindCls::type base_class;
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
                sequence_wrapper<value_type> helper(v);

                std::copy(helper.begin(), helper.end(), result.begin());
                return result;
            }), py::arg("v"));
            return cls;
        }

        template <class Th, class Ts1, class Ts0>
        decltype(auto) operator()(
            py::class_<sdsl::sd_vector<Th, Ts1, Ts0>>& cls
        )
        {
            typedef typename py::class_<sdsl::sd_vector<Th, Ts1, Ts0>> py_class;
            typedef typename py_class::type base_class;
            typedef typename base_class::value_type value_type;

            cls.def(py::init([](const py::sequence& v) {
                sequence_wrapper<value_type> helper(v);

                return base_class(helper.begin(), helper.end());
            }), py::arg("v"));
            return cls;
        }

        template <class B, class Tr, class Ts1, class Ts0>
        decltype(auto) operator()(
            py::class_<sdsl::wt_int<B, Tr, Ts1, Ts0>>& cls
        )
        {
            typedef typename py::class_<sdsl::wt_int<B, Tr, Ts1, Ts0>> py_class;
            typedef typename py_class::type base_class;
            typedef typename base_class::value_type value_type;

            cls.def(py::init([](const py::sequence& v) {
                sequence_wrapper<value_type> helper(v);

                base_class result;
                sdsl::int_vector<> temp(v.size());
                std::copy(helper.begin(), helper.end(), temp.begin());

                sdsl::construct_im(result, temp);

                return result;

            }), py::arg("v"), "\tInvolves intermediate IntVector creation");
            return cls;
        }

    };
}


template <typename... Ts, typename F>
inline
decltype(auto) for_each_in_tuple(const std::tuple<Ts...> &t, F f)
{
    return detail::for_each(t, f);
}


template <typename... Ts, typename F>
inline
decltype(auto) for_each_in_tuple(std::tuple<Ts...> &t, F f)
{
    return detail::for_each(t, f);
}


inline
auto make_pysequence_init_functor()
{
    return detail::add_pysequence_init_functor();
}


template <class... From>
inline
auto make_inits_many_functor(const std::tuple<From...>& from_each)
{
    return detail::add_many_inits_to_each<From...>(from_each);
}
