#pragma once

#include <algorithm>
#include <tuple>

#include <pybind11/pybind11.h>

#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "pysequence.hpp"
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
              typename Tag = typename BindCls::type::index_category,
              typename py_class = py::class_<typename BindCls::type>>
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
    class add_init_functor<py::class_<sdsl::int_vector<B>>, void>
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
                    std::copy(cbegin(from), cend(from), std::begin(result));
                    return result;
                }),
            py::arg("v"), py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

    private:
        BindCls& m_cls_to;
    };

    template <class BindCls>
    class add_init_functor<BindCls, sdsl::bv_tag>
    {
    public:
        add_init_functor(BindCls &cls_to_add_def) : m_cls_to(cls_to_add_def) {}

        auto operator()(const py::class_<sdsl::int_vector<1>>&)
        {
            m_cls_to.def(py::init([](const sdsl::bit_vector& from) {
                return typename BindCls::type(from);
            }), py::arg("v"), py::call_guard<py::gil_scoped_release>());
            return m_cls_to;
        }

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
                    sdsl::bit_vector temp(from.size());
                    std::copy(cbegin(from), cend(from), std::begin(temp));
                    return typename BindCls::type(temp);
                }),
                py::arg("v"), py::call_guard<py::gil_scoped_release>(),
                "\tInvolves intermediate BitVector creation"
            );
            return m_cls_to;
        }
    private:
        BindCls& m_cls_to;
    };

    template <class T>
    class add_init_functor<py::class_<T>, sdsl::wt_tag>
    {
    public:
        typedef py::class_<T> BindCls;

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
                    std::copy(cbegin(from), cend(from), std::begin(temp));
                    sdsl::construct_im(result, temp);
                    return result; }),
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
                    return result; }),
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

        template <uint8_t B>
        decltype(auto) operator()(py::class_<sdsl::int_vector<B>>& cls) {
            return for_each(m_from_each,
                            add_init_functor<py::class_<sdsl::int_vector<B>>,
                                             void>(cls)); }

        template <class BV, int ML>
        decltype(auto) operator()(py::class_<sdsl::dac_vector_dp<BV, ML>>& cls)
        {
            return for_each(
                m_from_each,
                add_init_functor<py::class_<sdsl::dac_vector_dp<BV, ML>>,
                                 void>(cls));
        }

        template <uint32_t K>
        decltype(auto) operator()(py::class_<sdsl::hyb_vector<K>>& cls) {
            return for_each(m_from_each,
                            add_init_functor<py::class_<sdsl::hyb_vector<K>>,
                                             sdsl::bv_tag>(cls)); }

    private:
        const std::tuple<From...>& m_from_each;
    };


    template <class T, typename Tag = typename T::index_category>
    class pysequence_init_functor
    {
        typedef typename T::value_type value_type;
    public:
        decltype(auto) operator()(py::class_<T>& cls) {
            return cls.def(py::init([](const py::sequence& v) {
                return T(sequence_wrapper<value_type>(v)); }), py::arg("v")); }
    };

    template <uint8_t width>
    class pysequence_init_functor<sdsl::int_vector<width>, void>
    {
        typedef py::class_<sdsl::int_vector<width>> BindCls;
        typedef typename BindCls::type base_class;
        typedef typename base_class::value_type value_type;

    public:
        decltype(auto) operator()(BindCls& cls)
        {
            return cls.def(py::init(
                [](const py::sequence& v) {
                    const auto vsize = v.size();
                    base_class result(vsize);
                    sequence_wrapper<value_type> helper(v);

                    std::copy(helper.begin(), helper.end(), result.begin());
                    return result; }),
                py::arg("v"));
        }
    };

    template <class T>
    class pysequence_init_functor<T, sdsl::wt_tag>
    {
        typedef py::class_<T> BindCls;
        typedef typename T::value_type value_type;
    public:
        decltype(auto) operator()(BindCls &cls)
        {
            return cls.def(py::init(
                [](const py::sequence &v) {
                    sequence_wrapper<value_type> helper(v);

                    T result;
                    sdsl::int_vector<> temp(v.size());
                    std::copy(helper.begin(), helper.end(), temp.begin());

                    sdsl::construct_im(result, temp);

                    return result; }),
                py::arg("v"), "\tInvolves intermediate IntVector creation");
        }
    };

    template <class T>
    class pysequence_init_functor<T, sdsl::csa_tag>
    {
        typedef py::class_<T> BindCls;
        typedef typename T::value_type value_type;
    public:
        decltype(auto) operator()(BindCls &cls)
        {
            cls.def(py::init(
                [](const py::sequence &v) {
                    sequence_wrapper<value_type> helper(v);

                    T result;
                    sdsl::int_vector<> temp(v.size());
                    std::copy(helper.begin(), helper.end(), temp.begin());

                    sdsl::construct_im(result, temp);

                    return result; }),
                py::arg("v"), "\tInvolves intermediate IntVector creation");
            cls.def_static(
                "construct_from_file",
                [] (const std::string& file_name) {
                    T result;
                    sdsl::construct(result, file_name);
                    return result; },
                py::arg("file_name"));
            return cls;
        }
    };

    class add_pysequence_init_functor
    {
    public:
        add_pysequence_init_functor() {}

        template <class BindCls,
                  class py_class = py::class_<typename BindCls::type>>
        decltype(auto) operator()(BindCls& cls) {
            return pysequence_init_functor<typename BindCls::type>()(cls); }

        template <uint8_t B>
        decltype(auto) operator()(py::class_<sdsl::int_vector<B>>& cls) {
            return pysequence_init_functor<sdsl::int_vector<B>, void>()(cls); }

        template <class BV, int ML>
        decltype(auto) operator()(py::class_<sdsl::dac_vector_dp<BV, ML>>& cls)
        {
            return pysequence_init_functor<sdsl::dac_vector_dp<BV, ML>,
                                           void>()(cls);
        }

        template <uint32_t K>
        decltype(auto) operator()(py::class_<sdsl::hyb_vector<K>>& cls) {
            return pysequence_init_functor<sdsl::hyb_vector<K>,
                                           sdsl::bv_tag>()(cls); }

        // template <class Th, class Ts1, class Ts0>
        // decltype(auto) operator()(
        //     py::class_<sdsl::sd_vector<Th, Ts1, Ts0>>& cls
        // )
        // {
        //     typedef typename py::class_<sdsl::sd_vector<Th, Ts1, Ts0>> py_class;
        //     typedef typename py_class::type base_class;
        //     typedef typename base_class::value_type value_type;

        //     cls.def(py::init([](const py::sequence& v) {
        //         sequence_wrapper<value_type> helper(v);

        //         return base_class(helper.begin(), helper.end());
        //     }), py::arg("v"));
        //     return cls;
        // }

    };
}


inline
auto make_pysequence_init_functor() {
    return detail::add_pysequence_init_functor(); }


template <class... From>
inline
auto make_inits_many_functor(const std::tuple<From...>& from_each) {
    return detail::add_many_inits_to_each<From...>(from_each); }
