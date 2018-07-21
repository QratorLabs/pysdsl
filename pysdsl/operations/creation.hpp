#pragma once

#include <algorithm>
#include <utility>
#include <tuple>
#include <type_traits>

#include <sdsl/vectors.hpp>
#include <sdsl/construct.hpp>

#include <pybind11/pybind11.h>

#include "operations/iteration.hpp"
#include "operations/sizes.hpp"
#include "types/pysequence.hpp"
#include "util/tupletricks.hpp"


namespace py = pybind11;


namespace detail
{

template <class T, typename value_type = typename T::value_type>
struct IntermediateVector { using type = sdsl::int_vector<>; };


template <class T>
struct IntermediateVector<T, bool> { using type = sdsl::int_vector<1>; };


template <
    typename T,
    typename = typename std::enable_if<
            !std::is_same<typename T::iterator,
                          typename T::const_iterator>::value
        >::type,
    typename = decltype(std::declval<T>().begin())>
std::true_type has_non_const_begin_impl(T *);

std::false_type has_non_const_begin_impl(...);


template <class T>
using has_non_const_begin = decltype(has_non_const_begin_impl(
        std::declval<T*>()));


struct construct_explicit { int value; };
struct construct_iter { int value; };
struct construct_copy_empty { int value; };
struct construct_copy_size { int value; };

}  // namespace detail


// The only version of construct_from if T can be constructed from From
template <class T, class From,
    typename /* direct construction */ = typename std::enable_if<
            std::is_constructible<T, const From&>::value
        >::type>
constexpr T construct_from(const From& obj,
                           detail::construct_explicit /* unused */ = {}) {
    return T(obj); }


// The only version of construct_from if T can be constructed from iterator
template <
    class T, class From,
    typename /* no direct construction */ = typename std::enable_if<
            !std::is_constructible<T, const From&>::value
        >::type,
    typename /* construct from iterator */ = typename std::enable_if<
            std::is_constructible<
                T,
                decltype(detail::cbegin(std::declval<From>())),
                decltype(detail::cend(std::declval<From>()))
            >::value
        >::type>
constexpr T construct_from(const From& obj,
                           detail::construct_iter /* unused */ = {}) {
    return T(detail::cbegin(obj), detail::cend(obj)); }


// The only version of construct_from if T can be filled from iterator and
// can only be created empty
template <
    class T, class From,
    typename /* no direct construction */ = typename std::enable_if<
            !std::is_constructible<T, const From&>::value
        >::type,
    typename /* no construction from iterator */ = typename std::enable_if<
            !std::is_constructible<
                T,
                decltype(detail::cbegin(std::declval<From>())),
                decltype(detail::cend(std::declval<From>()))>::value
        >::type,
    typename /* default constructable */ = typename std::enable_if<
            std::is_constructible<T>::value
        >::type,
    typename /* no construction with size */ = typename std::enable_if<
            !std::is_constructible<T, std::size_t>::value
        >::type,
    typename /* can be modified after construction */ = typename std::enable_if<
            detail::has_non_const_begin<T>::value
        >::type>
inline T construct_from(const From& obj,
                        detail::construct_copy_empty /* unused */ = {})
{
    T result;
    std::copy(detail::cbegin(obj), detail::cend(obj), result.begin());
    return result;
}


// The only version of construct_from if T can be filled from iterator and
// can allocate memory aforehand
template <
    class T, class From,
    typename /* no direct construction */ = typename std::enable_if<
            !std::is_constructible<T, const From&>::value
        >::type,
    typename /* no construction from iterator */ = typename std::enable_if<
            !std::is_constructible<
                T,
                decltype(detail::cbegin(std::declval<From>())),
                decltype(detail::cend(std::declval<From>()))>::value
        >::type,
    typename /* construction with known size */ = typename std::enable_if<
            std::is_constructible<
                T, decltype(detail::size(std::declval<From>()))>::value
        >::type,
    typename /* can be modified after construction */ = typename std::enable_if<
            detail::has_non_const_begin<T>::value
        >::type>
inline T construct_from(const From& obj,
                        detail::construct_copy_size /* unused */ = {})
{
    T result(detail::size(obj));
    std::copy(detail::cbegin(obj), detail::cend(obj), result.begin());
    return result;
}

// The only version of construct_from if T can only be constructed
// via sdsl::construct_im
template <
    class T, class From,
    class With = typename detail::IntermediateVector<T>::type,
    typename /* unused */ = typename std::enable_if<
            std::is_constructible<T>::value
        >::type,
    typename = typename std::enable_if<
            !std::is_constructible<T, const From&>::value
        >::type,
    typename = typename std::enable_if<
            !std::is_constructible<
                T,
                decltype(detail::cbegin(std::declval<From>())),
                decltype(detail::cend(std::declval<From>()))>::value
        >::type,
    typename = typename std::enable_if<
            !std::is_constructible<T, std::size_t>::value
        >::type,
    typename = typename std::enable_if<
            !detail::has_non_const_begin<T>::value
        >::type>
inline T construct_from(const From& obj)
{
    T result;
    sdsl::construct_im(result, construct_from<With>(obj));

    return result;
}


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
