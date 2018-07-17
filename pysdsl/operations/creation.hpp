#pragma once

#include <algorithm>
#include <utility>
#include <type_traits>

#include <sdsl/vectors.hpp>

#include "operations/iteration.hpp"
#include "operations/sizes.hpp"


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
