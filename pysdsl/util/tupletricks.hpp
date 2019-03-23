#pragma once

#include <tuple>
#include <utility>


namespace detail
{

template<typename P, typename Function, std::size_t... Is>
constexpr
decltype(auto) for_each_impl(P&& t, Function&& f, std::index_sequence<Is...>) {
    return std::make_tuple(f(std::get<Is>(t))...); }

template<typename P, typename Function, std::size_t... Is>
constexpr
decltype(auto) for_each_impl(P& t, Function&& f, std::index_sequence<Is...>) {
    return std::make_tuple(f(std::get<Is>(t))...); }

template<typename... T, typename Function>
constexpr
decltype(auto) for_each(const std::tuple<T...>& t, Function&& f) {
    return for_each_impl(t, f, std::index_sequence_for<T...>{}); }

template<typename... T, typename Function>
constexpr
decltype(auto) for_each(std::tuple<T...>& t, Function&& f) {
    return for_each_impl(t, f, std::index_sequence_for<T...>{}); }

template<typename P, typename Function, std::size_t... Is>
constexpr
decltype(auto) forward_each_impl(P&& t, Function&& f, std::index_sequence<Is...>) {
    return std::forward_as_tuple(f(std::get<Is>(t))...); }

template<typename P, typename Function, std::size_t... Is>
constexpr
decltype(auto) forward_each_impl(P& t, Function&& f, std::index_sequence<Is...>) {
    return std::forward_as_tuple(f(std::get<Is>(t))...); }

template<typename... T, typename Function>
constexpr
decltype(auto) forward_each(const std::tuple<T...>& t, Function&& f) {
    return forward_each_impl(t, f, std::index_sequence_for<T...>{}); }

template<typename... T, typename Function>
constexpr
decltype(auto) forward_each(std::tuple<T...>& t, Function&& f) {
    return forward_each_impl(t, f, std::index_sequence_for<T...>{}); }

}  // namespace detail


template <typename... Ts, typename F>
constexpr
decltype(auto) for_each_in_tuple(const std::tuple<Ts...> &t, F f) {
    return detail::for_each(t, f); }


template <typename... Ts, typename F>
constexpr
decltype(auto) for_each_in_tuple(std::tuple<Ts...> &t, F f) {
    return detail::for_each(t, f); }

template <typename... Ts, typename F>
constexpr
decltype(auto) forward_each_in_tuple(const std::tuple<Ts...> &t, F f) {
    return detail::forward_each(t, f); }


template <typename... Ts, typename F>
constexpr
decltype(auto) forward_each_in_tuple(std::tuple<Ts...> &t, F f) {
    return detail::forward_each(t, f); }
