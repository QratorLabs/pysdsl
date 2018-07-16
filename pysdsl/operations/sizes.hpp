#pragma once

#include <pybind11/pybind11.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/io.hpp>


namespace py = pybind11;


namespace detail
{
    struct no_max_size {};
    struct has_max_size: no_max_size {};
    template<typename> struct int_ { typedef int type; };

    template <class T, typename /* unused */ = decltype(T::max_size)>
    inline
    auto add_max_size_impl(py::class_<T>& cls, has_max_size /* unused */)
    {
        cls.def_property_readonly_static(
            "max_size",
            [](py::object /* self */) { return T::max_size(); },
            "Maximum size of the int_vector.");
        return cls;
    }

    template <class T>
    constexpr
    auto add_max_size_impl(py::class_<T>& cls, no_max_size /* unused */)
    { return cls; }

    template <class T> constexpr auto size(const T& seq) { return seq.size(); }
}  // namespace detail


template <class T>
constexpr auto add_max_size(py::class_<T>& cls)
{
    return detail::add_max_size_impl(cls, detail::has_max_size());
}


template <uint32_t B>
inline auto add_max_size(py::class_<sdsl::bit_vector_il<B>>& cls)
{
    cls.def_property_readonly_static(
        "max_size",
        [](py::object /* self */) {
            return sdsl::bit_vector::max_size();},
        "Maximum size of the bit_vector_il.");
    return cls;
}


template <class Sequence>
inline auto add_sizes(py::class_<Sequence>& cls)
{
    add_max_size(cls);

    auto size = [] (const Sequence& self) {
        return self.size(); };

    cls.def("__len__", size,
            "The number of elements in the container.");
    cls.def_property_readonly("size", size,
                              "The number of elements in the container.");
    return cls;
}
