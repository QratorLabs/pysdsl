#pragma once

#include "operations/sizes.hpp"
#include "util/indexiterator.hpp"


namespace detail
{
    struct no_iterator {};
    struct has_iterator: no_iterator {};

    template <class Sequence, typename /* unused */ = decltype(Sequence::begin)>
    constexpr
    auto cbegin_impl(const Sequence& sequence, has_iterator /* unused */)
    {
        return sequence.begin();
    }

    template <class Sequence, typename /* unused */ = decltype(Sequence::end)>
    constexpr
    auto cend_impl(const Sequence& sequence, has_iterator /* unused */)
    {
        return sequence.end();
    }

    template <class Sequence>
    constexpr
    auto cbegin_impl(const Sequence& sequence, no_iterator /* unused */)
    {
        return count_index_iterator<Sequence>(&sequence, 0);
    }

    template <class Sequence>
    constexpr
    auto cend_impl(const Sequence& sequence, no_iterator /* unused */)
    {
        return count_index_iterator<Sequence>(&sequence, size(sequence));
    }

    template <class Sequence>
    constexpr auto cbegin(const Sequence& sequence)
    {
        return cbegin_impl(sequence, has_iterator());
    }

    template <class Sequence>
    constexpr auto cend(const Sequence& sequence)
    {
        return cend_impl(sequence, has_iterator());
    }

}  // namespace detail


template <class Sequence>
inline auto add_iteration(py::class_<Sequence>& cls)
{
    return cls.def(
        "__iter__",
        [](const Sequence &sequence) {
            return py::make_iterator(detail::cbegin(sequence),
                                     detail::cend(sequence)); },
        py::keep_alive<0, 1>()
    );
}
