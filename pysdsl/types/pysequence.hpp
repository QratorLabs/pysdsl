#pragma once

#include <utility>

#include <pybind11/pybind11.h>

namespace py = pybind11;


namespace detail
{
    template <class T, class Base>
    class sequence_iterator_wrapper
    {
    public:
        using difference_type = std::size_t;
        using value_type = T;
        using pointer = T*;
        using reference = T&;
        using iterator_category = std::random_access_iterator_tag;

        sequence_iterator_wrapper(Base it): m_it(it) {}

        bool operator!=(const sequence_iterator_wrapper& other) const {
            return m_it != other.m_it; }

        bool operator==(const sequence_iterator_wrapper& other) const {
            return m_it == other.m_it; }

        value_type operator*() { return py::cast<T>(*m_it); }

        decltype(auto) operator++() {
            ++m_it;
            return *this;
        }

        decltype(auto) operator++(int) {
            return *sequence_iterator_wrapper<T, Base>(m_it++); }

        decltype(auto) operator-(difference_type step) const {
            return sequence_iterator_wrapper<T, Base>(m_it - step); }

        difference_type operator-(sequence_iterator_wrapper other) const {
            return m_it - other.m_it; }

        decltype(auto) operator=(sequence_iterator_wrapper other)
        {
            if (this != &other) {
                m_it = other.m_it; }
            return *this;
        }

    private:
        Base m_it;
    };
}  // namespace detail


template <class T>
class sequence_wrapper
{
private:
    //py::detail::sequence_iterator;
    using raw_iterator = decltype(std::declval<py::sequence>().begin());

public:
    using const_iterator = detail::sequence_iterator_wrapper<T, raw_iterator>;
    using value_type = T;
    using size_type = std::size_t;

    sequence_wrapper(const py::sequence& seq): m_seq(seq) {};

    bool empty() const { return m_seq.size() == 0; }

    const_iterator begin() const { return std::cbegin(m_seq); }
    const_iterator end() const { return std::cend(m_seq); }

    size_t size() const { return m_seq.size(); }

    decltype(auto)
    operator[] (const size_t i) const { return py::cast<T>(m_seq[i]); }

private:
    const py::sequence& m_seq;
};
