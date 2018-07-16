#pragma once

#include <cstddef>
#include <iterator>
#include <utility>

namespace {
    template<typename> struct int_ { using type = int; }; }


namespace detail
{
    template <class Container, typename T, typename TRef = T&>
    struct get_reference { using type = TRef; };

    template <class Container, typename T>
    struct get_reference<Container, T, typename Container::reference>
    { using type = typename Container::reference; };

    template <class Container,
              typename T = typename Container::value_type,
              typename S = typename Container::size_type,
              typename TRef = typename get_reference<Container, T>::type>
    class count_index_iterator
    {
    public:
        using difference_type = S;
        using size_type = S;
        using value_type = T;
        using pointer = T*;
        using reference = TRef;
        using const_reference = const T&;
        using iterator_category = std::random_access_iterator_tag;

        constexpr count_index_iterator() noexcept:
            m_parent(nullptr), m_index(0)
        {}

        constexpr count_index_iterator(const Container* container,
                                       S index) noexcept:
            m_parent(container),
            m_index(index)
        {}
        ~count_index_iterator() = default;

        constexpr count_index_iterator(
            const count_index_iterator& other
        ) noexcept = default;

        constexpr
        count_index_iterator(count_index_iterator&& other) noexcept = default;

        constexpr count_index_iterator&
        operator=(const count_index_iterator& other) noexcept = default;

        constexpr count_index_iterator&
        operator=(count_index_iterator&& other) noexcept = default;

        constexpr
        bool operator!=(const count_index_iterator& other) const noexcept {
            return (m_parent != other.m_parent) || (m_index != other.m_index); }

        constexpr
        bool operator==(const count_index_iterator& other) const noexcept {
            return (m_parent == other.m_parent) && (m_index == other.m_index); }

        constexpr
        bool operator>(const count_index_iterator& other) const noexcept {
            return m_index > other.m_index; }

        constexpr
        bool operator>=(const count_index_iterator& other) const noexcept {
            return m_index >= other.m_index; }

        constexpr
        bool operator<(const count_index_iterator& other) const noexcept {
            return m_index < other.m_index; }

        constexpr
        bool operator<=(const count_index_iterator& other) const noexcept {
            return m_index <= other.m_index; }

        value_type operator*() { return (*m_parent)[m_index]; }

        decltype(auto) operator++()
        {
            ++m_index;
            return *this;
        }

        const auto operator++(int) {
            return *count_index_iterator<Container, T, S>(m_parent,
                                                          m_index++); }

        decltype(auto) operator--()
        {
            --m_index;
            return *this;
        }

        const auto operator--(int) {
            return *count_index_iterator<Container, T, S>(m_parent,
                                                          m_index--); }

        constexpr
        difference_type
        operator-(const count_index_iterator& other) const noexcept {
            return m_index - other.m_index; }

        constexpr
        decltype(auto) operator-(const difference_type step) const noexcept {
            return count_index_iterator<Container, T, S>(m_parent,
                                                         m_index - step); }

        constexpr
        decltype(auto) operator+(const difference_type step) const noexcept {
            return count_index_iterator<Container, T, S>(m_parent,
                                                         m_index + step); }

        friend constexpr decltype(auto)
        operator+(const difference_type step,
                  const count_index_iterator &self) noexcept {
            return count_index_iterator<Container, T, S>(self.m_parent,
                                                         self.m_index + step); }

        friend constexpr decltype(auto)
        operator-(const difference_type step,
                  const count_index_iterator &self) noexcept {
            return count_index_iterator<Container, T, S>(self.m_parent,
                                                         self.m_index - step); }

        count_index_iterator& operator+=(const difference_type i)
        {
            m_index += i;
            return *this;
        }

        count_index_iterator& operator-=(difference_type i)
        {
            m_index -= i;
            return *this;
        }

        const_reference operator[](difference_type i) const {
            return (*m_parent)[m_index + i]; }

        void swap(count_index_iterator& other) noexcept
        {
            std::swap(m_parent, other.m_parent);
            std::swap(m_index, other.m_index);
        }

    private:
        const Container* m_parent;
        S m_index;
    };


    template <typename... P>
    void swap(count_index_iterator<P...>& first,
              count_index_iterator<P...>& second)
    {
        first.swap(second);
    }
}  // namespace detail
