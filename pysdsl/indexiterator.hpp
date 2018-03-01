#pragma once

#include <cstddef>
#include <iterator>
#include <utility>

#include <sdsl/wavelet_trees.hpp>


namespace detail
{

    template <class Container,
              typename T = typename Container::value_type,
              typename S = typename Container::size_type>
    class count_index_iterator
    {
    public:
        typedef S                                     difference_type;
        typedef S                                           size_type;
        typedef T                                          value_type;
        typedef T*                                            pointer;
        typedef T&                                          reference;
        typedef const T&                              const_reference;
        typedef std::random_access_iterator_tag     iterator_category;

        count_index_iterator(): m_parent(nullptr), m_index(0) {}

        count_index_iterator(const Container* container, S index):
            m_parent(container),
            m_index(index)
        {}

        count_index_iterator(const count_index_iterator& other) = default;
        count_index_iterator(count_index_iterator&& other) = default;

        count_index_iterator&
        operator=(const count_index_iterator& other) = default;

        count_index_iterator&
        operator=(count_index_iterator&& other) = default;

        bool operator!=(const count_index_iterator& other) const
        { return (m_parent != other.m_parent) || (m_index != other.m_index); }

        bool operator==(const count_index_iterator& other) const
        { return (m_parent == other.m_parent) && (m_index == other.m_index); }

        bool operator>(const count_index_iterator& other) const
        { return m_index > other.m_index; }

        bool operator>=(const count_index_iterator& other) const
        { return m_index >= other.m_index; }

        bool operator<(const count_index_iterator& other) const
        { return m_index < other.m_index; }

        bool operator<=(const count_index_iterator& other) const
        { return m_index <= other.m_index; }

        value_type operator*() { return (*m_parent)[m_index]; }

        decltype(auto) operator++() {
            ++m_index;
            return *this;
        }

        decltype(auto) operator++(int)
        { return *count_index_iterator<Container, T, S>(m_parent, m_index++); }

        decltype(auto) operator--() {
            --m_index;
            return *this;
        }

        decltype(auto) operator--(int)
        { return *count_index_iterator<Container, T, S>(m_parent, m_index--); }

        difference_type operator-(const count_index_iterator& other) const
        { return m_index - other.m_index; }

        decltype(auto) operator-(const difference_type step) const
        {
            return count_index_iterator<Container, T, S>(m_parent,
                                                         m_index - step);
        }

        decltype(auto) operator+(const difference_type step) const
        {
            return count_index_iterator<Container, T, S>(m_parent,
                                                         m_index + step);
        }

        friend decltype(auto) operator+(const difference_type step,
                                        const count_index_iterator &self)
        {
            return count_index_iterator<Container, T, S>(self.m_parent,
                                                         self.m_index + step);
        }

        friend decltype(auto) operator-(const difference_type step,
                                        const count_index_iterator &self)
        {
            return count_index_iterator<Container, T, S>(self.m_parent,
                                                         self.m_index - step);
        }

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

        const_reference operator[](difference_type i) const
        { return (*m_parent)[m_index + i]; }

        void swap(count_index_iterator& other)
        {
            std::swap(m_parent, other.m_parent);
            std::swap(m_index, other.m_index);
        }

    private:
        const Container* m_parent;
        S m_index;
    };


    template <class Container,
              typename T = typename Container::value_type,
              typename S = typename Container::size_type>
    void swap(count_index_iterator<Container, T, S>& first,
              count_index_iterator<Container, T, S>& second)
    {
        first.swap(second);
    }

}


template <class... T>
auto cbegin(const sdsl::wt_gmr<T...>& container) {
    return detail::count_index_iterator<sdsl::wt_gmr<T...>>(&container, 0);
}

template <class... T>
auto cend(const sdsl::wt_gmr<T...>& container) {
    return detail::count_index_iterator<sdsl::wt_gmr<T...>>(&container,
                                                            container.size());
}

template <class... T>
auto cbegin(const sdsl::wt_gmr_rs<T...>& container) {
    return detail::count_index_iterator<sdsl::wt_gmr_rs<T...>>(&container, 0);
}

template <class... T>
auto cend(const sdsl::wt_gmr_rs<T...>& container) {
    return detail::count_index_iterator<sdsl::wt_gmr_rs<T...>>(
        &container, container.size()
    );
}

template <class...T>
auto cbegin(const sdsl::wt_ap<T...>& container) {
    return detail::count_index_iterator<sdsl::wt_ap<T...>>(&container, 0);
}

template <class...T>
auto cend(const sdsl::wt_ap<T...>& container) {
    return detail::count_index_iterator<sdsl::wt_ap<T...>>(&container,
                                                           container.size());
}
