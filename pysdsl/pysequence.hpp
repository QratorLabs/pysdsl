#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;


namespace detail
{

    template <class T, class Base>
    class sequence_iterator_wrapper
    {
    public:
        sequence_iterator_wrapper(Base it): m_it(it) {}

        bool operator!=(const sequence_iterator_wrapper& other)
        {
            return m_it != other.m_it;
        }

        decltype(auto) operator*() { return py::cast<T>(*m_it); }

        decltype(auto) operator++() {
            ++m_it;
            return *this;
        }

        decltype(auto) operator++(int) {
            return *sequence_iterator_wrapper<T, Base>(m_it++);
        }

        decltype(auto) operator=(sequence_iterator_wrapper other)
        {
            if (this != &other) {
                m_it = other.m_it;
            }
            return *this;
        }

    private:
        Base m_it;
    };

}


template <class T>
class sequence_wrapper
{
private:
    // dummy function to help with decltype(m_seq)::const_iterator
    decltype(auto) _begin() const { return m_seq.begin(); }
    typedef decltype(_begin) raw_;
    typedef std::result_of_t<raw_> raw_iterator;

public:
    typedef detail::sequence_iterator_wrapper<T, raw_iterator> const_iterator;
    typedef T value_type;

    sequence_wrapper(const py::sequence& seq): m_seq(seq) {};

    bool empty() const { return m_seq.size() == 0; }

    const_iterator begin() const { return m_seq.begin(); }
    const_iterator end() const { return m_seq.end(); }

    size_t size() const { return m_seq.size(); }

    decltype(auto)
    operator[] (const size_t i) const { return py::cast<T>(m_seq[i]); }

private:
    const py::sequence& m_seq;
};
