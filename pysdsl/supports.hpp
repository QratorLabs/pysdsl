#pragma once

#include <string>

#include <sdsl/bit_vectors.hpp>

#include <pybind11/pybind11.h>

#include "docstrings.hpp"
#include "io.hpp"


namespace py = pybind11;


template <class T>
class support_helper
{
private:
    const sdsl::bit_vector& m_vec;
    const T m_support;
public:
    typedef T type;

    support_helper(const sdsl::bit_vector& vec, const T&& support):
        m_vec(vec),
        m_support(std::move(support))
    {}
    auto size() const { return m_vec.size(); }
    auto operator()(size_t idx) const { return m_support(idx); }

    operator const T&() const { return m_support; }
};


template <class Base>
inline
auto add_support_class(py::module &m,
                       const std::string&& name,
                       const std::string&& method_name,
                       const std::string&& doc_call,
                       const char* doc = nullptr)
{
    auto cls = py::class_<Base>(m, name.c_str());

    cls.def(
        method_name.c_str(),
        [](const Base& self, size_t idx)
        {
            if (idx >= self.size())
            {
                throw std::out_of_range(std::to_string(idx));
            }
            return self(idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        py::arg("idx"),
        doc_call.c_str()
    );
    cls.attr("__call__") = cls.attr(method_name.c_str());

    add_description(cls);

    if (doc) cls.doc() = doc;

     return cls;
}


template <class S, class T>
inline
decltype(auto) bind_support(const S*,
                            py::class_<T>& cls, const std::string& call_name,
                            const char* alt_name=nullptr)
{
    cls.def(
        call_name.c_str(),
        [](T& self)
        {
            S support;
            sdsl::util::init_support(support, &self);

            return support;
        },
        py::keep_alive<0, 1>()
    );

    if (alt_name) cls.attr(alt_name) = cls.attr(call_name.c_str());

    return cls;
}


template <class S, class T>
inline
decltype(auto) bind_support(const support_helper<S>*,
                            py::class_<T>& cls, const std::string& call_name,
                            const char* alt_name=nullptr)
{
    cls.def(
        call_name.c_str(),
        [](T& self)
        {
            S support;
            sdsl::util::init_support(support, &self);

            return support_helper<S>(self, std::move(support));
        },
        py::keep_alive<0, 1>()
    );

    if (alt_name) cls.attr(alt_name) = cls.attr(call_name.c_str());

    return cls;
}


template <class T,
          class R0=typename T::rank_0_type, class R1=typename T::rank_1_type>
inline
auto add_rank_support(py::module &m, py::class_<T>& cls,
                      const std::string& base_name,
                      const char* suffix = "",
                      bool defaults = true,
                      const std::string s0 = "0", const std::string s1 = "1",
                      const char* doc_rank = nullptr)
{
    add_support_class<R0>(m, base_name + "Rank" + suffix + "_" + s0, "rank",
                          "number of patterns `" + s0 + "` in the prefix "
                          "[0..idx) in supported vector", doc_rank);
    bind_support((R0 *)nullptr, cls,
                 std::string("init_rank") + suffix + "_" + s0);

    add_support_class<R1>(m, base_name + "Rank" + suffix + "_" + s1, "rank",
                          "number of patterns `" + s1 + "` in the prefix "
                          "[0..idx) in supported vector", doc_rank);
    bind_support((R1 *)nullptr, cls,
                 std::string("init_rank") + suffix + "_" + s1,
                 defaults ?
                    (std::string("init_rank") + suffix).c_str() :
                    nullptr);

    return cls;
}


template <class T,
          class S0=typename T::select_0_type,
          class S1=typename T::select_1_type>
inline
auto add_select_support(py::module &m, py::class_<T>& cls,
                        const std::string& base_name,
                        const char* suffix = "",
                        bool defaults = true,
                        const std::string s0 = "0", const std::string s1 = "1",
                        const char* doc_select = nullptr)
{
    add_support_class<S0>(m, base_name + "Select" + suffix + "_" + s0, "select",
                          "position of the idx-th pattern `" + s0 +
                          "` in supported vector", doc_select);
    bind_support((S0 *)nullptr,
                 cls, std::string("init_select") + suffix + "_" + s0);

    add_support_class<S1>(m, base_name + "Select" + suffix + "_" + s1, "select",
                          "position of the idx-th pattern `" + s1 +
                          "` in supported vector", doc_select);
    bind_support((S1 *)nullptr, cls,
                 std::string("init_select") + suffix + "_" + s1,
                 defaults ?
                    (std::string("init_select") + suffix).c_str() :
                    nullptr);

    return cls;
}


inline
void add_bitvector_supports(py::module& m, py::class_<sdsl::bit_vector>& cls)
{
    add_rank_support<sdsl::bit_vector,
                    sdsl::rank_support_v<0, 1>,
                    sdsl::rank_support_v<1, 1>>(
        m, cls, "_BitVector", "V", true, "0", "1", doc_rank_v
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_v<00, 2>,
                     sdsl::rank_support_v<01, 2>>(
        m, cls, "_BitVector", "V", false, "00", "01", doc_rank_v
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_v<10, 2>,
                     sdsl::rank_support_v<11, 2>>(
        m, cls, "_BitVector", "V", false, "10", "11", doc_rank_v
    );
    add_rank_support<sdsl::bit_vector,
                    sdsl::rank_support_v5<0, 1>,
                    sdsl::rank_support_v5<1, 1>>(
        m, cls, "_BitVector", "V5", false, "0", "1", doc_rank_v5
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_v5<00, 2>,
                     sdsl::rank_support_v5<01, 2>>(
        m, cls, "_BitVector", "V5", false, "00", "01", doc_rank_v5
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_v5<10, 2>,
                     sdsl::rank_support_v5<11, 2>>(
        m, cls, "_BitVector", "V5", false, "10", "11", doc_rank_v5
    );
    cls.attr("init_rank") = cls.attr("init_rankV");
    cls.attr("init_rank_0") = cls.attr("init_rankV_0");
    cls.attr("init_rank_1") = cls.attr("init_rankV_1");

    add_select_support<sdsl::bit_vector,
                       support_helper<sdsl::select_support_mcl<0, 1>>,
                       support_helper<sdsl::select_support_mcl<1, 1>>>(
        m, cls, "_BitVector", "MCL", true, "0", "1", doc_select_mcl
    );
    add_select_support<sdsl::bit_vector,
                       support_helper<sdsl::select_support_mcl<10, 2>>,
                       support_helper<sdsl::select_support_mcl<11, 2>>>(
        m, cls, "_BitVector", "MCL", false, "10", "11", doc_select_mcl
    );
    cls.attr("init_select") = cls.attr("init_selectMCL");
    cls.attr("init_select_0") = cls.attr("init_selectMCL_0");
    cls.attr("init_select_1") = cls.attr("init_selectMCL_1");

    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_scan<0, 1>,
                     sdsl::rank_support_scan<1, 1>>(
        m, cls, "_BitVector", "Scan", false, "0", "1", doc_rank_scan
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_scan<00, 2>,
                     sdsl::rank_support_scan<01, 2>>(
        m, cls, "_BitVector", "Scan", false, "00", "01", doc_rank_scan
    );
    add_rank_support<sdsl::bit_vector,
                     sdsl::rank_support_scan<10, 2>,
                     sdsl::rank_support_scan<11, 2>>(
        m, cls, "_BitVector", "Scan", false, "10", "11", doc_rank_scan
    );

    add_select_support<sdsl::bit_vector,
                       support_helper<sdsl::select_support_scan<0, 1>>,
                       support_helper<sdsl::select_support_scan<1, 1>>>(
        m, cls, "_BitVector", "Scan", false, "0", "1", doc_select_scan
    );
    add_select_support<sdsl::bit_vector,
                       support_helper<sdsl::select_support_scan<10, 2>>,
                       support_helper<sdsl::select_support_scan<01, 2>>>(
        m, cls, "_BitVector", "Scan", false, "10", "01", doc_select_scan
    );
}
