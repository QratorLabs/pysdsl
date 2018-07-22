#pragma once

#include <fstream>

#include <sdsl/io.hpp>

#include <pybind11/pybind11.h>

#include "operations/iteration.hpp"


namespace py = pybind11;


template <class T> class support_helper;


template <class T, typename value_type = typename T::value_type>
decltype(auto) to_string(const T &self, const size_t max_elements=100,
                         const char* sep=", ", const char* start="[",
                                               const char* ends="]")
{
    std::ostringstream fout;
    fout.exceptions(std::ostringstream::failbit | std::ostringstream::badbit);

    fout << start;
    size_t count = 0;
    for (auto i = detail::cbegin(self); i != detail::cend(self); i++) {
        if (count) fout << sep;

        const value_type value = *i;

        fout << value;

        if (max_elements > 0 && count >= max_elements) {
            fout << sep << "...(" << self.size() << " elements)";
            break; }
        count++; }
    fout << ends;

    return fout.str();
}


template <class T, typename S = typename T::value_type>
inline auto add_to_string(py::class_<T>& cls)
{
    cls.def("__str__", [](const T& self) { return to_string<T, S>(self); });

    cls.def(
        "to_string", &to_string<T, S>,
        py::arg("max_elements") = 0, py::arg("sep") = ", ",
        py::arg("begin") = "[", py::arg("end") = "]");

    const auto name = py::cast<std::string>(cls.attr("__name__"));

    cls.def(
        "__repr__",
        [name] (const T& self) {
            return "<" + name +
                to_string<T, S>(self, 100, ", ", " [", "]>"); });

    return cls;
}


template <class T>
inline auto add_serialization(py::class_<T>& cls)
{
    cls.def(py::pickle(
        [](const T& self){
            std::stringstream fout;
            self.serialize(fout);
            return py::bytes(fout.str()); },
        [](const py::bytes& serialized){
            T result;
            std::stringstream fin(serialized);
            result.load(fin);
            return result; }));
    cls.def(
        "store_to_file",
        [](const T &self, const std::string& file_name) {
            return sdsl::store_to_file(self, file_name); },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>());

    cls.def_static(
        "load_from_file",
        [](const std::string& file_name) {
            T self;
            if (sdsl::load_from_file(self, file_name)) {
                return self; }
            throw std::exception(); },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>());

    cls.def(
        "store_to_checked_file",
        [](const T &self, const std::string& file_name) {
            return sdsl::store_to_checked_file(self, file_name); },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>() );

    cls.def_static(
        "load_from_checkded_file",
        [](const std::string& file_name) {
            T self;
            if (sdsl::load_from_checked_file(self, file_name)) {
                return self; }
            throw std::exception(); },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>());
    return cls;
}


template <class X, class T = typename X::type>
inline auto add_description(X& cls)
{
    typedef typename X::type P;
    cls.def(
        "write_structure_json",
        [](const P& self, const std::string& file_name) {
            std::ofstream fout;
            fout.open(file_name, std::ios::out | std::ios::binary);
            if (!fout.good()) throw std::runtime_error("Can't write to file");
            sdsl::write_structure<sdsl::JSON_FORMAT, T>(self, fout);
            if (!fout.good()) throw std::runtime_error("Error during write");
            fout.close(); },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>());
    cls.def(
        "write_structure_html",
        [](const P& self, const std::string& file_name) {
            std::ofstream fout;
            fout.open(file_name, std::ios::out | std::ios::binary);
            if (!fout.good()) throw std::runtime_error("Can't write to file");
            sdsl::write_structure<sdsl::HTML_FORMAT, T>(self, fout);
            if (!fout.good()) throw std::runtime_error("Error during write");
            fout.close(); },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>());

    cls.def_property_readonly(
        "structure_json",
        [](const P& self) {
            std::ostringstream fout;
            fout.exceptions(std::ostringstream::failbit |
                            std::ostringstream::badbit);

            sdsl::write_structure<sdsl::JSON_FORMAT, T>(self, fout);
            return fout.str(); },
        py::call_guard<py::gil_scoped_release>());

    cls.def_property_readonly(
        "structure_html",
        [](const P& self) {
            std::ostringstream fout;
            fout.exceptions(std::ostringstream::failbit |
                            std::ostringstream::badbit);

            sdsl::write_structure<sdsl::HTML_FORMAT, T>(self, fout);
            return fout.str();},
        py::call_guard<py::gil_scoped_release>());

    cls.def_property_readonly(
        "structure",
        [](const P& self) {
            std::ostringstream fout;
            fout.exceptions(std::ostringstream::failbit |
                            std::ostringstream::badbit);

            sdsl::write_structure<sdsl::JSON_FORMAT, T>(self, fout);
            auto json = py::module::import("json");
            return json.attr("loads")(fout.str()); });
    cls.def_property_readonly(
        "size_in_mega_bytes",
        [](const P &self) { return sdsl::size_in_mega_bytes<T>(self); });
    cls.def_property_readonly(
        "size_in_bytes",
        [](const P &self) { return sdsl::size_in_bytes<T>(self); });

    return cls;
}


template <class T>
inline auto add_description(py::class_<support_helper<T>>& cls)
{
    return add_description<py::class_<support_helper<T>>, T>(cls);
}
