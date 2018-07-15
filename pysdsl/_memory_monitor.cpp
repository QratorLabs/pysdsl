#include <string>
#include <sstream>

#include <pybind11/pybind11.h>

#include <sdsl/memory_management.hpp>

namespace py = pybind11;


PYBIND11_MODULE(_memory_monitor, m)
{
    m.doc() = "Internals of memory monitor module";

    m.def("start", [] () { return sdsl::memory_monitor::start(); });
    m.def("stop", [] () { return sdsl::memory_monitor::stop(); });
    m.def(
        "report",
        [] ()
        {
            std::stringstream fout;
            sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(fout);
            auto json = py::module::import("json");
            return json.attr("loads")(fout.str());
        }
    );
    m.def(
        "report_json",
        [] ()
        {
            std::stringstream fout;
            sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(fout);
            return fout.str();
        }
    );
    m.def(
        "report_html",
        [] ()
        {
            std::stringstream fout;
            sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(fout);
            return fout.str();
        }
    );
    m.def(
        "report_html",
        [](const std::string& file_name)
        {
            std::ofstream fout;
            fout.open(file_name, std::ios::out | std::ios::binary);
            if (!fout.good()) throw std::runtime_error("Can't write to file");
            sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(fout);
            if (!fout.good()) throw std::runtime_error("Error during write");
            fout.close();
        },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>()
    );
    m.def(
        "report_json",
        [](const std::string& file_name)
        {
            std::ofstream fout;
            fout.open(file_name, std::ios::out | std::ios::binary);
            if (!fout.good()) throw std::runtime_error("Can't write to file");
            sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(fout);
            if (!fout.good()) throw std::runtime_error("Error during write");
            fout.close();
        },
        py::arg("file_name"),
        py::call_guard<py::gil_scoped_release>()
    );

}
