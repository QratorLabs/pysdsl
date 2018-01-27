/*cppimport
<%
cfg['compiler_args'] = ['-std=c++14']
cfg['include_dirs'] = ['sdsl-lite/include']
%>
*/

#include "sdsl/vectors.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;


PYBIND11_MODULE(sdsl, m) {
}
