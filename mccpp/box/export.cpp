#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "box.hpp"

PYBIND11_MODULE(sss_cpp, m)
{
  m.doc() = "This is example c++ module called from python";

  pybind11::class_<Box>(m, "Box")
    .def(pybind11::init<>())
    .def_readwrite("box_dims", &Box::box_dims, "box dimensions")
    .def_readwrite("coordinates", &Box::coordinates, "particle coordinates")
    .def("volume", &Box::volume)
    .def("coord_wrap", &Box::coord_wrap)
    .def("minimum_image_dist", &Box::minimum_image_dist);
}
