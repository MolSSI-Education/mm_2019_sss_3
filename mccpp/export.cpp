#include "pairwise.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>


PYBIND11_MODULE(pairwise, m)
{
    pybind11::class_<PairwisePotentials> pairwise(m, "PairwisePotentials");
    pairwise
    .def("potential", &PairwisePotentials::potential);

    pybind11::class_<LJ>(m, "LJ", pairwise)
    .def(pybind11::init<double, double, double>(),
        pybind11::arg("sigma")=1, pybind11::arg("epsilon")=1,
        pybind11::arg("cutoff")=2.6)
    .def("potential", &LJ::potential)
    .def("__call__", &LJ::potential)
    .def("cutoff_correction", &LJ::cutoff_correction);

/* 
    pybind11::class_<Box>(m, "Box")
    .def(pybind11::init<>())
    .def_readwrite("box_dims", &Box::box_dims, "box dimensions")
    .def_readwrite("coordinates", &Box::coordinates, "particle coordinates")
    .def("volume", &Box::volume)
    .def("coord_wrap", &Box::coord_wrap)
    .def("minimum_image_dist", &Box::minimum_image_dist);
    */
}
