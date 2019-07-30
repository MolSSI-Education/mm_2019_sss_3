#include "potentials.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>


PYBIND11_MODULE(potentials_cpp, m)
{
    m.doc() = "Pairwiswe potential energy at considered distance, implemented in C++";
    
    m.def("LJ", LJ, "Pairwise potential and correction energy by Lennard-Jones potential.");
    m.def("cutoff_correction", cutoff_correction, "The function corrects interaction energy from energy cutoff.");
}