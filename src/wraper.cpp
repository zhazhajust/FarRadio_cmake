#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <cstring>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include "Tracer/Tracer.hpp"
#include "Field/Field.hpp"
#include "Detector/Detector.hpp"

#define CHECK_ERROR(err, str) \
    if (err < 0) { \
        std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << std::endl; \
        std::cerr << str << std::endl; \
        return; \
    } \

namespace py = pybind11;

PYBIND11_MODULE(faradio, m) {
    m.doc() = "Faradio python interface";
    py::class_<Field3D, std::shared_ptr<Field3D>>(m, "Field3D", py::buffer_protocol())
    .def(py::init<unsigned int, unsigned int, unsigned int>())
    .def("get_dim", &Field3D::get_dim)
    .def("__call__", &Field3D::operator())
    .def_buffer([](Field3D &m) -> py::buffer_info {
        return py::buffer_info(
            m.get_data(),                                                  /* Pointer to buffer */
            sizeof(double),                                                /* Size of one scalar */
            py::format_descriptor<double>::format(),                       /* Python struct-style format descriptor */
            3,                                                             /* Number of dimensions */
            { m.get_dim(0), m.get_dim(1), m.get_dim(2) },                  /* Buffer dimensions */
            { sizeof(double) * m.get_dim(1) * m.get_dim(2),                /* Strides (in bytes) for each index */
              sizeof(double) * m.get_dim(2),
              sizeof(double)}
        );
    })
    .def("to_memoryview", [](Field3D& m){
        int dim[3];
        dim[0] = m.get_dim(0);
        dim[1] = m.get_dim(1);
        dim[2] = m.get_dim(2);
        return py::memoryview::from_buffer(
            m.get_data(),               // buffer pointer
            {dim[0], dim[1], dim[2]},                 // shape (rows, cols)
            {sizeof(double) * dim[1] * dim[2], sizeof(double) * dim[2], sizeof(double)} // strides
        );
    });

    py::class_<SpheDetector>(m, "SpheDetector")
    .def(py::init<std::vector<double>, std::vector<double>, std::vector<int>>())
    .def("get_emf", &SpheDetector::get_emf)
    .def("get_screen_potisions", &SpheDetector::get_screen_potisions)
    .def("get_dmax", &SpheDetector::get_dmax)
    .def("get_dmin", &SpheDetector::get_dmin)
    .def("get_nf", &SpheDetector::get_nf)
#ifndef NONMPI
    .def("get_mpi", &SpheDetector::get_mpi)
#endif
    .def("set_approx", &SpheDetector::set_approx)
    .def("cmp_emf", &SpheDetector::cmp_emf)
    .def("reduce", &SpheDetector::reduce);
    
#ifndef NONMPI
    py::class_<FaradioMPI, std::shared_ptr<FaradioMPI>>(m, "FaradioMPI")
    .def(py::init<>())
    .def("getFaradioMPI", &FaradioMPI::getFaradioMPI)
    .def("rank", &FaradioMPI::rank)
    .def("size", &FaradioMPI::size);
#endif
    py::class_<Particle>(m, "Particle")
    .def(py::init<std::vector<std::vector<double>>, std::vector<std::vector<double>>,
        std::vector<std::vector<double>>, std::vector<std::vector<double>>,  
        double, double, double, int>())
    .def("get_position", &Particle::get_position)
    .def("get_beta", &Particle::get_beta)
    .def("get_time", &Particle::get_time);
};