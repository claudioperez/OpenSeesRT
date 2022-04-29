#include <pybind11/pybind11.h>
#include <HystereticBackbone.h>
#include <ManderBackbone.h>

namespace py = pybind11;

class PyHystereticBackbone : public HystereticBackbone {
public:
    /* Inherit the constructors */
    using HystereticBackbone::HystereticBackbone;

    /* Trampoline (need one for each virtual function) */
    double getStress(double strain) override {
        PYBIND11_OVERRIDE_PURE(
            double,                 /* Return type */
            HystereticBackbone,     /* Parent class */
            getStress,    /* Name of function in C++ (must match Python name) */
            strain       /* Argument(s) */
        );
    }
};


PYBIND11_MODULE(libOpenSeesRT, m) {
    py::class_<HystereticBackbone, PyHystereticBackbone>(m, "HystereticBackbone");
       // .def(py::init<>())
       // .def("getStress", &HystereticBackbone::getStress);

    py::class_<ManderBackbone, HystereticBackbone>(m, "PopovicsBackbone")
      .def(py::init<int, double, double, double>(),
           py::arg("tag"), py::arg("f"), py::arg("e"), py::arg("E")
      )
      .def("getStress", &ManderBackbone::getStress)
      ;

}

