#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <elementAPI.h> // G3_getRuntime/SafeBuilder
#include "TclSafeBuilder.h"

#include <Vector.h>
#include <Node.h>
#include <Element.h>
#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>
#include <HystereticBackbone.h>

#include <ManderBackbone.h>

std::unique_ptr<TclSafeBuilder, py::nodelete> 
get_builder(py::object interpaddr) {
      void *interp_addr;
      interp_addr = (void*)PyLong_AsVoidPtr(interpaddr.ptr());
      void *builder_addr = G3_getSafeBuilder(G3_getRuntime((Tcl_Interp*)interp_addr));
      return std::unique_ptr<TclSafeBuilder, py::nodelete>((TclSafeBuilder*)builder_addr);
} // , py::return_value_policy::reference


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

py::array_t<double>
copy_vector(Vector vector)
{
  py::array_t<double> array(vector.Size());
  double *ptr = static_cast<double*>(array.request().ptr);
  for (int i=0; i<vector.Size(); i++)
    ptr[i] = vector(i);
  return array;
}


void
init_obj_module(py::module &m)
{
    py::class_<Vector, std::unique_ptr<Vector, py::nodelete>>(m, "Vector", py::buffer_protocol())
        .def (py::init([](
             py::array_t<double, py::array::c_style|py::array::forcecast> array
        )->Vector{
             bool verbose = true;
             py::buffer_info info = array.request();
             if (verbose){
                 py::print("ptr\t",info.ptr);
                 py::print("itemsize\t", info.itemsize);
                 py::print("format\t", info.format);
                 py::print("ndim\t", info.ndim);
                 py::print("shape\t", py::cast(info.shape));
                 py::print("strides\t", py::cast(info.strides));
                 py::print("array\t", array);
                 printf("%lf\n", *((double*)info.ptr));
             }
             return Vector(static_cast<double*>(info.ptr),(int)info.shape[0]);

        }))
        /* Allow reference by numpy array; requires access to Vector.theData */
        // .def_buffer([](Vector& v) -> py::buffer_info{
        //       return py::buffer_info(
        //       );
        // })

        /* pyg3.Vector(array:Seq, assert_size:int) */
        .def (py::init([](
             py::array_t<double, py::array::c_style|py::array::forcecast> array,
             int assert_size
        ) -> Vector {
             bool verbose = true;
             py::buffer_info info = array.request();
             if (verbose){
                 py::print("ptr\t",info.ptr);
                 py::print("itemsize\t", info.itemsize);
                 py::print("format\t", info.format);
                 py::print("ndim\t", info.ndim);
                 py::print("shape\t", py::cast(info.shape));
                 py::print("strides\t", py::cast(info.strides));
             }
             if (info.shape[0] != assert_size)
                 throw std::runtime_error("Incompatible buffer dimension.");
             return Vector(
                 static_cast<double*>(info.ptr), 
                 // static_cast<int>(info.itemsize)
                 static_cast<int>(info.shape[0])
             );
        }))
    ;
    py::class_<Matrix, std::unique_ptr<Matrix, py::nodelete>>(m, "Matrix", py::buffer_protocol())
        .def (py::init([](
             py::array_t<double, py::array::c_style|py::array::forcecast> array,
             int assert_size
        ) -> Matrix {
             bool verbose = true;
             py::buffer_info info = array.request();
             if (verbose){
                 py::print("ptr\t",info.ptr);
                 py::print("itemsize\t", info.itemsize);
                 py::print("format\t", info.format);
                 py::print("ndim\t", info.ndim);
                 py::print("shape\t", py::cast(info.shape));
                 py::print("strides\t", py::cast(info.strides));
             }
             if (info.shape[0] != assert_size)
                 throw std::runtime_error("Incompatible buffer dimension.");
             return Matrix(
                 static_cast<double*>(info.ptr), 
                 static_cast<int>(info.shape[0]),
                 static_cast<int>(info.shape[1])
             );
        }))
    ;
    
    py::class_<Node,    std::unique_ptr<Node,py::nodelete>>(m, "_Node")
    ;
    py::class_<Element, std::unique_ptr<Element,py::nodelete>>(m, "_Element")
      .def ("commitState",           &Element::commitState)
      .def ("revertToStart",         &Element::revertToStart)
      .def ("revertToLastCommit",    &Element::revertToLastCommit)
    ;
    py::class_<SectionForceDeformation, std::unique_ptr<SectionForceDeformation, py::nodelete> > (m, "_SectionForceDeformation")
      //.def ("setTrialSectionDeformation", &SectionForceDeformation::setTrialSectionDeformation)
      .def ("getSectionDeformation", &SectionForceDeformation::getSectionDeformation)
      .def ("getStressResultant",    &SectionForceDeformation::getStressResultant)
      .def ("commitState",           &SectionForceDeformation::commitState)
      .def ("revertToStart",         &SectionForceDeformation::revertToStart)
      .def ("revertToLastCommit",    &SectionForceDeformation::revertToLastCommit)
    ;
    py::class_<UniaxialMaterial, std::unique_ptr<UniaxialMaterial, py::nodelete>>(m, "_UniaxialMaterial")
      .def ("getStress",             &UniaxialMaterial::getStress)
      .def ("getTangent",            &UniaxialMaterial::getTangent)
      .def ("getDampTangent",        &UniaxialMaterial::getDampTangent)
      .def ("getStrainRate",         &UniaxialMaterial::getStrainRate)
      .def ("commitState",           &UniaxialMaterial::commitState)
      .def ("revertToStart",         &UniaxialMaterial::revertToStart)
      .def ("revertToLastCommit",    &UniaxialMaterial::revertToLastCommit)
    ;
    py::class_<HystereticBackbone, PyHystereticBackbone>(m, "HystereticBackbone")
      .def("getStress", &HystereticBackbone::getStress);
    ;

    py::class_<ManderBackbone, HystereticBackbone>(m, "PopovicsBackbone")
      .def(py::init<int, double, double, double>(),
           py::arg("tag"), py::arg("f"), py::arg("e"), py::arg("E")
      )
      .def("getStress", &ManderBackbone::getStress)
    ;




    py::class_<TclSafeBuilder, std::unique_ptr<TclSafeBuilder, py::nodelete> >(m, "TclTclSafeBuilder")  // 
                                                                                                        //
      .def (py::init([](py::object interpaddr)->std::unique_ptr<TclSafeBuilder, py::nodelete>{
            void *interp_addr;
            interp_addr = (void*)PyLong_AsVoidPtr(interpaddr.ptr());
            void *builder_addr = Tcl_GetAssocData((Tcl_Interp*)interp_addr, "OPS::theTclSafeBuilder", NULL);
            return std::unique_ptr<TclSafeBuilder, py::nodelete>((TclSafeBuilder*)builder_addr);
        }) // , py::return_value_policy::reference
      )
      .def ("getSection", [](TclSafeBuilder& builder, py::str id){
          return builder.getSection(id);
      })
      .def ("getUniaxialMaterial", [](TclSafeBuilder& builder, py::str id){
          return builder.getUniaxialMaterial(id);
      })
      //.def ("getUniaxialMaterial", &TclSafeBuilder::getUniaxialMaterial)
        //.def ("getDomain", &TclSafeBuilder::getDomainPtr);
    ;

    m.def ("get_builder", &get_builder
    );

    /*
    m.def ("builder_addr", [](py::object rtaddr)->std::unique_ptr<TclSafeBuilder, py::nodelete>{
        void *rt_addr;
        rt_addr = (void*)PyLong_AsVoidPtr(rtaddr.ptr());
        void *builder_addr = G3_getSafeBuilder((G3_Runtime*)rt_addr);
        return std::unique_ptr<TclSafeBuilder, py::nodelete>((TclSafeBuilder*)builder_addr);
      }
    );
    */
}
/*
void
init_builder_module(py::module &m)
{
}
*/
PYBIND11_MODULE(libOpenSeesRT, m) {
  init_obj_module(m);
}

