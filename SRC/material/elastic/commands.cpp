/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
#include <tcl.h>
#include <string.h>
#include <BasicModelBuilder.h>

#include <ElasticMaterial.h>
#include <ElasticIsotropic.h>
#include <ElasticIsotropicMaterial.h>
#include <ElasticOrthotropicMaterial.h>
#include <material/Solid/ElasticIsotropicThreeDimensional.h>
#include <material/nD/ElasticIsotropicAxiSymm.h>
#include <material/nD/ElasticIsotropicMaterial.h>
#include <material/nD/ElasticIsotropicMaterialThermal.h>
#include <material/nD/ElasticIsotropicPlateFiber.h>

#include <material/nD/ElasticIsotropic3DThermal.h>
#include <material/Plane/ElasticIsotropicPlaneStress2D.h>
#include <material/Plane/ElasticIsotropicPlaneStrain2D.h>
#include <material/Frame/Fiber/ElasticIsotropicBeamFiber2d.h>
#include <material/Frame/Fiber/ElasticIsotropicBeamFiber.h>

#include <ElasticIsotropicMaterial.h>
// #include <ElasticCrossAnisotropic.h>
#include <PlaneStressMaterial.h>
#include <PlateFiberMaterial.h>
#include <BeamFiberMaterial.h>


enum class MaterialSymmetry {
  Triclinic     , // 21
  Monoclinic    , // 13
  Orthorhombic  , //  9
  Tetragonal    , //  7
//Tetragonal    , //  6
  Rhombohedral  , //  7
//Rhombohedral  , //  6
  Hexagonal     , //  5
  Cubic         , //  3
  Isotropic       //  2
};


int
TclCommand_newElasticMaterial(ClientData clientData, Tcl_Interp* interp, int argc, const char**const argv)
{
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  MaterialSymmetry symm = MaterialSymmetry::Isotropic;
  PlaneType type = PlaneType::None;

  if ((strcmp(argv[1], "ElasticIsotropic") == 0) || 
      (strcmp(argv[1], "Elastic") == 0) ||
      (strcmp(argv[1], "ElasticIsotropic3D") == 0))

  {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial ElasticIsotropic tag? E? v? <rho?>" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double E, v;
    double rho = 0.0;

    int loc = 2;
    if (Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
      opserr << "WARNING invalid ElasticIsotropic tag" << "\n";
      return TCL_ERROR;
    }
    loc++;

    if (Tcl_GetDouble(interp, argv[loc], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      return TCL_ERROR;
    }
    loc++;

    if (Tcl_GetDouble(interp, argv[loc], &v) != TCL_OK) {
      opserr << "WARNING invalid v\n";
      return TCL_ERROR;
    }
    loc++;

    if (argc > loc && Tcl_GetDouble(interp, argv[loc], &rho) != TCL_OK) {
      opserr << "WARNING invalid rho\n";
      return TCL_ERROR;
    }
    loc++;

    if (argc > loc && strcmp(argv[loc], "-plane-strain") ==0) {
      type = PlaneType::Strain;
    }
    loc++;

    switch (type) {
      case PlaneType::None:
        builder->addTaggedObject<UniaxialMaterial>(*new ElasticMaterial(tag, E, 0.0, E));
        builder->addTaggedObject<NDMaterial>(*new ElasticIsotropicMaterial(tag, E, v, rho));
        builder->addTaggedObject<Mate<3>>   (*new ElasticIsotropic<3>(tag, E, v, rho));
        break;
      case PlaneType::Strain:
        builder->addTaggedObject<NDMaterial>(*new ElasticIsotropicPlaneStrain2D(tag, E, v, rho));
        builder->addTaggedObject<Mate<2>>   (*new ElasticIsotropic<2,PlaneType::Strain>(tag, E, v, rho));
        break;
      case PlaneType::Stress:
        builder->addTaggedObject<NDMaterial>(*new ElasticIsotropicPlaneStress2D(tag, E, v, rho));
        builder->addTaggedObject<Mate<2>>   (*new ElasticIsotropic<2,PlaneType::Stress>(tag, E, v, rho));
        break;
    }
    return TCL_OK;
  }

  else if (strcmp(argv[1], "ElasticCrossAnisotropic") == 0)
  {
    if (argc < 8) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial ElasticCrossAnisotropic tag? Ehh? Ehv? nuhv? nuvv? Ghv? <rho?>"
             << "\n";
      return TCL_ERROR;
    }

    int tag;
    double Eh, Ev, nuhv, nuhh, Ghv;
    double rho = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid ElasticCrossAnisotropic tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &Eh) != TCL_OK) {
      opserr << "WARNING invalid Eh\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &Ev) != TCL_OK) {
      opserr << "WARNING invalid Ev\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &nuhv) != TCL_OK) {
      opserr << "WARNING invalid nuhv\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &nuhh) != TCL_OK) {
      opserr << "WARNING invalid nuhh\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &Ghv) != TCL_OK) {
      opserr << "WARNING invalid Ghv\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (argc > 8 && Tcl_GetDouble(interp, argv[8], &rho) != TCL_OK) {
      opserr << "WARNING invalid rho\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    // theMaterial = new ElasticCrossAnisotropic(tag, Eh, Ev, nuhv, nuhh, Ghv, rho);
  }
}

