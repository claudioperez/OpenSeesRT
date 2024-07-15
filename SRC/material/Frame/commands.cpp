//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#include <tcl.h>
#include <ArgumentTracker.h>
#include <G3_Logging.h>
#include <string.h>
#include <runtime/BasicModelBuilder.h>
#include <ElasticLinearFrameSection3d.h>

int
validate(FrameSectionConstants& section, int ndm, int shear)
{
  int status = 0;
//if (section.E <= 0.0)  {
//  status -= 1;
//  //opserr << "Input E <= 0.0\n";
//}
  
  if (section.A <= 0.0) {
    status -= 1;
    //opserr << "Input A <= 0.0\n";
  }
  
  if (section.Iz <= 0.0) {
    status -= 1;
    //opserr << "Input Iz <= 0.0\n";
  }
  
  
//if (shear)
//  if (section.G <= 0.0) {
//    status -= 1;
//    //opserr << "Input G <= 0.0\n";
//  }
  
  if (ndm == 3) {
//  if (section.J <= 0.0) {
//    //opserr << "Input J <= 0.0\n";
//    status -= 1;
//  }
    if (section.Iy <= 0.0)  {
      //opserr << "Input Iy <= 0.0\n";
      status -= 1;
    }
  }
  return status;
}


// section ElasticFrame tag E? A? Iz? <Iy? G? J?>
//                          E  A  I
//                          E  A  I  G  J
//    -E     $E
//    -G     $G
//
//    -A     $A
//          {$A $Ay}      if ndm == 2
//          {$A $Ay $Az}
//    -Ay    $Ay
//    -Az    $Az
//
//    -I/B   $Iz
//          {$Iy $Iz}
//          {$Iy $Iz $Iyz}
//
//    -J     $J
//    -C     $Cw
//          {$Cw $Ca}
//
//    -Q    {$Qy $Qz <$Qyx $Qyz>}
//
//    -R    {$Qy $Qz}
//
int
TclCommand_newElasticFrame(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  FrameSectionConstants consts;

  int ndm = 0;
  int have_shear = 0;

  for (int i=0; i<argc; i++) {

    if (strcmp(argv[i], "-E") == 0) {
    }

    else if (strcmp(argv[i], "-A") == 0) {
      ++i;
      if (i >= argc) {
        return TCL_ERROR;
      }

      int listc;
      const char **list;
      if (Tcl_GetDouble(interp, argv[i], &consts.A) == TCL_OK) {
        continue;
      } else if (Tcl_SplitList(interp, argv[i], &listc, &list) == TCL_OK) {
      } else {
      }
    }
  }

  if (validate(consts, ndm, have_shear) != 0)
    return TCL_ERROR;
  
  return TCL_OK;

}

int
TclCommand_newElasticSection(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
    BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
    FrameSectionConstants consts{};

    enum class PositionalArgument : int {
      Tag, E, A, Iz, Iy, G, J, End
    };
    ArgumentTracker<PositionalArgument> tracker;

    FrameSection* theSection = nullptr;

    if (argc < 5) {
        opserr << OpenSees::PromptParseError << "insufficient arguments\n";
      //opserr << "Want: section Elastic tag? EA? EIz? <EIy? GJ?>" << "\n";
        opserr << "Want: section Elastic tag? E? A? Iz? <Iy? G? J?>" << "\n";
        return TCL_ERROR;
    }

    int tag;
    double E, G, J;

    bool use_mass = false;
    double mass=0.0;
    std::set<int> positional;

    int i;
    for (i=2; i<argc; i++) {

      if (strcmp(argv[i], "-mass") == 0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &mass) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid mass" << "\n";
          return TCL_ERROR;
        }
        use_mass = true;
      }

      else if ((strcmp(argv[i], "-youngs-modulus") == 0) ||
               (strcmp(argv[i], "-E") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &E) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Young's modulus." << "\n";
          return TCL_ERROR;
        }
        tracker.consume(PositionalArgument::E);
      }

      else if ((strcmp(argv[i], "-shear-modulus") == 0) ||
               (strcmp(argv[i], "-G") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &G) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid shear modulus." << "\n";
          return TCL_ERROR;
        }
        tracker.consume(PositionalArgument::G);
      }

      else if ((strcmp(argv[i], "-area") == 0) ||
               (strcmp(argv[i], "-A") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.A) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid area." << "\n";
          return TCL_ERROR;
        }
        tracker.consume(PositionalArgument::A);
      }

      else if ((strcmp(argv[i], "-inertia") == 0) ||
               (strcmp(argv[i], "-I") == 0) ||
               (strcmp(argv[i], "-Iz") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Iz) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid inertia." << "\n";
          return TCL_ERROR;
        }
        tracker.consume(PositionalArgument::Iz);
      }

      else if ((strcmp(argv[i], "-inertia-y") == 0) ||
               (strcmp(argv[i], "-Iy") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Iy) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid inertia." << "\n";
          return TCL_ERROR;
        }
        tracker.consume(PositionalArgument::Iy);
      }

      else
        positional.insert(i);

    }

    // positional arguments
    for (int i : positional) {
      switch (tracker.current()) {
        case PositionalArgument::Tag :
          if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid section Elastic tag" << "\n";
              return TCL_ERROR;           
          } else {
            tracker.increment();
            break;
          }

        case PositionalArgument::E:
          if (Tcl_GetDouble (interp, argv[i], &E) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid E" << "\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case PositionalArgument::A:
          if (Tcl_GetDouble (interp, argv[i], &consts.A) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid A" << "\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case PositionalArgument::Iz:
          if (Tcl_GetDouble (interp, argv[i], &consts.Iz) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid Iz" << "\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case PositionalArgument::Iy:
          if (Tcl_GetDouble (interp, argv[i], &consts.Iy) != TCL_OK || consts.Iy < 0) {
              opserr << OpenSees::PromptParseError << "invalid Iy" << "\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case PositionalArgument::G:
          if (Tcl_GetDouble (interp, argv[i], &G) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid G" << "\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case PositionalArgument::J:
          if (Tcl_GetDouble (interp, argv[i], &J) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid J" << "\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case PositionalArgument::End:
          opserr << OpenSees::PromptParseError << "unexpected argument" << argv[i] << "\n";
          return TCL_ERROR;
      }
    }

    if (tracker.current() != PositionalArgument::End) {
      opserr << OpenSees::PromptParseError
             << "missing required positional arguments\n";
      return TCL_ERROR;
    }

    if (argc > 8) {
      theSection = new ElasticLinearFrameSection3d(tag, E, consts.A, consts.Iz, consts.Iy, G, J, mass, use_mass);
    }
    else
      theSection = new ElasticLinearFrameSection3d(tag, E, consts.A, consts.Iz, mass, use_mass);

    if (theSection == nullptr || builder->addTaggedObject<FrameSection>(*theSection) < 0) {
      if (theSection != nullptr)
        delete theSection;
      return TCL_ERROR;
    } else
      return TCL_OK;

}
