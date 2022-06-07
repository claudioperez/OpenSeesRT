

#include <g3_api.h>


#include <SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.h>
void *OPS_ASDAbsorbingBoundary3D(void)
{
  static bool first_done = false;
  if (!first_done) {
    opserr << "Using ASDAbsorbingBoundary3D - Developed by: Massimo Petracca, "
              "Guido Camata, ASDEA Software Technology\n";
    first_done = true;
  }

  const char *descr =
      "Want: element ASDAbsorbingBoundary3D $tag $n1 $n2 $n3 $n4 $n5 $n6 $n7 "
      "$n8 $G $v $rho $btype <-fx $tsxTag> <-fy $tsyTag> <-fz $tszTag>\n";

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 13) {
    opserr << "ASDAbsorbingBoundary3D ERROR : Few arguments:\n" << descr;
    return 0;
  }

  // int parameters
  int iData[9];
  int numData = 9;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "ASDAbsorbingBoundary3D ERROR: Invalid integer mandatory values: "
              "element ASDAbsorbingBoundary3D wants 9 integer parameters\n"
           << descr;
    return 0;
  }

  // double parameters
  double dData[3];
  numData = 3;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "ASDAbsorbingBoundary3D ERROR: Invalid double mandatory values: "
              "element ASDAbsorbingBoundary3D wants 3 double parameters\n"
           << descr;
    return 0;
  }

  // string parameter
  const char *btype = OPS_GetString();
  int bflag = BND_NONE;
  if (strstr(btype, "B"))
    bflag |= BND_BOTTOM;
  if (strstr(btype, "L"))
    bflag |= BND_LEFT;
  if (strstr(btype, "R"))
    bflag |= BND_RIGHT;
  if (strstr(btype, "F"))
    bflag |= BND_FRONT;
  if (strstr(btype, "K"))
    bflag |= BND_BACK;
  if (bflag == BND_NONE) {
    opserr
        << "ASDAbsorbingBoundary3D ERROR: Invalid string mandatory value: the "
           "$btype "
           "argument should contain at least one of the following characters:\n"
           "'B', 'L', 'R', 'F', 'K'.\n"
        << descr;
    return 0;
  }

  // optional time series
  TimeSeries *fx = nullptr;
  TimeSeries *fy = nullptr;
  TimeSeries *fz = nullptr;
  // only on bottom boundaries
  if (bflag & BND_BOTTOM) {
    numData = 1;
    int tsTag = 0;
    // util: get x
    auto get_fx = [&numData, &tsTag, &fx, descr]() -> bool {
      if (OPS_GetInt(&numData, &tsTag) != 0) {
        opserr << "ASDAbsorbingBoundary3D ERROR: Invalid integer for -fx "
                  "optional time series.\n"
               << descr;
        return false;
      }
      fx = OPS_getTimeSeries(tsTag);
      if (fx == nullptr) {
        opserr << "ASDAbsorbingBoundary3D ERROR: Cannot find -fx time series "
                  "with id = "
               << tsTag << ".\n"
               << descr;
        return false;
      }
      return true;
    };
    // util: get y
    auto get_fy = [&numData, &tsTag, &fy, descr]() -> bool {
      if (OPS_GetInt(&numData, &tsTag) != 0) {
        opserr << "ASDAbsorbingBoundary3D ERROR: Invalid integer for -fy "
                  "optional time series.\n"
               << descr;
        return false;
      }
      fy = OPS_getTimeSeries(tsTag);
      if (fy == nullptr) {
        opserr << "ASDAbsorbingBoundary3D ERROR: Cannot find -fy time series "
                  "with id = "
               << tsTag << ".\n"
               << descr;
        return false;
      }
      return true;
    };
    // util: get z
    auto get_fz = [&numData, &tsTag, &fz, descr]() -> bool {
      if (OPS_GetInt(&numData, &tsTag) != 0) {
        opserr << "ASDAbsorbingBoundary3D ERROR: Invalid integer for -fz "
                  "optional time series.\n"
               << descr;
        return false;
      }
      fz = OPS_getTimeSeries(tsTag);
      if (fz == nullptr) {
        opserr << "ASDAbsorbingBoundary3D ERROR: Cannot find -fz time series "
                  "with id = "
               << tsTag << ".\n"
               << descr;
        return false;
      }
      return true;
    };
    // parse first
    numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs > 1) {
      const char *key = OPS_GetString();
      if (strcmp(key, "-fx") == 0) {
        if (!get_fx())
          return 0;
      } else if (strcmp(key, "-fy") == 0) {
        if (!get_fy())
          return 0;
      } else if (strcmp(key, "-fz") == 0) {
        if (!get_fz())
          return 0;
      } else {
        opserr << "ASDAbsorbingBoundary3D ERROR: Invalid optional flag \""
               << key << "\".\n"
               << descr;
        return 0;
      }
    }
    // parse second and third
    for (int i = 0; i < 2; ++i) {
      numArgs = OPS_GetNumRemainingInputArgs();
      if (numArgs > 1) {
        const char *key = OPS_GetString();
        if (strcmp(key, "-fx") == 0) {
          if (fx) {
            opserr
                << "ASDAbsorbingBoundary3D ERROR: -fx flag specified twice!.\n"
                << descr;
            return 0;
          }
          if (!get_fx())
            return 0;
        } else if (strcmp(key, "-fy") == 0) {
          if (fy) {
            opserr
                << "ASDAbsorbingBoundary3D ERROR: -fy flag specified twice!.\n"
                << descr;
            return 0;
          }
          if (!get_fy())
            return 0;
        } else if (strcmp(key, "-fz") == 0) {
          if (fz) {
            opserr
                << "ASDAbsorbingBoundary3D ERROR: -fz flag specified twice!.\n"
                << descr;
            return 0;
          }
          if (!get_fz())
            return 0;
        } else {
          opserr << "ASDAbsorbingBoundary3D ERROR: Invalid optional flag \""
                 << key << "\".\n"
                 << descr;
          return 0;
        }
      }
    }
  }

  // done
  return new ASDAbsorbingBoundary3D(
      iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], iData[6],
      iData[7], iData[8], dData[0], dData[1], dData[2], bflag, fx, fy, fz);
}
