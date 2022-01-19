

#include <g3_api.h>


#include <SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.h>
void *OPS_ASDAbsorbingBoundary2D(void)
{
  static bool first_done = false;
  if (!first_done) {
    opserr << "Using ASDAbsorbingBoundary2D - Developed by: Massimo Petracca, "
              "Guido Camata, ASDEA Software Technology\n";
    first_done = true;
  }

  const char *descr =
      "Want: element ASDAbsorbingBoundary2D $tag $n1 $n2 $n3 $n4 $G $v $rho "
      "$thickness $btype <-fx $tsxTag> <-fy $tsyTag>\n";

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 10) {
    opserr << "ASDAbsorbingBoundary2D ERROR : Few arguments:\n" << descr;
    return 0;
  }

  // int parameters
  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "ASDAbsorbingBoundary2D ERROR: Invalid integer mandatory values: "
              "element ASDAbsorbingBoundary2D wants 5 integer parameters\n"
           << descr;
    return 0;
  }

  // double parameters
  double dData[4];
  numData = 4;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "ASDAbsorbingBoundary2D ERROR: Invalid double mandatory values: "
              "element ASDAbsorbingBoundary2D wants 4 double parameters\n"
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
  if (bflag == BND_NONE) {
    opserr
        << "ASDAbsorbingBoundary2D ERROR: Invalid string mandatory value: the "
           "$btype "
           "argument should contain at least one of the following characters:\n"
           "'B', 'L', 'R'.\n"
        << descr;
    return 0;
  }

  // optional time series
  TimeSeries *fx = nullptr;
  TimeSeries *fy = nullptr;
  // only on bottom boundaries
  if (bflag & BND_BOTTOM) {
    numData = 1;
    int tsTag = 0;
    // util: get x
    auto get_fx = [&numData, &tsTag, &fx, descr]() -> bool {
      if (OPS_GetInt(&numData, &tsTag) != 0) {
        opserr << "ASDAbsorbingBoundary2D ERROR: Invalid integer for -fx "
                  "optional time series.\n"
               << descr;
        return false;
      }
      fx = OPS_getTimeSeries(tsTag);
      if (fx == nullptr) {
        opserr << "ASDAbsorbingBoundary2D ERROR: Cannot find -fx time series "
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
        opserr << "ASDAbsorbingBoundary2D ERROR: Invalid integer for -fy "
                  "optional time series.\n"
               << descr;
        return false;
      }
      fy = OPS_getTimeSeries(tsTag);
      if (fy == nullptr) {
        opserr << "ASDAbsorbingBoundary2D ERROR: Cannot find -fy time series "
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
      } else {
        opserr << "ASDAbsorbingBoundary2D ERROR: Invalid optional flag \""
               << key << "\".\n"
               << descr;
        return 0;
      }
    }
    // parse second
    numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs > 1) {
      const char *key = OPS_GetString();
      if (strcmp(key, "-fx") == 0) {
        if (fx) {
          opserr << "ASDAbsorbingBoundary2D ERROR: -fx flag specified twice!.\n"
                 << descr;
          return 0;
        }
        if (!get_fx())
          return 0;
      } else if (strcmp(key, "-fy") == 0) {
        if (fy) {
          opserr << "ASDAbsorbingBoundary2D ERROR: -fy flag specified twice!.\n"
                 << descr;
          return 0;
        }
        if (!get_fy())
          return 0;
      } else {
        opserr << "ASDAbsorbingBoundary2D ERROR: Invalid optional flag \""
               << key << "\".\n"
               << descr;
        return 0;
      }
    }
  }

  // done
  return new ASDAbsorbingBoundary2D(iData[0], iData[1], iData[2], iData[3],
                                    iData[4], dData[0], dData[1], dData[2],
                                    dData[3], bflag, fx, fy);
}
