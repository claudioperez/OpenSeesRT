
int
defaultUnits(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  if (argc < 7) {
    opserr << "defaultUnits - missing a unit type want: defaultUnits -Force "
              "type? -Length type? -Time type?\n";
    return -1;
  }

  const char *force = 0;
  const char *length = 0;
  const char *time = 0;
  const char *temperature = "N/A";

  int count = 1;
  while (count < argc) {
    if ((strcmp(argv[count], "-force") == 0) ||
        (strcmp(argv[count], "-Force") == 0) ||
        (strcmp(argv[count], "-FORCE") == 0)) {
      force = argv[count + 1];
    } else if ((strcmp(argv[count], "-length") == 0) ||
               (strcmp(argv[count], "-Length") == 0) ||
               (strcmp(argv[count], "-LENGTH") == 0)) {
      length = argv[count + 1];
    } else if ((strcmp(argv[count], "-time") == 0) ||
               (strcmp(argv[count], "-Time") == 0) ||
               (strcmp(argv[count], "-TIME") == 0)) {
      time = argv[count + 1];
    } else if ((strcmp(argv[count], "-temperature") == 0) ||
               (strcmp(argv[count], "-Temperature") == 0) ||
               (strcmp(argv[count], "-TEMPERATURE") == 0) ||
               (strcmp(argv[count], "-temp") == 0) ||
               (strcmp(argv[count], "-Temp") == 0) ||
               (strcmp(argv[count], "-TEMP") == 0)) {
      temperature = argv[count + 1];
    } else {
      opserr << "defaultUnits - unrecognized unit: " << argv[count]
             << " want: defaultUnits -Force type? -Length type? -Time type?\n";
      return -1;
    }
    count += 2;
  }

  if (length == 0 || force == 0 || time == 0) {
    opserr << "defaultUnits - missing a unit type want: defaultUnits -Force "
              "type? -Length type? -Time type?\n";
    return -1;
  }

  double lb, kip, n, kn, mn, kgf, tonf;
  double in, ft, mm, cm, m;
  double sec, msec;

  if ((strcmp(force, "lb") == 0) || (strcmp(force, "lbs") == 0)) {
    lb = 1.0;
  } else if ((strcmp(force, "kip") == 0) || (strcmp(force, "kips") == 0)) {
    lb = 0.001;
  } else if ((strcmp(force, "N") == 0)) {
    lb = 4.4482216152605;
  } else if ((strcmp(force, "kN") == 0) || (strcmp(force, "KN") == 0) ||
             (strcmp(force, "kn") == 0)) {
    lb = 0.0044482216152605;
  } else if ((strcmp(force, "mN") == 0) || (strcmp(force, "MN") == 0) ||
             (strcmp(force, "mn") == 0)) {
    lb = 0.0000044482216152605;
  } else if ((strcmp(force, "kgf") == 0)) {
    lb = 4.4482216152605 / 9.80665;
  } else if ((strcmp(force, "tonf") == 0)) {
    lb = 4.4482216152605 / 9.80665 / 1000.0;
  } else {
    lb = 1.0;
    opserr << "defaultUnits - unknown force type, valid options: lb, kip, N, "
              "kN, MN, kgf, tonf\n";
    return TCL_ERROR;
  }

  if ((strcmp(length, "in") == 0) || (strcmp(length, "inch") == 0)) {
    in = 1.0;
  } else if ((strcmp(length, "ft") == 0) || (strcmp(length, "feet") == 0)) {
    in = 1.0 / 12.0;
  } else if ((strcmp(length, "mm") == 0)) {
    in = 25.4;
  } else if ((strcmp(length, "cm") == 0)) {
    in = 2.54;
  } else if ((strcmp(length, "m") == 0)) {
    in = 0.0254;
  } else {
    in = 1.0;
    opserr << "defaultUnits - unknown length type, valid options: in, ft, mm, "
              "cm, m\n";
    return TCL_ERROR;
  }

  if ((strcmp(time, "sec") == 0) || (strcmp(time, "Sec") == 0)) {
    sec = 1.0;
  } else if ((strcmp(time, "msec") == 0) || (strcmp(time, "mSec") == 0)) {
    sec = 1000.0;
  } else {
    sec = 1.0;
    opserr << "defaultUnits - unknown time type, valid options: sec, msec\n";
    return TCL_ERROR;
  }

  kip = lb / 0.001;
  n = lb / 4.4482216152605;
  kn = lb / 0.0044482216152605;
  mn = lb / 0.0000044482216152605;
  kgf = lb / (4.4482216152605 / 9.80665);
  tonf = lb / (4.4482216152605 / 9.80665 / 1000.0);

  ft = in * 12.0;
  mm = in / 25.4;
  cm = in / 2.54;
  m = in / 0.0254;

  msec = sec * 0.001;

  char string[50];

  sprintf(string, "set lb %.18e", lb);
  Tcl_Eval(interp, string);
  sprintf(string, "set lbf %.18e", lb);
  Tcl_Eval(interp, string);
  sprintf(string, "set kip %.18e", kip);
  Tcl_Eval(interp, string);
  sprintf(string, "set N %.18e", n);
  Tcl_Eval(interp, string);
  sprintf(string, "set kN %.18e", kn);
  Tcl_Eval(interp, string);
  sprintf(string, "set Newton %.18e", n);
  Tcl_Eval(interp, string);
  sprintf(string, "set kNewton %.18e", kn);
  Tcl_Eval(interp, string);
  sprintf(string, "set MN %.18e", mn);
  Tcl_Eval(interp, string);
  sprintf(string, "set kgf %.18e", kgf);
  Tcl_Eval(interp, string);
  sprintf(string, "set tonf %.18e", tonf);
  Tcl_Eval(interp, string);

  sprintf(string, "set in %.18e", in);
  Tcl_Eval(interp, string);
  sprintf(string, "set inch %.18e", in);
  Tcl_Eval(interp, string);
  sprintf(string, "set ft %.18e", ft);
  Tcl_Eval(interp, string);
  sprintf(string, "set mm %.18e", mm);
  Tcl_Eval(interp, string);
  sprintf(string, "set cm %.18e", cm);
  Tcl_Eval(interp, string);
  sprintf(string, "set m  %.18e", m);
  Tcl_Eval(interp, string);
  sprintf(string, "set meter  %.18e", m);
  Tcl_Eval(interp, string);

  sprintf(string, "set sec %.18e", sec);
  Tcl_Eval(interp, string);
  sprintf(string, "set msec %.18e", msec);
  Tcl_Eval(interp, string);

  double g = 32.174049 * ft / (sec * sec);
  sprintf(string, "set g %.18e", g);
  Tcl_Eval(interp, string);
  sprintf(string, "set kg %.18e", n * sec * sec / m);
  Tcl_Eval(interp, string);
  sprintf(string, "set Mg %.18e", 1e3 * n * sec * sec / m);
  Tcl_Eval(interp, string);
  sprintf(string, "set slug %.18e", lb * sec * sec / ft);
  Tcl_Eval(interp, string);
  sprintf(string, "set Pa %.18e", n / (m * m));
  Tcl_Eval(interp, string);
  sprintf(string, "set kPa %.18e", 1e3 * n / (m * m));
  Tcl_Eval(interp, string);
  sprintf(string, "set MPa %.18e", 1e6 * n / (m * m));
  Tcl_Eval(interp, string);
  sprintf(string, "set psi %.18e", lb / (in * in));
  Tcl_Eval(interp, string);
  sprintf(string, "set ksi %.18e", kip / (in * in));
  Tcl_Eval(interp, string);
  sprintf(string, "set psf %.18e", lb / (ft * ft));
  Tcl_Eval(interp, string);
  sprintf(string, "set ksf %.18e", kip / (ft * ft));
  Tcl_Eval(interp, string);
  sprintf(string, "set pcf %.18e", lb / (ft * ft * ft));
  Tcl_Eval(interp, string);
  sprintf(string, "set in2 %.18e", in * in);
  Tcl_Eval(interp, string);
  sprintf(string, "set ft2 %.18e", ft * ft);
  Tcl_Eval(interp, string);
  sprintf(string, "set mm2 %.18e", mm * mm);
  Tcl_Eval(interp, string);
  sprintf(string, "set cm2 %.18e", cm * cm);
  Tcl_Eval(interp, string);
  sprintf(string, "set m2 %.18e", m * m);
  Tcl_Eval(interp, string);
  sprintf(string, "set in4 %.18e", in * in * in * in);
  Tcl_Eval(interp, string);
  sprintf(string, "set ft4 %.18e", ft * ft * ft * ft);
  Tcl_Eval(interp, string);
  sprintf(string, "set mm4 %.18e", mm * mm * mm * mm);
  Tcl_Eval(interp, string);
  sprintf(string, "set cm4 %.18e", cm * cm * cm * cm);
  Tcl_Eval(interp, string);
  sprintf(string, "set m4 %.18e", m * m * m * m);
  Tcl_Eval(interp, string);
  sprintf(string, "set pi %.18e", 2.0 * asin(1.0));
  Tcl_Eval(interp, string);
  sprintf(string, "set PI %.18e", 2.0 * asin(1.0));
  Tcl_Eval(interp, string);

  int res = simulationInfo.setForceUnit(force);
  res += simulationInfo.setLengthUnit(length);
  res += simulationInfo.setTimeUnit(time);
  res += simulationInfo.setTemperatureUnit(temperature);

  return res;
}


