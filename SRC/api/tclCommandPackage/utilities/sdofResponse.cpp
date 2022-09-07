
#include <iostream>
//
// https://portwooddigital.com/2021/02/14/how-many-clicks-does-it-take/
//

int
sdfResponse(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 9) {
    opserr << "Insufficient arguments to sdfResponse" << endln;
    return TCL_ERROR;
  }

  double m, zeta, k, Fy, alpha, dtF, dt;

  if (Tcl_GetDouble(interp, argv[1], &m) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read mass \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &zeta) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read zeta \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &k) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read k \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &Fy) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read Fy \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &alpha) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read alpha \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[6], &dtF) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read dtF \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[8], &dt) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read dt \n";
    return TCL_ERROR;
  }
  double uresidual = 0.0;
  double umaxprev = 0.0;
  if (argc > 9) {
    if (Tcl_GetDouble(interp, argv[9], &uresidual) != TCL_OK) {
      opserr << "WARNING sdfResponse -- could not read uresidual \n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &umaxprev) != TCL_OK) {
      opserr << "WARNING sdfResponse -- could not read umaxprev \n";
      return TCL_ERROR;
    }
  }

  double gamma = 0.5;
  double beta = 0.25;
  double tol = 1.0e-8;
  int maxIter = 10;

  std::ifstream infile(argv[7]);

  double c = zeta * 2 * sqrt(k * m);
  double Hkin = alpha / (1.0 - alpha) * k;

  double p0 = 0.0;
  double u0 = uresidual;
  double v0 = 0.0;
  double fs0 = 0.0;
  double a0 = (p0 - c * v0 - fs0) / m;

  double a1 = m / (beta * dt * dt) + (gamma / (beta * dt)) * c;
  double a2 = m / (beta * dt) + (gamma / beta - 1.0) * c;
  double a3 = (0.5 / beta - 1.0) * m + dt * (0.5 * gamma / beta - 1.0) * c;

  double au = 1.0 / (beta * dt * dt);
  double av = 1.0 / (beta * dt);
  double aa = 0.5 / beta - 1.0;

  double vu = gamma / (beta * dt);
  double vv = 1.0 - gamma / beta;
  double va = dt * (1 - 0.5 * gamma / beta);

  double kT0 = k;

  double umax = fabs(umaxprev);
  double amax = 0.0;
  double tamax = 0.0;
  double up = uresidual;
  double up0 = up;
  int i = 0;
  double ft, u, du, v, a, fs, zs, ftrial, kT, kTeff, dg, phat, R, R0, accel;
  while (infile >> ft) {
    i++;

    u = u0;

    fs = fs0;
    kT = kT0;
    up = up0;

    phat = ft + a1 * u0 + a2 * v0 + a3 * a0;

    R = phat - fs - a1 * u;
    R0 = R;
    if (R0 == 0.0) {
      R0 = 1.0;
    }

    int iter = 0;

    while (iter < maxIter && fabs(R / R0) > tol) {
      iter++;

      kTeff = kT + a1;

      du = R / kTeff;

      u = u + du;

      fs = k * (u - up0);
      zs = fs - Hkin * up0;
      ftrial = fabs(zs) - Fy;
      if (ftrial > 0) {
        dg = ftrial / (k + Hkin);
        if (fs < 0) {
          fs = fs + dg * k;
          up = up0 - dg;
        } else {
          fs = fs - dg * k;
          up = up0 + dg;
        }
        kT = k * Hkin / (k + Hkin);
      } else {
        kT = k;
      }

      R = phat - fs - a1 * u;
    }

    v = vu * (u - u0) + vv * v0 + va * a0;
    a = au * (u - u0) - av * v0 - aa * a0;

    u0 = u;
    v0 = v;
    a0 = a;
    fs0 = fs;
    kT0 = kT;
    up0 = up;

    if (fabs(u) > umax) {
      umax = fabs(u);
    }
    if (fabs(a) > amax) {
      amax = fabs(a);
      tamax = iter * dt;
    }
  }

  infile.close();

  char buffer[80];
  sprintf(buffer, "%f %f %f %f %f", umax, u, up, amax, tamax);

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}
