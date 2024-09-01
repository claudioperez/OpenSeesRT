template <int ndm>
int GPYS(const Matrix& GPYSC, const VectorND<ndm>& x, ScVec, 
         double* f, VectorND<ndm>& g, MatrixND<ndm,ndm>& h)
{
//GPYS function value, gradient and Hessian of polynomial yield surface    
//  [F,G,H] = GPYS (GPYSC,XYZ,SCVEC)
//  the function determines the value F(X,Y,Z), the gradient G(X,Y,Z), and
//  the Hessian matrix (2nd derivative) H(X,Y,Z) of F at a point XYZ
//  for a general polynomial yield surface with coefficients GPYSC
//  SCVEC is a scale vector for the variables X, Y, and Z
//   G = the gradient (normal)   of the yield surface = [dF/dX; dF/dY; dF/dZ]
//   H = the Hessian (2nd deriv) of the yield surface = dG/dXYZ
//     = [d2(F)/d(X)^2    d2(F)/d(X)d(Y) d2(F)/d(X)d(Z);
//        d2(F)/d(Y)d(X)  d2(F)/d(Y)^2   d2(F)/d(Y)d(Z);
//        d2(F)/d(Z)d(X)  d2(F)/d(Z)d(Y) d2(F)/d(Z)^2]
//
//  The coefficients of the polynomial yield surface are specified as follows
//  general 3d case
//           GPYSC = [d1 a1 b1 c1; 
//                    d2 a2 b2 c2; 
//                    d3 a3 b3 c3; 
//                           ...]
//          for  F = Sum_i (di*(X^ai)*(Y^bi)*(Z^ci))
//  general 2d case
//           GPYSC = [c1 a1 b1; c2 a2 b2; c3 a3 b3; ...]
//       for     F = Sum_i (ci*(X^ai)*(Y^bi))

//  Examples
//    for 2d         GPYSC   = [1.5 2 0  ; 1 0 3 ; 3 2 4     ; -1.2 0 0]
//              for F(X,Y)   =  1.5*X^2   + Y^3   + 3*X^2*Y^4  -1.2
//
//    for 3d         GPYSC   = [1.5 2 0 0 ; 1 0 3 0 ; 3 2 2 1    ; -1.2 0 0 0]
//              for F(X,Y,Z) =  1.5*X^2    +   Y^3 +  3*X^2*Y^2*Z  -1.2
//
//  =========================================================================================
//  function  written by Chin-Long Lee(c)                                             03-2006
//  function  extended to 3d case by Svetlana Seizovic                                11-2007
//  absolute  values for x,y,z and several sign(x) changes by Thanh Ngoc Do           08-2016
//  -----------------------------------------------------------------------------------------
//

  if (nargin>2)
    xyz = xyz./ScVec;

  // Perturbation of Hessian when x = 0 or y = 0 raised to power < 2
  if any(GPYSC(1:3,2:3) < 2) {
    xyz(xyz == 0) = 1e-6;
  }

  fscalar(GPYSC, x, f);
  gvector(GPYSC, x, g);
  hmatrix(GPYSC, x, h);

  if (nargin>2) {
    J = diag(1./ScVec);
    g = J*g;
    h = J*h*J;
  }
}

////  function fscalar -----------------------------------------------------------------------
template <int ndm> int
fscalar(const Matrix& GPYSC, VectorND<ndm> x, double* f)
{
  // yield condition function f
  *f = 0.0;
  for (int i=0; i<GPYSC.noRows(); i++) {
    double c = 1.0;
    for (int j=0; j<ndm; j++)
      c *= GPYSC(i, j)*x[j];
    *f += c;
  }
}

////  function gvector -----------------------------------------------------------------------
double* 
gvector(GPYSC,xyz)
{
  // plastic deformation vector g (gradient of f normal to yield surface)

  if (size(GPYSC,2) == 3)  //2-dim case
    double x = std::abs(xyz[0]); 
    double y = std::abs(xyz[1]); 
    // df/dx
    // =====
    // look for terms with x
    I = find(GPYSC(:,2));
    double dfdx = 0.0;
    for (auto& row : GPYSC) {
      double coeff = row[0]*row[1];
      dfdx += coeff * pow(x, row[1]-1) * pow(y,row[2]) * sign(xyz[0]);
    }
    
    // df/dy
    // =====
    // look for terms with y
    I = find(GPYSC(:,3));
    coeff = GPYSC(I,1).*GPYSC(I,3);
    dfdy = sum(coeff.*(x.^GPYSC(I,2)).*(y.^(GPYSC(I,3)-1))) * sign(xyz(2));
    
    // define g-vector
    g = [dfdx; dfdy];
    
  else if (size(GPYSC,2)==4) {//3-dim case
    x = abs(xyz(1)); 
    y = abs(xyz(2)); 
    z = abs(xyz(3));
    // df/dx
    // =====
    // look for terms with x
    I = find(GPYSC(:,2));
    coeff = GPYSC(I,1).*GPYSC(I,2);
    dfdx = sum(coeff.*(x.^(GPYSC(I,2)-1)).*(y.^GPYSC(I,3)).*(z.^GPYSC(I,4))) * sign(xyz(1));
    
    // df/dy
    // =====
    // look for terms with y
    I = find(GPYSC(:,3));
    coeff = GPYSC(I,1).*GPYSC(I,3);
    dfdy = sum(coeff.*(x.^GPYSC(I,2)).*(y.^(GPYSC(I,3)-1)).*(z.^GPYSC(I,4))) * sign(xyz(2));
    
    // df/dz
    // =====
    // look for terms with z
    I = find(GPYSC(:,4));
    coeff = GPYSC(I,1).*GPYSC(I,4);
    dfdz = sum(coeff.*(x.^GPYSC(I,2)).*(y.^GPYSC(I,3)).*(z.^(GPYSC(I,4)-1))) * sign(xyz(3));
    
    // define g-vector
    g = [dfdx; dfdy; dfdz];
  }
}

////  function hmatrix -----------------------------------------------------------------------
double* hmatrix(std::vector<VectorND<3>> &GPYSC, VectorND<2> &xyz)
{
  // H = dg/dx matrix (Hessian matrix of yield function f)
    //2-dim case
    // extract variables
    double x = std::abs(xyz[0]); 
    double y = std::abs(xyz[1]);
    
    // set up df/dx coefficients
    I = find(GPYSC(:,2));
    c = GPYSC(I,1).*GPYSC(I,2);
    GPYSCx = [c GPYSC(I,2)-1 GPYSC(I,3)]; // df/dx coefficient
    
    // set up d2f/dxdx coefficients
    I = find(GPYSCx(:,2));
    c = GPYSCx(I,1).*GPYSCx(I,2);
    GPYSCxx = [c GPYSCx(I,2)-1 GPYSCx(I,3)]; // d2f/dx2 coefficient
    
    // set up df/dy coefficients
    I = find(GPYSC(:,3));
    c = GPYSC(I,1).*GPYSC(I,3);
    GPYSCy = [c GPYSC(I,2) GPYSC(I,3)-1];
    
    // set up d2f/dydy coefficients
    I = find(GPYSCy(:,3));
    c = GPYSCy(I,1).*GPYSCy(I,3);
    GPYSCyy = [c GPYSCy(I,2) GPYSCy(I,3)-1]; // d2f/dy2 coeffs
    
    // set up d2f/dxdy coefficients
    I = find(GPYSCx(:,3));
    c = GPYSCx(I,1).*GPYSCx(I,3);
    GPYSCxy = [c GPYSCx(I,2) GPYSCx(I,3)-1];
    
    d2fdx2  = sum(GPYSCxx(:,1).*(x.^GPYSCxx(:,2)).*(y.^GPYSCxx(:,3)));
    d2fdy2  = sum(GPYSCyy(:,1).*(x.^GPYSCyy(:,2)).*(y.^GPYSCyy(:,3)));
    d2fdxdy = sum(GPYSCxy(:,1).*(x.^GPYSCxy(:,2)).*(y.^GPYSCxy(:,3))) * sign(xyz(1)) * sign(xyz(2));
    
    // define h-matrix
    h  = [d2fdx2 d2fdxdy; d2fdxdy  d2fdy2];
}

double* 
hmatrix(std::vector<VectorND<4>> &GPYSC, VectorND<3> &xyz)
{
    //3-dim case
    // extract variables
    double x = std::abs(xyz[0]); 
    double y = std::abs(xyz[1]); 
    double z = std::abs(xyz[2]);
    
    // set up df/dx coefficients
    I = find(GPYSC(:,2));
    c = GPYSC(I,1).*GPYSC(I,2);
    GPYSCx = [c GPYSC(I,2)-1 GPYSC(I,3) GPYSC(I,4)]; // df/dx coefficient
    
    // set up d2f/dxdx coefficients
    I = find(GPYSCx(:,2));
    c = GPYSCx(I,1).*GPYSCx(I,2);
    GPYSCxx = [c GPYSCx(I,2)-1 GPYSCx(I,3) GPYSCx(I,4)]; // d2f/dx2 coefficient
    
    // set up df/dy coefficients
    I = find(GPYSC(:,3));
    c = GPYSC(I,1).*GPYSC(I,3);
    GPYSCy = [c GPYSC(I,2) GPYSC(I,3)-1 GPYSC(I,4)];
    
    // set up d2f/dydy coefficients
    I = find(GPYSCy(:,3));
    c = GPYSCy(I,1).*GPYSCy(I,3);
    GPYSCyy = [c GPYSCy(I,2) GPYSCy(I,3)-1 GPYSCy(I,4)]; // d2f/dy2 coeffs
    
    // set up df/dz coefficients
    I = find(GPYSC(:,4));
    c = GPYSC(I,1).*GPYSC(I,4);
    GPYSCz = [c GPYSC(I,2) GPYSC(I,3) GPYSC(I,4)-1];
    
    // set up d2f/dzdz coefficients
    I = find(GPYSCz(:,4));
    c = GPYSCz(I,1).*GPYSCz(I,4);
    GPYSCzz = [c GPYSCz(I,2) GPYSCz(I,3) GPYSCz(I,4)-1]; // d2f/dz2 coeffs
    
    // set up d2f/dxdy coefficients
    I = find(GPYSCx(:,3));
    c = GPYSCx(I,1).*GPYSCx(I,3);
    GPYSCxy = [c GPYSCx(I,2) GPYSCx(I,3)-1 GPYSCx(I,4)];
    
    // set up d2f/dxdz coefficients
    I = find(GPYSCx(:,4));
    c = GPYSCx(I,1).*GPYSCx(I,4);
    GPYSCxz = [c GPYSCx(I,2) GPYSCx(I,3) GPYSCx(I,4)-1];
    
    // set up d2f/dydz coefficients
    I = find(GPYSCy(:,4));
    c = GPYSCy(I,1).*GPYSCy(I,4);
    GPYSCyz = [c GPYSCy(I,2) GPYSCy(I,3) GPYSCy(I,4)-1];
    
    d2fdx2  = sum(GPYSCxx(:,1).*(x.^GPYSCxx(:,2)).*(y.^GPYSCxx(:,3)).*(z.^GPYSCxx(:,4)));
    d2fdy2  = sum(GPYSCyy(:,1).*(x.^GPYSCyy(:,2)).*(y.^GPYSCyy(:,3)).*(z.^GPYSCyy(:,4)));
    d2fdz2  = sum(GPYSCzz(:,1).*(x.^GPYSCzz(:,2)).*(y.^GPYSCzz(:,3)).*(z.^GPYSCzz(:,4)));
    d2fdxdy = sum(GPYSCxy(:,1).*(x.^GPYSCxy(:,2)).*(y.^GPYSCxy(:,3)).*(z.^GPYSCxy(:,4))) * sign(xyz(1)) * sign(xyz(2));
    d2fdxdz = sum(GPYSCxz(:,1).*(x.^GPYSCxz(:,2)).*(y.^GPYSCxz(:,3)).*(z.^GPYSCxz(:,4))) * sign(xyz(1)) * sign(xyz(3));
    d2fdydz = sum(GPYSCyz(:,1).*(x.^GPYSCyz(:,2)).*(y.^GPYSCyz(:,3)).*(z.^GPYSCyz(:,4))) * sign(xyz(2)) * sign(xyz(3));
    
    // define h-matrix
    h  = [d2fdx2  d2fdxdy d2fdxdz; 
          d2fdxdy  d2fdy2 d2fdydz; 
          d2fdxdz d2fdydz d2fdz2];
}
