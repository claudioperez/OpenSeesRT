#include <elementAPI.h>
#include "BWBN.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>

static inline double 
signum(double value)
{
  if (value > 0.0)
    return  1.0;
  else
    return -1.0;
}


void * OPS_ADD_RUNTIME_VPV(OPS_BWBN)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData1[1];
  double dData[13];
  int    iData2[1];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData1) != 0) {
    opserr << "WARNING invalid uniaxialMaterial BWBN tag" << endln;
    return 0;
  }

  numData = 13;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Double Values\n";
    return 0;    
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData2) != 0) {
    opserr << "WARNING invalid maxNumIter" << endln;
    return 0;
  }

  theMaterial = new BWBN(iData1[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
             dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],iData2[0]);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type BWBN\n";
    return 0;
  }

  return theMaterial;
}



BWBN::BWBN(int tag, 
       double p_alpha,
       double p_ko,
       double p_n,
       double p_gamma,
       double p_beta,
       double p_Ao,
       double p_q,
       double p_zetas,
       double p_p,
       double p_Shi,
       double p_deltaShi,
       double p_lamda,
       double ptolerance,
       int pMaxNumIter)
 :UniaxialMaterial(tag,MAT_TAG_BWBN),
  alpha(p_alpha), ko(p_ko), n(p_n), gamma(p_gamma), beta(p_beta), Ao(p_Ao), q(p_q), 
  zetas(p_zetas), p(p_p), Shi(p_Shi), deltaShi(p_deltaShi), lamda(p_lamda), tolerance(ptolerance),
  maxNumIter(pMaxNumIter)
{
  // Initialize variables
  this->revertToStart();
}


BWBN::~BWBN()
{
  // does nothing
}


int 
BWBN::setTrialStrain (double strain, double strainRate)
{
    // Set trial strain and compute strain increment
    Tstrain = strain;
    const double dStrain = Tstrain - Cstrain;
    const double sgn  = signum(dStrain);

    // Newton-Raphson scheme to solve for z_{i+1} := z1
    double startPoint = 0.01;
    Tz = startPoint;

    double Tzold = startPoint;
    double Tznew = 1.0;  
    int    count = 0;
    while ( ( fabs(Tzold-Tznew) > tolerance ) && count<maxNumIter) {
        Te = Ce + (1.0-alpha)*ko*dStrain*Tz;
        double Tzeta1 = zetas*(1.0 - exp(-p*Te));
        double Tzeta2 = (Shi+deltaShi*Te)*(lamda+Tzeta1);
        double zu = pow(1/(beta+gamma),1/n);
        double h = 1.0 - Tzeta1*exp(-pow(Tz*sgn - q*zu, 2)/(Tzeta2*Tzeta2));
        double Psi = gamma + beta*signum(dStrain*Tz);
        double Phi = Ao - pow(fabs(Tz),n)*Psi;
        double f = Tz - Cz - Phi*h*dStrain;

        // Evaluate function derivative f' (underscore:=prime)
        double Te_ = (1.0-alpha)*ko*dStrain;
        double Tzeta1_ = zetas*p*exp(-p*Te)*Te_;
        double Tzeta2_ = Shi*Tzeta1_ + lamda*deltaShi*Te_ + deltaShi*Te*Tzeta1_ + deltaShi*Te_*Tzeta1;
        //h_ = -exp(-pow(Tz*sgn-q*zu,2)/(Tzeta2*Tzeta2))*(Tzeta1_-Tzeta1*2*(Tz*sgn-q*zu)/(Tzeta2*Tzeta2)+Tzeta2_*Tzeta1*2*pow((Tz*sgn-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2));
        double h_ = -exp(-pow(Tz*sgn-q*zu,2)/(Tzeta2*Tzeta2))*(
                 Tzeta1_-Tzeta1*2*(Tz*sgn-q*zu)*sgn/(Tzeta2*Tzeta2)
                +Tzeta2_*Tzeta1*2*pow((Tz*sgn-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2)
              );

        double pow1 = (Tz==0.0)? 0.0 : pow(fabs(Tz),(n-1));
        double Phi_ = - n*pow1*signum(Tz)*Psi;
        double f_ = 1.0 - (Phi_*h+Phi*h_)*dStrain;


        // Issue warning if derivative is zero
        if ( fabs(f_) < 1.0e-10 ) {
            opserr << "WARNING: BWBN::setTrialStrain() -- zero derivative " << endln
                << " in Newton-Raphson scheme" << endln;
        }

        // Take a Newton step
        Tznew = Tz - f/f_;

        // Update the root (but the keep the old for convergence check)
        Tzold = Tz; 
        Tz = Tznew;

        // Update counter
        count++;

        // Issue warning if we didn't converge
        if (count == maxNumIter) {
            opserr << "WARNING: BWBN::setTrialStrain() -- did not" << endln
                << " find the root z_{i+1}, after " << maxNumIter << " iterations" << endln
                << " and norm: " << fabs(Tzold-Tznew) << endln;
        }

        // Compute stress
        Tstress = alpha*ko*Tstrain + (1-alpha)*ko*Tz;

        // Compute deterioration parameters
        Te     = Ce + (1-alpha)*ko*dStrain*Tz;
        Tzeta1 = zetas*(1-exp(-p*Te));
        Tzeta2 = (Shi+deltaShi*Te)*(lamda+Tzeta1);        

        // Compute tangent
        if (Tz != 0.0) {
            Psi = gamma + beta*signum(dStrain*Tz);
            Phi = Ao - pow(fabs(Tz),n)*Psi;
            double b1 = (1-alpha)*ko*Tz;
            double b2 = zetas*p*exp(-p*Te)*b1;
            double b3 = Shi*b2+lamda*deltaShi*b1+deltaShi*Te*b2+deltaShi*b1*Tzeta1;
            double b4 = -exp(-pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2))*(b2+b3*Tzeta1*2*pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2)); 
            h = 1.0-Tzeta1*exp(-pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2));
            
            double b5 = (1.0-alpha)*ko*dStrain;
            double b6 = zetas*p*exp(-p*Te)*b5;
            double b7 = Shi*b6+lamda*deltaShi*b5+deltaShi*Te*b6+deltaShi*b5*Tzeta1;
            //b8 = -exp(-pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2))*(b6-Tzeta1*2*(Tz*signum(dStrain)-q*zu)/(Tzeta2*Tzeta2)+b7*Tzeta1*2*pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2));

            double b8 = -exp(-pow((Tz*sgn-q*zu),2)/(Tzeta2*Tzeta2))*(b6-Tzeta1*2*(Tz*sgn-q*zu)*sgn/(Tzeta2*Tzeta2)+b7*Tzeta1*2*pow((Tz*sgn-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2));
            pow1 = pow(fabs(Tz),(n-1));
            double b9 = - n*pow1*signum(Tz)*Psi;
            double DzDeps = (h*Phi-b4*Phi)/(1.0 - (b9*h+Phi*b8)*dStrain);
            Ttangent = alpha*ko + (1-alpha)*ko*DzDeps;
            //Ttangent = Tstress/Tstrain;
        } else {
            Ttangent = alpha*ko + (1-alpha)*ko;
        }

    }

    return 0;
}

double 
BWBN::getStress(void)
{
    return Tstress;
}

double 
BWBN::getInitialTangent(void)
{
  return ( alpha*ko + (1-alpha)*ko*Ao );
}


double 
BWBN::getTangent(void)
{
  return Ttangent;
}

double 
BWBN::getStrain(void)
{
  return Tstrain;
}

int 
BWBN::commitState(void)
{
  // Commit trial history variables
  Cstrain = Tstrain;
  Cz = Tz;
  Ce = Te;

  return 0;
}

int 
BWBN::revertToLastCommit(void)
{
    // Nothing to do here
    return 0;
}

int 
BWBN::revertToStart(void)
{
    Tstrain = 0.0;
    Cstrain = 0.0;
    Tz = 0.0;
    Cz = 0.0;
    Te = 0.0;
    Ce = 0.0;
    Tstress = 0.0;
    Ttangent = alpha*ko + (1-alpha)*ko*Ao;

    return 0;
}

UniaxialMaterial *
BWBN::getCopy(void)
{
    BWBN *theCopy =
    new BWBN(this->getTag(), alpha, ko, n, gamma,
                        beta, Ao, q, zetas, p, Shi, deltaShi, lamda,tolerance,maxNumIter);
        
    theCopy->Tstrain = Tstrain;
    theCopy->Cstrain = Cstrain;
    theCopy->Tz = Tz;
    theCopy->Cz = Cz;
    theCopy->Te = Te;
    theCopy->Ce = Ce;
    theCopy->Tstress = Tstress;
    theCopy->Ttangent = Ttangent;

    return theCopy;
}

int 
BWBN::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
BWBN::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

void 
BWBN::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": " << "\"BWBN\", ";
    s << "\"alpha\": " << alpha << ", ";
    s << "\"ko\": " << ko << ", ";
    s << "\"n\": " << n << ", ";
    s << "\"gamma\": " << gamma << ", ";
    s << "\"beta\": " << beta << ", ";
    s << "\"Ao\": " << Ao << ", ";
    s << "\"q\": " << q << ", ";
    s << "\"deltaA\": " << zetas << ", ";
    s << "\"deltaNu\": " << p << ", ";
    s << "\"deltaEta\": " << Shi << ", ";
    s << "\"deltaNu\": " << deltaShi << ", ";
    s << "\"deltaEta\": " << lamda ;
    s << "}";

  } else {
    s << "BWBN, tag: " << this->getTag() << endln;
    s << "  alpha: " << alpha << endln;
    s << "  ko: " << ko << endln;
    s << "  n: " << n << endln;
    s << "  gamma: " << gamma << endln;
    s << "  beta: " << beta << endln;
    s << "  Ao: " << Ao << endln;
    s << "  q: " << q << endln;
    s << "  deltaA: " << zetas << endln;
    s << "  deltaNu: " << p << endln;
    s << "  deltaEta: " << Shi << endln;
    s << "  deltaNu: " << deltaShi << endln;
    s << "  deltaEta: " << lamda << endln;
  }
}


