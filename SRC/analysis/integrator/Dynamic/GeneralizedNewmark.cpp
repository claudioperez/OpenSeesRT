//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation of the GeneralizedNewmark class.
//
// Written : cmp
// Created : 06/2024
// Revision: A
//
#include <stdexcept>
#include <GeneralizedNewmark.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <NodeIter.h>
#include <Domain.h>
// for sensitivity
#include <Node.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Parameter.h>
#include <ParameterIter.h>


GeneralizedNewmark::GeneralizedNewmark(double _gamma, double _beta, 
                                       double alphaF, double alphaM,
                                       int uFlag, int iFlag, bool aflag)
    : TransientIntegrator(0),
      gamma(_gamma), beta(_beta), 
      alphaF(1.0), alphaM(1.0), 
      unknown(uFlag), unknown_initialize(iFlag),
      step(0), dt(0.0),
      c1(0.0), c2(0.0), c3(0.0), 
      Ut(nullptr), Utdot(nullptr), Utdotdot(nullptr), 
      U(nullptr),  Udot(nullptr),  Udotdot(nullptr),
      Ua(nullptr), Uadot(nullptr), Uadotdot(nullptr), 
      determiningMass(false),
      sensitivityFlag(0), gradNumber(0), 
      massMatrixMultiplicator(0),
      dampingMatrixMultiplicator(0), assemblyFlag(aflag), independentRHS(),
      dUn(), dVn(), dAn()
{

}

GeneralizedNewmark::~GeneralizedNewmark()
{
    // clean up the memory created
    if (Ut != nullptr)
        delete Ut;
    if (Utdot != nullptr)
        delete Utdot;
    if (Utdotdot != nullptr)
        delete Utdotdot;
    if (Ua != nullptr)
        delete Ua;
    if (Uadot != nullptr)
        delete Uadot;
    if (Uadotdot != nullptr)
        delete Uadotdot;
    if (U != nullptr)
        delete U;
    if (Udot != nullptr)
        delete Udot;
    if (Udotdot != nullptr)
        delete Udotdot;

    // clean up sensitivity
    if (massMatrixMultiplicator != nullptr)
      delete massMatrixMultiplicator;
    
    if (dampingMatrixMultiplicator != nullptr)
      delete dampingMatrixMultiplicator;
}


int
GeneralizedNewmark::newStep(double deltaT)
{
    if (deltaT <= 0.0)  {
        opserr << "GeneralizedNewmark::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;  
    }
    
    if (U == nullptr)  {
      throw std::invalid_argument( "domainChange failed or not called");
      return -3;
    }

    // mark step as bootstrap or not
    if ( deltaT != dt )
        step = 0;
    else
        step++;

    dt = deltaT;

    // Set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;        
    (*Utdot) = *Udot;  
    (*Utdotdot) = *Udotdot;

    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    
    // set the constants
    switch (unknown) {
    case Displacement:
      if (beta == 0)  {
          opserr << "GeneralizedNewmark::newStep() - error in variable\n";
          opserr << "gamma = " << gamma << " beta = " << beta << endln;
          return -1;
      }
      c1 = 1.0;
      c2 = gamma/(beta*deltaT);
      c3 = 1.0/(beta*deltaT*deltaT);
      break;

    case Velocity:
      if (gamma == 0)  {
          opserr << "GeneralizedNewmark::newStep() - error in variable\n";
          opserr << "gamma = " << gamma << " beta = " << beta << endln;
          return -1;
      }
      c1 = deltaT*beta/gamma;
      c2 = 1.0;
      c3 = 1.0/(gamma*deltaT);
      break;

    case Acceleration:
      c1 = beta*deltaT*deltaT;
      c2 = gamma*deltaT;
      c3 = 1.0;
      break;
    }

    //
    // Set initial guesses for {u,v,a}_{t + dt}
    //
    int init = (step < 2 && unknown==Displacement) ? unknown : unknown_initialize;
//  int init = unknown_initialize;

    if (unknown == Displacement) {
      // determine new velocities and accelerations at t+deltaT
      double buu =  0.0;
      double buv =  0.0;
      double bua =  0.0;

      double bvu = -gamma/(beta*deltaT);
      double bvv = (1.0 - gamma/beta); 
      double bva = deltaT*(1.0 - 0.5*gamma/beta);

      double bau =  1/(beta*deltaT*deltaT);
      double bav = -1.0/(beta*deltaT);
      double baa =  1.0 + 0.5/beta;

      switch (init) {// unknown_initialize) {
        case Displacement:
          Udot->addVector(bvv, *Utdotdot, bva);

          Udotdot->addVector(baa-1/beta, *Utdot, bav);
//        Udotdot->addVector(baa, *Utdot, bav);
          break;

        case Velocity:
          // TODO: Check
          U   ->addVector(0.0, *Ut,          buu      +c1*(    - bvu)/c2);
          U   ->addVector(1.0, *Utdot,       buv      +c1*(1.0 - bvv)/c2);
          U   ->addVector(1.0, *Utdotdot,    bua      +c1*(    - bva)/c2);

          Udotdot->addVector(0.0, *Ut,       bau      +c3*(    - bvu)/c2);
          Udotdot->addVector(1.0, *Utdot,    bav      +c3*(1.0 - bvv)/c2);
          Udotdot->addVector(1.0, *Utdotdot, baa      +c3*(    - bva)/c2);
          break;

        case Acceleration:
          U   ->addVector(0.0, *Ut,          buu      +c1*(    - bau)/c3);
          U   ->addVector(1.0, *Utdot,       buv      +c1*(    - bav)/c3);
          U   ->addVector(1.0, *Utdotdot,    bua      +c1*(1.0 - baa)/c3);

          Udot->addVector(0.0, *Ut,          bvu      +c2*(    - bau)/c3);
          Udot->addVector(1.0, *Utdot,       bvv      +c2*(    - bav)/c3);
          Udot->addVector(1.0, *Utdotdot,    bva      +c2*(1.0 - baa)/c3);
          break;

      }

    } else if (unknown == Velocity) {
      double buu = 1.0;
      double buv = -deltaT*beta/gamma*(1 - gamma/beta);
      double bua = deltaT*deltaT*beta/gamma*(gamma*0.5/beta - 1.0);

      double bvu = 0.0;
      double bvv = 0.0;
      double bva = 0.0;

      double bau = 0.0;
      double bav = -1/(gamma*deltaT);
      double baa =  1 - 1/gamma;

      switch (init) {// unknown_initialize) {
        case Displacement:
//        if (step < 2)  { // trapezoidal
//          //c1 = 1.0;
//          //c2 = 2.0/deltaT;
//          //c3 = 4.0/(deltaT*deltaT);

//            (*Udot) *= -1.0;

//            double a3 = -4.0/deltaT;
//            double a4 = -1;
//            Udotdot->addVector(a4, *Utdot, a3);
//        } else {
          Udot->addVector(0.0, *Ut,          bvu      +c2*(1.0 - buu)/c1);
          Udot->addVector(1.0, *Utdot,       bvv      +c2*(    - buv)/c1);
          Udot->addVector(1.0, *Utdotdot,    bva      +c2*(    - bua)/c1);

          // a += c3*a_{n+1}
          Udotdot->addVector(0.0, *Ut,       bau      +c3*(1.0 - buu)/c1);
          Udotdot->addVector(1.0, *Utdot,    bav      +c3*(    - buv)/c1);
          Udotdot->addVector(1.0, *Utdotdot, baa      +c3*(    - bua)/c1);
//        }
          break;

        case Velocity:
          // TODO: Check
          U   ->addVector(0.0, *Ut,          buu      +c1*(    - bvu)/c2);
          U   ->addVector(1.0, *Utdot,       buv      +c1*(1.0 - bvv)/c2);
          U   ->addVector(1.0, *Utdotdot,    bua      +c1*(    - bva)/c2);

          Udotdot->addVector(0.0, *Ut,       bau      +c3*(    - bvu)/c2);
          Udotdot->addVector(1.0, *Utdot,    bav      +c3*(1.0 - bvv)/c2);
          Udotdot->addVector(1.0, *Utdotdot, baa      +c3*(    - bva)/c2);
          break;

        case Acceleration:
          U   ->addVector(0.0, *Ut,          buu      +c1*(    - bau)/c3);
          U   ->addVector(1.0, *Utdot,       buv      +c1*(    - bav)/c3);
          U   ->addVector(1.0, *Utdotdot,    bua      +c1*(1.0 - baa)/c3);

          Udot->addVector(0.0, *Ut,          bvu      +c2*(    - bau)/c3);
          Udot->addVector(1.0, *Utdot,       bvv      +c2*(    - bav)/c3);
          Udot->addVector(1.0, *Utdotdot,    bva      +c2*(1.0 - baa)/c3);
          break;

      }

    } else {
      double buu = 1.0;
      double buv = deltaT;
      double bua = deltaT*deltaT*(0.5 - beta);

      double bvu = 0.0;
      double bvv = 1.0;
      double bva = deltaT*(1.0 - gamma);

      double bau = 0.0;
      double bav = 0.0;
      double baa = 0.0;


      // Choose how to initialize state
      switch (init) { // unknown_initialize) {
        case Displacement:
          // Initialize: U == Ut
          // implying   Da = -vc/(beta dt) - ac/(2 beta)

          // u += c1*Da
//        U->addVector(*Utdot,                   -c1*(      buv/c1));
//        U->addVector(*Utdotdot,                -c1*(1.0 + bua/c1));

          Udot->addVector(0.0, *Ut,          bvu      +c2*(1.0 - buu)/c1); // 0
          Udot->addVector(1.0, *Utdot,       bvv      +c2*(    - buv)/c1); // (beta*deltaT));
          Udot->addVector(1.0, *Utdotdot,    bva      +c2*(    - bua)/c1);

          Udotdot->addVector(0.0, *Ut,       bau      +c3*(1.0 - buu)/c1); // 0
          Udotdot->addVector(1.0, *Utdot,    bav      +c3*(    - buv)/c1);
          Udotdot->addVector(1.0, *Utdotdot, baa      +c3*(    - bua)/c1);
          break; 

        case Velocity:
          // TODO: Check
          U   ->addVector(0.0, *Ut,          buu      +c1*(    - bvu)/c2);
          U   ->addVector(1.0, *Utdot,       buv      +c1*(1.0 - bvv)/c2);
          U   ->addVector(1.0, *Utdotdot,    bua      +c1*(    - bva)/c2);


          // a += c3*a_{n+1}
          Udotdot->addVector(0.0, *Ut,       bau      +c3*(    - bvu)/c2);
          Udotdot->addVector(1.0, *Utdot,    bav      +c3*(1.0 - bvv)/c2);
          Udotdot->addVector(1.0, *Utdotdot, baa      +c3*(    - bva)/c2);
          break;

        case Acceleration:
          // implying          Da = 0
          U->addVector(buu, *Utdot,        buv);
          U->addVector(1.0, *Utdotdot,     bua + c1);

          Udot->addVector(0.0, *Ut,        bvu);
          Udot->addVector(1.0, *Utdot,     bvv);
          Udot->addVector(1.0, *Utdotdot,  bva + c2); // deltaT
          break;
      }
    }

    //
    // set the trial response quantities
    //
    // determine the displacements at t+alphaF*deltaT
    (*Ua) = *Ut;
    Ua->addVector((1.0-alphaF), *U, alphaF);

    // determine the velocities at t+alphaF*deltaT
    (*Uadot) = *Utdot;
    Uadot->addVector((1.0-alphaF), *Udot, alphaF);

    // determine the velocities at t+alphaM*deltaT
    (*Uadotdot) = *Utdotdot;
    Uadotdot->addVector((1.0-alphaM), *Udotdot, alphaM);
    
    theModel->setResponse(*Ua, *Uadot, *Uadotdot);

    //
    // increment the time and apply the load
    //
    double time = theModel->getCurrentDomainTime();
    time += alphaF*deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "GeneralizedNewmark::newStep() - failed to update the domain\n";
        return -4;
    }

    return 0;
}


int
GeneralizedNewmark::update(const Vector &deltaX)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == nullptr)  {
        opserr << "WARNING GeneralizedNewmark::update() - no AnalysisModel set\n";
        return -1;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == nullptr)  {
        opserr << "WARNING GeneralizedNewmark::update() - domainChange() failed or not called\n";
        return -2;
    }  
    
    // check deltaX is of correct size
    if (deltaX.Size() != U->Size())  {
        opserr << "WARNING GeneralizedNewmark::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaX.Size() << endln;
        return -3;
    }
    
    //  determine the response at t+deltaT
    switch (unknown) {
    case Displacement:
      (*U) += deltaX;
      Udot->addVector(1.0, deltaX, c2);
      Udotdot->addVector(1.0, deltaX, c3);
      break;

    case Velocity:
      U->addVector(1.0, deltaX, c1);
      (*Udot) += deltaX;
      Udotdot->addVector(1.0, deltaX, c3);
      break;

    case Acceleration:
      U->addVector(1.0, deltaX, c1);
      Udot->addVector(1.0, deltaX, c2);        
      (*Udotdot) += deltaX;
      break;
    }

    // determine displacement and velocity at t+alphaF*deltaT
    (*Ua) = *Ut;
    Ua->addVector((1.0-alphaF), *U, alphaF);

    (*Uadot) = *Utdot;
    Uadot->addVector((1.0-alphaF), *Udot, alphaF);

    // determine the velocities at t+alphaM*deltaT
    (*Uadotdot) = *Utdotdot;
    Uadotdot->addVector((1.0-alphaM), *Udotdot, alphaM);

    // update the response at the DOFs
    theModel->setResponse(*Ua,*Uadot,*Uadotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "GeneralizedNewmark::update - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}    


const Vector &
GeneralizedNewmark::getVel()
{
  return *Udot;
}

int
GeneralizedNewmark::revertToLastStep()
{
  // set response at t+deltaT to be that at t .. for next newStep
  if (U != nullptr)  {
    (*U) = *Ut;        
    (*Udot) = *Utdot;  
    (*Udotdot) = *Utdotdot;  
  }

  return 0;
}


int
GeneralizedNewmark::formEleTangent(FE_Element *theEle)
{
    if (determiningMass == true)
        return 0;

    theEle->zeroTangent();
    
    switch (statusFlag) {
    case CURRENT_TANGENT:
        theEle->addKtToTang(alphaF*c1);
        theEle->addCtoTang(alphaF*c2);
        theEle->addMtoTang(alphaM*c3);
        break;
    case INITIAL_TANGENT:
        theEle->addKiToTang(alphaF*c1);
        theEle->addCtoTang(alphaF*c2);
        theEle->addMtoTang(alphaM*c3);
        break;
    case HALL_TANGENT:
        theEle->addKtToTang(c1*cFactor);
        theEle->addKiToTang(c1*iFactor);
        theEle->addCtoTang(c2);
        theEle->addMtoTang(c3);
        break;
    }
    
    return 0;
}

int
GeneralizedNewmark::formNodTangent(DOF_Group *theDof)
{
    if (determiningMass == true)
        return 0;
    
    theDof->zeroTangent();
    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(alphaM*c3);
    
    return 0;
}    


int
GeneralizedNewmark::domainChanged()
{
    AnalysisModel *myModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();

    // create the new Vector objects
    if (Ut == nullptr || Ut->Size() != size)  {
        
        // delete the old
        if (Ut != nullptr)
            delete Ut;
        if (Utdot != nullptr)
            delete Utdot;
        if (Utdotdot != nullptr)
            delete Utdotdot;
        if (Ua != nullptr)
            delete Ua;
        if (Uadot != nullptr)
            delete Uadot;
        if (Uadotdot != nullptr)
            delete Utdotdot;
        if (U != nullptr)
            delete U;
        if (Udot != nullptr)
            delete Udot;
        if (Udotdot != nullptr)
            delete Udotdot;
        
        // perform the allocations
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        Ua = new Vector(size);
        Uadot = new Vector(size);
        Uadotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        dUn.resize(size); 
        dUn.Zero();
        dVn.resize(size); 
        dVn.Zero();
        dAn.resize(size); 
        dAn.Zero(); 
    }        
    
    // now go through and populate U, Udot and Udotdot by iterating through
    // the DOF_Groups and getting the last committed velocity and accel
    DOF_GrpIter &theDOFs = myModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != nullptr)  {
      const ID &id = dofPtr->getID();
      int idSize = id.Size();
      
      int i;
      const Vector &disp = dofPtr->getCommittedDisp();  
      for (i=0; i < idSize; i++)  {
          int loc = id(i);
          if (loc >= 0)  {
              (*U)(loc) = disp(i);    
          }
      }
      
      const Vector &vel = dofPtr->getCommittedVel();
      for (i=0; i < idSize; i++)  {
          int loc = id(i);
          if (loc >= 0) {
              (*Udot)(loc) = vel(i);
          }
      }
      
      const Vector &accel = dofPtr->getCommittedAccel();  
      for (i=0; i < idSize; i++)  {
          int loc = id(i);
          if (loc >= 0) {
              (*Udotdot)(loc) = accel(i);
          }
      }

      // The remaining get**Sensitivity methods cause seg faults with Lagrange constraint
      // handler in dynamic (transient) analysis even when there is no sensitivity algorithm.
      // However, I don't think these methods need to be called in domainChanged -- MHS
      continue;
      
      const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);  
      for (i=0; i < idSize; i++) {
          int loc = id(i);
          if (loc >= 0) {
            dUn(loc) = dispSens(i);    
          }
      }

      const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
      for (i=0; i < idSize; i++) {
          int loc = id(i);
          if (loc >= 0) {
            dVn(loc) = velSens(i);
          }
      }

      const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);  
      for (i=0; i < idSize; i++) {
          int loc = id(i);
          if (loc >= 0) {
            dAn(loc) = accelSens(i);
          }
      }
    }    
    
    return 0;
}


int
GeneralizedNewmark::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(3);
    data(0) = gamma;
    data(1) = beta;
    data(2) = unknown;

    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING GeneralizedNewmark::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int
GeneralizedNewmark::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(3);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING GeneralizedNewmark::recvSelf() - could not receive data\n";
        gamma = 0.5;
        beta = 0.25; 
        return -1;
    }
    
    gamma  = data(0);
    beta   = data(1);
    unknown  = data(2);

    return 0;
}


void
GeneralizedNewmark::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != nullptr) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t GeneralizedNewmark - currentTime: " << currentTime;
    }
    s << "\t gamma: " << gamma << "  beta: " << beta << "\n";
    s << "\t alphaF: " << alphaF << "  alphaM: " << alphaM << "\n";
    s << "\t unknown: " << unknown << "  initialization: " << unknown_initialize << "\n";
}



int
GeneralizedNewmark::revertToStart()
{
    if (Ut != nullptr) 
        Ut->Zero();
    if (Utdot != nullptr) 
        Utdot->Zero();
    if (Utdotdot != nullptr) 
        Utdotdot->Zero();
    if (U != nullptr) 
        U->Zero();
    if (Udot != nullptr) 
        Udot->Zero();
    if (Udotdot != nullptr) 
        Udotdot->Zero();
    
    return 0;
}

int
GeneralizedNewmark::formEleResidual(FE_Element* theEle)
{
  if (sensitivityFlag == 0) {  // no sensitivity
      this->TransientIntegrator::formEleResidual(theEle);

  } else {
  
      theEle->zeroResidual();

      // Compute the time-stepping parameters on the form
      // udotdot = a1*ui+1 + a2*ui + a3*udoti + a4*udotdoti
      // udot    = a5*ui+1 + a6*ui + a7*udoti + a8*udotdoti
      // (see p. 166 of Chopra)

      // The constants are:
      // a1 = 1.0/(beta*dt*dt)
      // a2 = -1.0/(beta*dt*dt)
      // a3 = -1.0/beta*dt
      // a4 = 1.0 - 1.0/(2.0*beta)
      // a5 = gamma/(beta*dt)
      // a6 = -gamma/(beta*dt)
      // a7 = 1.0 - gamma/beta
      // a8 = 1.0 - gamma/(2.0*beta)

      // We can make use of the data members c2 and c3 of this class. 
      // As long as disp==true, they are defined as:
      // c2 = gamma/(beta*dt)
      // c3 = 1.0/(beta*dt*dt)

      // So, the constants can be computed as follows:
      if (unknown != Displacement) {
          opserr << "ERROR: GeneralizedNewmark::formEleResidual() -- the implemented"
           << " scheme only works if the displ variable is set to true." << endln;
      }

      double a2 = -c3;
      double a3 = -c2/gamma;
      double a4 = 1.0 - 1.0/(2.0*beta);
      double a6 = -c2;
      double a7 = 1.0 - gamma/beta;
      double dt = gamma/(beta*c2);
      double a8 = dt*(1.0 - gamma/(2.0*beta));

      // Pre-compute the vectors involving a2, a3, etc.
      //Vector tmp1 = V*a2 + Vdot*a3 + Vdotdot*a4;
      int vectorSize = U->Size();
      Vector dUn(vectorSize);
      Vector dVn(vectorSize);
      Vector dAn(vectorSize);
      int loc;

      AnalysisModel *myModel = this->getAnalysisModel();
      DOF_GrpIter &theDOFs = myModel->getDOFs();
      DOF_Group *dofPtr;
      while ((dofPtr = theDOFs()) != nullptr) {

        const ID &id = dofPtr->getID();
        int idSize = id.Size();
        const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
        for (int i = 0; i < idSize; i++) {
          loc = id(i);
          if (loc >= 0) {
            dUn(loc) = dispSens(i);
          }
        }

        const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
        for (int i = 0; i < idSize; i++) {
          loc = id(i);
          if (loc >= 0) {
            dVn(loc) = velSens(i);
          }
        }

        const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
        for (int i = 0; i < idSize; i++) {
          loc = id(i);
          if (loc >= 0) {
            dAn(loc) = accelSens(i);
          }
        }
      }

      // Pre-compute the vectors involving a2, a3, etc.
      // Vector tmp1 = V*a2 + Vdot*a3 + Vdotdot*a4;
      Vector tmp1(vectorSize);
      tmp1.addVector(0.0, dUn, a2);
      tmp1.addVector(1.0, dVn, a3);
      tmp1.addVector(1.0, dAn, a4);
      //Vector tmp2 = V*a6 + Vdot*a7 + Vdotdot*a8;
      Vector tmp2(vectorSize);
      tmp2.addVector(0.0, dUn, a6);
      tmp2.addVector(1.0, dVn, a7);
      tmp2.addVector(1.0, dAn, a8);

      if (massMatrixMultiplicator == 0)
          massMatrixMultiplicator = new Vector(tmp1.Size());

      if (dampingMatrixMultiplicator == 0)
          dampingMatrixMultiplicator = new Vector(tmp2.Size());

      (*massMatrixMultiplicator) = tmp1;
      (*dampingMatrixMultiplicator) = tmp2;


      // Now we're ready to make calls to the FE Element:

      // The term -dPint/dh|u fixed
      theEle->addResistingForceSensitivity(gradNumber); 

      // The term -dM/dh*acc
      theEle->addM_ForceSensitivity(gradNumber, *Udotdot, -1.0);

      // The term -M*(a2*v + a3*vdot + a4*vdotdot)
      theEle->addM_Force(*massMatrixMultiplicator,-1.0);

      // The term -C*(a6*v + a7*vdot + a8*vdotdot)
      theEle->addD_Force(*dampingMatrixMultiplicator,-1.0);

      // The term -dC/dh*vel
      theEle->addD_ForceSensitivity(gradNumber, *Udot,-1.0);
        
    }

    return 0;
}

int
GeneralizedNewmark::formNodUnbalance(DOF_Group *theDof)
{

    if (sensitivityFlag == 0) {

      this->TransientIntegrator::formNodUnbalance(theDof);

    }
    else {
      // Assemble sensitivity residual

      theDof->zeroUnbalance();

      // The term -M*(a2*v + a3*vdot + a4*vdotdot)
      theDof->addM_Force(*massMatrixMultiplicator,-1.0);

      // The term -dM/dh*acc
      theDof->addM_ForceSensitivity(*Udotdot, -1.0);

      // The term -C*(a6*v + a7*vdot + a8*vdotdot)
      theDof->addD_Force(*dampingMatrixMultiplicator,-1.0);

      // The term -dC/dh*vel
      theDof->addD_ForceSensitivity(*Udot,-1.0);

      // In case of random loads (have already been formed by 'applyLoadSensitivity')
      theDof->addPtoUnbalance();

    }


    return 0;
}

int 
GeneralizedNewmark::formSensitivityRHS(int passedGradNumber)
{
    // Set a couple of data members
    sensitivityFlag = 1;
    gradNumber = passedGradNumber;


    LinearSOE *theSOE = this->getLinearSOE();

    // Possibly set the independent part of the RHS
    if (assemblyFlag != 0)
      theSOE->setB(independentRHS);

    // Get the analysis model
    AnalysisModel *theModel = this->getAnalysisModel();

    //
    // Randomness in external load (including randomness in time series)
    //

    Domain *theDomain = theModel->getDomainPtr();

    // Loop through nodes to zero the unbalaced load
    Node *nodePtr;
    NodeIter &theNodeIter = theDomain->getNodes();
    while ((nodePtr = theNodeIter()) != nullptr)
      nodePtr->zeroUnbalancedLoad();

    // Loop through load patterns to add external load sensitivity
    LoadPattern *loadPatternPtr;
    LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    double time;
    while((loadPatternPtr = thePatterns()) != nullptr) {
      time = theDomain->getCurrentTime();
      loadPatternPtr->applyLoadSensitivity(time);
    }

    // Randomness in element/material contributions
    // Loop through FE elements
    FE_Element *elePtr;
    FE_EleIter &theEles = theModel->getFEs();    
    while((elePtr = theEles()) != nullptr) {
      theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
    }

    // Loop through DOF groups (IT IS IMPORTANT THAT THIS IS DONE LAST!)
    DOF_Group *dofPtr;
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    while((dofPtr = theDOFs()) != nullptr) {
      theSOE->addB(  dofPtr->getUnbalance(this),  dofPtr->getID()  );
    }

    // Reset the sensitivity flag
    sensitivityFlag = 0;

    return 0;
}

int 
GeneralizedNewmark::formIndependentSensitivityRHS()
{
    // For now; don't use this
/*
  sensitivityFlag = 2; // Tell subsequent methods what to be assembled

  // Get pointer to the SOE
  LinearSOE *theSOE = this->getLinearSOEPtr();


  // Get the analysis model
  AnalysisModel *theModel = this->getAnalysisModelPtr();

  
  // Loop through FE elements
  FE_Element *elePtr;
  FE_EleIter &theEles = theModel->getFEs();    
  while((elePtr = theEles()) != 0) {
  theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
  }


  // Loop through DOF groups (IT IS IMPORTANT THAT THIS IS DONE LAST!)
  DOF_Group *dofPtr;
  DOF_GrpIter &theDOFs = theModel->getDOFs();
  while((dofPtr = theDOFs()) != 0) {
  theSOE->addB(  dofPtr->getUnbalance(this),  dofPtr->getID()  );
  }


  // Set the data member of this class
  independentRHS = theSOE->getB();


  // Reset the sensitivity flag
  sensitivityFlag = 0;
*/

    return 0;
}

int 
GeneralizedNewmark::saveSensitivity(const Vector & vNew,int gradNum,int numGrads)
{

    // Compute GeneralizedNewmark parameters in general notation
    double a1 = c3;
    double a2 = -c3;
    double a3 = -c2/gamma;
    double a4 = 1.0 - 1.0/(2.0*beta);
    double a5 = c2;
    double a6 = -c2;
    double a7 = 1.0 - gamma/beta;
    double dt = gamma/(beta*c2);
    double a8 = dt*(1.0 - gamma/(2.0*beta));



  // Obtain sensitivity vectors from previous step modified by lei July 2018
  int vectorSize = U->Size();
  Vector dUn(vectorSize);
  Vector dVn(vectorSize);
  Vector dAn(vectorSize);
  int i, loc;

  AnalysisModel *myModel = this->getAnalysisModel();
  DOF_GrpIter &theDOFs = myModel->getDOFs();
  DOF_Group *dofPtr;
  while ((dofPtr = theDOFs()) != 0) {

    const ID &id = dofPtr->getID();
    int idSize = id.Size();
    const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
    for (i = 0; i < idSize; i++) {
      loc = id(i);
      if (loc >= 0) {
        dUn(loc) = dispSens(i);
      }
    }

    const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
    for (i = 0; i < idSize; i++) {
      loc = id(i);
      if (loc >= 0) {
        dVn(loc) = velSens(i);
      }
    }

    const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
    for (i = 0; i < idSize; i++) {
      loc = id(i);
      if (loc >= 0) {
        dAn(loc) = accelSens(i);
      }
    }
  }



    // Compute new acceleration and velocity vectors:
    Vector vdotNew(vectorSize);
    Vector vdotdotNew(vectorSize);
    //(*vdotdotNewPtr) = vNew*a1 + V*a2 + Vdot*a3 + Vdotdot*a4;
    vdotdotNew.addVector(0.0, vNew, a1);
    vdotdotNew.addVector(1.0, dUn, a2);
    vdotdotNew.addVector(1.0, dVn, a3);
    vdotdotNew.addVector(1.0, dAn, a4);
    
    //(*vdotNewPtr) = vNew*a5 + V*a6 + Vdot*a7 + Vdotdot*a8;
    vdotNew.addVector(0.0, vNew, a5);
    vdotNew.addVector(1.0, dUn, a6);
    vdotNew.addVector(1.0, dVn, a7);
    vdotNew.addVector(1.0, dAn, a8);

    // update
    dUn = vNew;
    dVn = vdotNew;
    dAn = vdotdotNew;

    // Now we can save vNew, vdotNew and vdotdotNew
    //AnalysisModel *myModel = this->getAnalysisModel();
    DOF_GrpIter &theDOFGrps = myModel->getDOFs();
    DOF_Group   *dofPtr1;
    while ( (dofPtr1 = theDOFGrps() ) != 0)  {
  dofPtr1->saveSensitivity(vNew,vdotNew,vdotdotNew,gradNum,numGrads);
    }
  
    return 0;
}

int 
GeneralizedNewmark::commitSensitivity(int gradNum, int numGrads)
{
    // Loop through the FE_Elements and set unconditional sensitivities
    AnalysisModel *theAnalysisModel = this->getAnalysisModel();
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != nullptr) {
      elePtr->commitSensitivity(gradNum, numGrads);
    }

    return 0;
}


double
GeneralizedNewmark::getCFactor() 
{
  return c2;
}


int 
GeneralizedNewmark::computeSensitivities()
{

  LinearSOE *theSOE = this->getLinearSOE();

  // Zero out the old right-hand side of the SOE
  theSOE->zeroB();
  
  // Form the part of the RHS which are indepent of parameter
  this->formIndependentSensitivityRHS();
  AnalysisModel *theModel = this->getAnalysisModel();  //Abbas 
  Domain *theDomain=theModel->getDomainPtr();//Abbas
  ParameterIter &paramIter = theDomain->getParameters();

  Parameter *theParam;
  // De-activate all parameters
  while ((theParam = paramIter()) != nullptr)
    theParam->activate(false);
  
  // Now, compute sensitivity wrt each parameter
  int numGrads = theDomain->getNumParameters();

  paramIter = theDomain->getParameters();
  
  while ((theParam = paramIter()) != nullptr) {
    
    // Activate this parameter
    theParam->activate(true);
    
    // Zero the RHS vector
    theSOE->zeroB();
    
    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();

    // Form the RHS
    this->formSensitivityRHS(gradIndex);
    
    // Solve for displacement sensitivity 
    theSOE->solve();

    // Save sensitivity to nodes
    this->saveSensitivity( theSOE->getX(), gradIndex, numGrads );
    

    // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
    this->commitSensitivity(gradIndex, numGrads);
    
    // De-activate this parameter for next parameter sensitivity
    theParam->activate(false);
  }
  
  return 0;
}

