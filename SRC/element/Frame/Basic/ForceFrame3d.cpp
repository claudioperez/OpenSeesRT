//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
/*
 * References
 *
 *  Element State Determination Algorithm
 *  ---
 *  Neuenhofer, A. and F. C. Filippou (1997). "Evaluation of Nonlinear Frame Finite
 *  Element Models." Journal of Structural Engineering, 123(7):958-966.
 *
 *  Spacone, E., V. Ciampi, and F. C. Filippou (1996). "Mixed Formulation of
 *  Nonlinear Beam Finite Element." Computers and Structures, 58(1):71-83.
 *
 *  Response Sensitivity
 *  ---
 *  Scott, M. H., P. Franchin, G. L. Fenves, and F. C. Filippou (2004).
 *  "Response Sensitivity for Nonlinear Beam-Column Elements."
 *  Journal of Structural Engineering, 130(9):1281-1288.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <array>

#include <Information.h>
#include <Parameter.h>
#include <ForceFrame3d.h>
#include <BeamIntegration.h>
#include <FrameSection.h>
#include <interpolate/cbdi.h>
#include <Node.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <ElementIter.h>
#include <Matrix.h>

#define ELE_TAG_ForceFrame3d 0 // TODO

Matrix ForceFrame3d::theMatrix(12, 12);
Vector ForceFrame3d::theVector(12);

VectorND<ForceFrame3d::nq>                    ForceFrame3d::es_trial[maxNumSections];
MatrixND<ForceFrame3d::nq,ForceFrame3d::nq> ForceFrame3d::Fs_trial[maxNumSections];
VectorND<ForceFrame3d::nq>                    ForceFrame3d::sr_trial[maxNumSections];

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
ForceFrame3d::ForceFrame3d()
 : BasicFrame3d(0, ELE_TAG_ForceFrame3d),
   beamIntegr(nullptr),
   numSections(0), sections(nullptr),
   maxIters(0), tol(0.0),
   initialFlag(0),
   fs(nullptr), es(nullptr), Ssr(nullptr), es_save(nullptr),
   density(0), mass_flag(0),
   mass_initialized(false),
   Ki(0),
   parameterID(0)
{
  kv.zero();
  K_save.zero();
  q_save.zero();
  q_pres.zero();
}

// Constructor which takes the unique element tag, sections,
// and the node ID's of its nodal end points.
// allocates the necessary space needed by each object
ForceFrame3d::ForceFrame3d(int tag, std::array<int, 2>& nodes, int numSec,
                           FrameSection** sec, BeamIntegration& bi,
                           FrameTransform3d& coordTransf, 
                           double massDensPerUnitLength, int maxNumIters,
                           double tolerance, int mass_flag_)
 : BasicFrame3d(tag, ELE_TAG_ForceFrame3d, nodes, coordTransf),
   beamIntegr(nullptr),
   numSections(0), sections(nullptr),
   maxIters(maxNumIters), tol(tolerance),
   initialFlag(0),
   fs(nullptr), es(nullptr), Ssr(nullptr), es_save(nullptr),
   Ki(0),
   density(massDensPerUnitLength), mass_flag(mass_flag_),
   mass_initialized(false),
   parameterID(0)
{
  kv.zero();
  K_save.zero();
  q_save.zero();
  q_pres.zero();

  beamIntegr = bi.getCopy();
  if (beamIntegr == nullptr) {
    opserr << "Error: ForceFrame3d::ForceFrame3d: could not create copy of beam integration object"
           << "\n";
    exit(-1);
  }

  this->setSectionPointers(numSec, sec);
}

// Destructor
ForceFrame3d::~ForceFrame3d()
{
  if (sections != 0) {
    for (int i = 0; i < numSections; i++)
      if (sections[i] != 0)
        delete sections[i];
    delete[] sections;
  }

  if (fs != nullptr)
    delete[] fs;

  if (es != nullptr)
    delete[] es;

  if (Ssr != nullptr)
    delete[] Ssr;

  if (es_save != nullptr)
    delete[] es_save;

  if (beamIntegr != nullptr)
    delete beamIntegr;

  if (Ki != nullptr)
    delete Ki;
}


int
ForceFrame3d::setNodes()
{
  // call the DomainComponent class method
  this->BasicFrame3d::setNodes();

  double L = this->getLength(State::Init);

  if (L == 0.0)
    return -1;

  beamIntegr->getSectionLocations(numSections, L, xi);
  beamIntegr->getSectionWeights(numSections, L, wt);


  if (initialFlag == 0)
    this->initializeSectionHistoryVariables();

//this->revertToStart();

  return 0;
}

int
ForceFrame3d::getIntegral(Field field, State state, double& total)
{

  if (this->setState(State::Init) != 0)
    return -1;

  total = 0.0;
  switch (field) {

    // Integrate density to compute total mass
    case Field::Density: {
      double value = 0.0;
      for (int i=0; i< numSections; i++) {
        // First try using section's internal density
        if (sections[i]->getIntegral(Field::Density, state, value) == 0) {
          total += wt[i]*value;
        }
        // if that didnt work, just multiply by our density
        else if (sections[i]->getIntegral(Field::Unit, state, value) == 0) {
          total += wt[i]*density;
        }
        else {
          ; // TODO: This should be written to a log
        }
      }
      return 0;
    }

    case Field::PolarInertia: {
      for (int i=0; i< numSections; i++) {
        double A;
        if (sections[i]->getIntegral(Field::UnitYY, state, A) != 0)
          continue;

        // Get \int \rho y^2
        double Iz;
        if (sections[i]->getIntegral(Field::DensityYY, state, Iz) != 0) {
          // Section does not allow integrating density; try
          // integrating product of inertia and multiplying by rho
          if (sections[i]->getIntegral(Field::UnitYY, state, Iz) == 0)
            Iz *= density/A;
          else
            continue;
        }
        // Get \int \rho z^2
        double Iy;
        if (sections[i]->getIntegral(Field::DensityZZ, state, Iy) != 0) {
          if (sections[i]->getIntegral(Field::UnitZZ, state, Iy) == 0)
            Iy *= density/A;
          else
            continue;
        }
        total += wt[i]*(Iy + Iz);
      }
      return 0;
    }

    default:
      return -1;
  }
}

int
ForceFrame3d::commitState()
{
  int err = 0;
  int i   = 0;

  // Call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0)
    opserr << "ForceFrame3d::commitState () - failed in base class";

  do {

    es_save[i] = es[i];
    err         = sections[i++]->commitState();

  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // commit the transformation between coord. systems
  if ((err = theCoordTransf->commitState()) != 0)
    return err;

  // commit the element variables state
  K_save = kv;
  q_save = q_pres;

  //   initialFlag = 0;  fmk - commented out, see what happens to Example3.1.tcl if uncommented
  //                         - i have not a clue why, ask remo if he ever gets in contact with us again!

  return err;
}

int
ForceFrame3d::revertToLastCommit()
{
  int err;
  int i = 0;
  do {
    es[i] = es_save[i];
    err   = sections[i]->revertToLastCommit();

    sections[i]->setTrialState<nsr,scheme>(es[i]);
    Ssr[i] = sections[i]->getResultant<nsr,scheme>();
    fs[i]  = sections[i]->getFlexibility<nsr,scheme>();

    i++;

  } while (err == 0 && i < numSections);


  if (err)
    return err;

  // Revert the transformation to last commit
  if ((err = theCoordTransf->revertToLastCommit()) != 0)
    return err;

  // Revert the element state to last commit
  q_pres = q_save;
  kv = K_save;

  initialFlag = 0;

  return err;
}


int
ForceFrame3d::revertToStart()
{
  // revert the sections state to start
  int err;
  int i = 0;
  do {
    fs[i].zero();
    es[i].zero();
    Ssr[i].zero();
    err = sections[i++]->revertToStart();

  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // revert the transformation to start
  if ((err = theCoordTransf->revertToStart()) != 0)
    return err;

  // revert the element state to start
  q_save.zero();
  K_save.zero();

  q_pres.zero();
  kv.zero();

  initialFlag = 0;
  // this->update();
  return err;
}

VectorND<6>&
ForceFrame3d::getBasicForce()
{
  return q_pres;
}

MatrixND<6, 6>&
ForceFrame3d::getBasicTangent(State state, int rate)
{
  return kv;
}

const Matrix &
ForceFrame3d::getMass()
{
    if (!mass_initialized) {
      if (this->getIntegral(Field::Density, State::Init, total_mass) != 0)
        ;
      if (this->getIntegral(Field::PolarInertia, State::Init, twist_mass) != 0)
        ;
      mass_initialized = true;
    }

    if (total_mass == 0.0) {

        thread_local MatrixND<12,12> M{0.0};
        thread_local Matrix Wrapper{M};
        return Wrapper;

    } else if (mass_flag == 0)  {

        thread_local MatrixND<12,12> M{0.0};
        thread_local Matrix Wrapper{M};
        // lumped mass matrix
        double m = 0.5*total_mass;
        M(0,0) = m;
        M(1,1) = m;
        M(2,2) = m;
        M(6,6) = m;
        M(7,7) = m;
        M(8,8) = m;
        return Wrapper;

    } else {
        // consistent (cubic) mass matrix

        // get initial element length
        double L  = this->getLength(State::Init);
        double m  = total_mass/420.0;
        double mx = twist_mass;
        thread_local MatrixND<12,12> ml{0.0};

        ml(0,0) = ml(6,6) = m*140.0;
        ml(0,6) = ml(6,0) = m*70.0;

        ml(3,3) = ml(9,9) = mx/3.0; // Twisting
        ml(3,9) = ml(9,3) = mx/6.0;

        ml( 2, 2) = ml( 8, 8) =  m*156.0;
        ml( 2, 8) = ml( 8, 2) =  m*54.0;
        ml( 4, 4) = ml(10,10) =  m*4.0*L*L;
        ml( 4,10) = ml(10, 4) = -m*3.0*L*L;
        ml( 2, 4) = ml( 4, 2) = -m*22.0*L;
        ml( 8,10) = ml(10, 8) = -ml(2,4);
        ml( 2,10) = ml(10, 2) =  m*13.0*L;
        ml( 4, 8) = ml( 8, 4) = -ml(2,10);

        ml( 1, 1) = ml( 7, 7) =  m*156.0;
        ml( 1, 7) = ml( 7, 1) =  m*54.0;
        ml( 5, 5) = ml(11,11) =  m*4.0*L*L;
        ml( 5,11) = ml(11, 5) = -m*3.0*L*L;
        ml( 1, 5) = ml( 5, 1) =  m*22.0*L;
        ml( 7,11) = ml(11, 7) = -ml(1,5);
        ml( 1,11) = ml(11, 1) = -m*13.0*L;
        ml( 5, 7) = ml( 7, 5) = -ml(1,11);

        // transform local mass matrix to global system
        return theCoordTransf->getGlobalMatrixFromLocal(ml);
    }
}

/*
const Matrix &
ForceFrame3d::getInitialStiff()
{
  // check for quick return
  if (Ki != 0)
    return *Ki;

  static Matrix f(nq,nq);   // element flexibility matrix  
  this->getInitialFlexibility(f);
    
  // calculate element stiffness matrix
  static Matrix kvInit(nq, nq);
  if (f.Invert(kvInit) < 0)
    opserr << "ForceFrame3d::getInitialStiff -- could not invert flexibility";

  Ki = new Matrix(theCoordTransf->getInitialGlobalStiffMatrix(kvInit));

  return *Ki;
}
*/



void
ForceFrame3d::initializeSectionHistoryVariables()
{
  for (int i = 0; i < numSections; i++) {

    fs[i]      = MatrixND<nsr, nsr>{};
    es[i]      = VectorND<nsr>{};
    Ssr[i]     = VectorND<nsr>{};
    es_save[i] = VectorND<nsr>{};
  }
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS *********************/
int
ForceFrame3d::update()
{
  this->BasicFrame3d::update();

  // if have completed a recvSelf() - do a revertToLastCommit
  // to get Ssr, etc. set correctly
  if (initialFlag == 2)
    this->revertToLastCommit();

  // update the transformation
  theCoordTransf->update();

  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  // get basic displacements and increments
  const Vector& v = theCoordTransf->getBasicTrialDisp();

  thread_local VectorND<nq> dv;
  dv = theCoordTransf->getBasicIncrDeltaDisp();

  if (initialFlag != 0 && (dv.norm() <= DBL_EPSILON) && (eleLoads.size()==0))
    return 0;

  thread_local VectorND<nq> Dv;
  Dv  = v;
  Dv -= dv;


  thread_local VectorND<nq> vr;        // element residual displacements
  thread_local MatrixND<nq, nq> F;   // element flexibility matrix
  thread_local VectorND<nq> dvToDo, dvTrial;

  dvToDo  = dv;
  dvTrial = dvToDo;


  int numSubdivide = 1;
  bool converged   = false;

  const double factor = 10.0;
  double dW;             // section strain energy (work) norm
  double dW0  = 0.0;


  //
  // fmk - modification to get compatible ele forces and deformations
  //   for a change in deformation dV we try first a Newton iteration, if
  //   that fails we try an initial flexibility iteration on first iteration
  //   and then regular Newton, if that fails we use the initial flexiblity
  //   for all iterations.
  //
  //   if they both fail we subdivide dV & try to get compatible forces
  //   and deformations. if they work and we have subdivided we apply
  //   the remaining dV.
  //
  while (converged == false && numSubdivide <= maxSubdivisions) {

    // try regular Newton (if l==0), or
    // initial tangent iterations (if l==1), or
    // initial tangent on first iteration then regular Newton (if l==2)

    for (int l = 0; l < 3; l++) {

      for (int i = 0; i < numSections; i++) {
        es_trial[i]  = es[i];
        Fs_trial[i]  = fs[i];
        sr_trial[i] = Ssr[i];
      }

      if (initialFlag == 2)
        continue;

      VectorND<nq>       q_trial = q_pres;
      MatrixND<nq, nq> K_trial = kv;

      // calculate nodal force increments and update nodal forces
      VectorND<nq> dqe = kv * dv;
      q_trial += dqe;

      // Allow 10 times more iterations for initial tangent strategy
      int numIters = (l==1) ? 10*maxIters : maxIters;

      for (int j = 0; j < numIters; j++) {
        F.zero();
        vr.zero();

        Matrix f(F);
        if (beamIntegr->addElasticFlexibility(L, f) < 0) {
          vr[0] += F(0, 0) * q_trial[0];
          vr[1] += F(1, 1) * q_trial[1] + F(1, 2) * q_trial[2];
          vr[2] += F(2, 1) * q_trial[1] + F(2, 2) * q_trial[2];
          vr[3] += F(3, 3) * q_trial[3] + F(3, 4) * q_trial[4];
          vr[4] += F(4, 3) * q_trial[3] + F(4, 4) * q_trial[4];
          vr[5] += F(5, 5) * q_trial[5];
        }


        // Add effects of element loads
        double v0[5];
        v0[0] = v0[1] = v0[2] = v0[3] = v0[4] = 0.0;
        for (auto[load, factor] : eleLoads)
          beamIntegr->addElasticDeformations(load, factor, L, v0);

        vr[0] += v0[0];
        vr[1] += v0[1];
        vr[2] += v0[2];
        vr[3] += v0[3];
        vr[4] += v0[4];

        //
        // Gauss Loop
        //
        for (int i = 0; i < numSections; i++) {
            double xL  = xi[i];
            double xL1 = xL - 1.0;
            double wtL = wt[i] * L;

            auto& Fs = Fs_trial[i];
            auto& sr = sr_trial[i];

            //
            // a. Calculate interpolated section force
            //
            //    si = b*q + bp*w;

            // Interpolation of q_trial
            //    b*q_trial
            VectorND<6> si {
                 q_trial[0],                           // N
                 oneOverL * (q_trial[1] + q_trial[2]), // VY
                 oneOverL * (q_trial[3] + q_trial[4]), // VZ
                 q_trial[5],                           // T
                 xL1 * q_trial[3] + xL * q_trial[4],   // MY
                 xL1 * q_trial[1] + xL * q_trial[2],   // MZ
            };
            // Add the effects of element loads
            // si += bp*w
            if (eleLoads.size() != 0)
              this->addLoadAtSection(si, i);


            //
            // b. Compute section strain es_trial
            //
            //    es += Fs * ( si - sr(e) );
            //
            if (initialFlag != 0) {

              VectorND<nsr> ds; //, des;

              // Form stress increment ds
              // ds = si - sr(e);
              ds = si;
              ds.addVector(1.0, sr, -1.0);

              // Add strain increment
              //    es += Fs * ds;
              if (l == 0) {
                  //  regular Newton
                  es_trial[i].addMatrixVector(1.0, Fs, ds, 1.0);

              } else if (l == 2) {
                  //  Newton with initial tangent if first iteration
                  //  otherwise regular Newton
                  if (j == 0) {
                    MatrixND<nsr,nsr> Fs0 = sections[i]->getFlexibility<nsr,scheme>(State::Init);
                    es_trial[i].addMatrixVector(1.0, Fs0, ds, 1.0);
                  } else
                    es_trial[i].addMatrixVector(1.0, Fs, ds, 1.0);

              } else {
                  //  Newton with initial tangent
                  MatrixND<nsr,nsr> Fs0 = sections[i]->getFlexibility<nsr,scheme>(State::Init);
                  es_trial[i].addMatrixVector(1.0, Fs0, ds, 1.0);
              }
            }


            //
            // c. Set trial section state and get response
            //
            if (sections[i]->setTrialState<nsr,scheme>(es_trial[i]) < 0) {
              opserr << "ForceFrame3d::update() - section failed in setTrial\n";
              return -1;
            }

            sr_trial[i] = sections[i]->getResultant<nsr, scheme>();
            Fs_trial[i] = sections[i]->getFlexibility<nsr, scheme>();

            //
            // d. Integrate element flexibility matrix
            //
            //    F += (B' * Fs * B) * wi * L;
            {
              MatrixND<nsr,nq> FsB;
              FsB.zero();
              for (int jj = 0; jj < nsr; jj++) {
                  // SECTION_RESPONSE_P:
                  FsB(jj, 0) += Fs_trial[i](jj, 0) * wtL;

                  // SECTION_RESPONSE_VY:
                  FsB(jj, 1) += Fs(jj, 1) * wtL * oneOverL;
                  FsB(jj, 2) += Fs(jj, 1) * wtL * oneOverL;

                  // SECTION_RESPONSE_VZ:
                  FsB(jj, 3) += Fs(jj, 2) * wtL * oneOverL;
                  FsB(jj, 4) += Fs(jj, 2) * wtL * oneOverL;

                  // SECTION_RESPONSE_T:
                  FsB(jj, 5) += Fs(jj, 3) * wtL;

                  // SECTION_RESPONSE_MY:
                  FsB(jj, 3) += xL1 * Fs(jj, 4)*wtL;
                  FsB(jj, 4) += xL  * Fs(jj, 4)*wtL;

                  // SECTION_RESPONSE_MZ:
                  FsB(jj, 1) += xL1 * Fs(jj, 5)*wtL;
                  FsB(jj, 2) += xL  * Fs(jj, 5)*wtL;
              }

              for (int jj = 0; jj < nq; jj++) {
                  double tmp;
                  // SECTION_RESPONSE_P:
                  F(0, jj) += FsB( 0, jj);

                  // SECTION_RESPONSE_VY:
                  tmp = oneOverL * FsB( 1, jj);
                  F(1, jj) += tmp;
                  F(2, jj) += tmp;

                  // SECTION_RESPONSE_VZ:
                  tmp = oneOverL * FsB( 2, jj);
                  F(3, jj) += tmp;
                  F(4, jj) += tmp;

                  // SECTION_RESPONSE_T:
                  F(5, jj) += FsB( 3, jj);

                  // SECTION_RESPONSE_MY:
                  F(3, jj) += xL1 * FsB( 4, jj);
                  F(4, jj) += xL  * FsB( 4, jj);

                  // SECTION_RESPONSE_MZ:
                  F(1, jj) += xL1 * FsB( 5, jj);
                  F(2, jj) += xL  * FsB( 5, jj);
              }
            }


            //
            // e. Integrate residual deformations
            //
            //    vr += (B' * (es + des)) * wi * L;
            {
              VectorND<nsr> des, ds;
              // calculate section residual deformations
              // des = Fs * ds,  with  ds = si - sr[i];
              ds = si;
              ds.addVector(1.0, sr, -1.0);

              des.addMatrixVector(0.0, Fs, ds, 1.0);
              des.addVector(1.0, es_trial[i], 1.0);

              // SECTION_RESPONSE_P:
              vr[0] += des[0]*wtL;
              // SECTION_RESPONSE_VY:
              vr[1] += oneOverL*des[1]*wtL;
              vr[2] += oneOverL*des[1]*wtL;
              // SECTION_RESPONSE_VZ:
              vr[3] += oneOverL*des[2]*wtL;
              vr[4] += oneOverL*des[2]*wtL;
              // SECTION_RESPONSE_T:
              vr[5] += des[3]*wtL;
              // SECTION_RESPONSE_MY:
              vr[3] += xL1 * des[4]*wtL;
              vr[4] += xL * des[4]*wtL;
              // SECTION_RESPONSE_MZ:
              vr[1] += xL1 * des[5]*wtL;
              vr[2] += xL * des[5]*wtL;

            }

        } // Gauss loop


        //
        // Finalize trial element state 
        //
        //    K_trial  = int(F)
        //    q_trial += K * (Dv + dv_trial - vr)
        //
        if (F.invert(K_trial) < 0)
          opserr << "ForceFrame3d::update -- could not invert flexibility\n";


        // dv = Dv + dvTrial  - vr
        dv = Dv;
        dv += dvTrial;
        dv -= vr;

        // dqe = kv * dv;
        dqe.addMatrixVector(0.0, K_trial, dv, 1.0);

        dW = dqe.dot(dv);
        if (dW0 == 0.0)
          dW0 = dW;

        q_trial += dqe;


        //
        // check for convergence of this interval
        //
        if (fabs(dW) < tol) {

          // set the target displacement
          dvToDo -= dvTrial;
          Dv  += dvTrial;

          // check if we have got to where we wanted
          if (dvToDo.norm() <= DBL_EPSILON) {
            converged = true;

          } else { // we convreged but we have more to do
            // reset variables for start of next subdivision
            dvTrial      = dvToDo;
            // NOTE setting subdivide to 1 again maybe too much
            numSubdivide = 1; 
          }

          // set kv, es and q_pres values
          kv = K_trial;
          q_pres = q_trial;

          for (int k = 0; k < numSections; k++) {
            es[k]  = es_trial[k];
            fs[k]  = Fs_trial[k];
            Ssr[k] = sr_trial[k];
          }

          // break out of j & l loops
          j = numIters + 1;
          l = 4;

        }
        else { //  if (fabs(dW) < tol) {

          // if we have failed to converge for all of our Newton schemes
          // - reduce step size by the factor specified
          if (j == (numIters - 1) && (l == 2)) {
            dvTrial /= factor;
            numSubdivide++;
          }

        }
      } // for (int j=0; j < numIters; j++)
    }   // for (int l=0; l < 3; l++)
  }     // while (converged == false)

  // if fail to converge we return an error flag & print an error message

  if (converged == false) {
    opserr << "WARNING - ForceFrame3d::update - failed to get compatible ";
    opserr << "element forces & deformations for element: ";
    opserr << this->getTag() << "; dW: << " << dW << ", dW0: " << dW0 << ")\n";

    return -1;
  }

  initialFlag = 1;

  return 0;
}


void
ForceFrame3d::addLoadAtSection(VectorND<nsr>& sp, int isec)
{

  double L = theCoordTransf->getInitialLength();

  double x = xi[isec] * L;

  for (auto[load, loadFactor] : eleLoads) {

    int type;
    const Vector& data = load->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wa = data(2) * loadFactor; // Axial
      double wy = data(0) * loadFactor; // Transverse
      double wz = data(1) * loadFactor; // Transverse

      for (int ii = 0; ii < nsr; ii++) {
        switch (scheme[ii]) {
        case SECTION_RESPONSE_P:  sp[ii] += wa * (L - x); break;
        case SECTION_RESPONSE_VY: sp[ii] += wy * (x - 0.5 * L); break;
        case SECTION_RESPONSE_VZ: sp[ii] += wz * (0.5 * L - x); break;
        case SECTION_RESPONSE_MZ: sp[ii] += wy * 0.5 * x * (x - L); break;
        case SECTION_RESPONSE_MY: sp[ii] += wz * 0.5 * x * (L - x); break;
        default:                  break;
        }
      }

    } else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
      double wa = data(2) * loadFactor; // Axial
      double wy = data(0) * loadFactor; // Transverse
      double wz = data(1) * loadFactor; // Transverse
      double a  = data(3) * L;
      double b  = data(4) * L;

      double Fa  = wa * (b - a); // resultant axial load
      double Fy  = wy * (b - a); // resultant transverse load
      double Fz  = wz * (b - a); // resultant transverse load
      double c   = a + 0.5 * (b - a);
      double VyI = Fy * (1 - c / L);
      double VyJ = Fy * c / L;
      double VzI = Fz * (1 - c / L);
      double VzJ = Fz * c / L;

      for (int ii = 0; ii < nsr; ii++) {

        if (x <= a) {
          switch (scheme[ii]) {
          case SECTION_RESPONSE_P:  sp[ii] += Fa; break;
          case SECTION_RESPONSE_VY: sp[ii] -= VyI; break;
          case SECTION_RESPONSE_VZ: sp[ii] += VzI; break;
          case SECTION_RESPONSE_MZ: sp[ii] -= VyI * x; break;
          case SECTION_RESPONSE_MY: sp[ii] += VzI * x; break;
          }
        } else if (x >= b) {
          switch (scheme[ii]) {
          case SECTION_RESPONSE_VY: sp[ii] += VyJ; break;
          case SECTION_RESPONSE_VZ: sp[ii] -= VzJ; break;
          case SECTION_RESPONSE_MZ: sp[ii] += VyJ * (x - L); break;
          case SECTION_RESPONSE_MY: sp[ii] -= VzJ * (x - L); break;
          }
        } else {
          switch (scheme[ii]) {
          case SECTION_RESPONSE_P:  sp[ii] +=  Fa  - wa * (x - a); break;
          case SECTION_RESPONSE_VY: sp[ii] += -VyI + wy * (x - a); break;
          case SECTION_RESPONSE_VZ: sp[ii] -= -VzI + wz * (x - a); break;
          case SECTION_RESPONSE_MY: sp[ii] +=  VzI * x - 0.5 * wz * x * x - wz * a * (0.5 * a - x);  break;
          case SECTION_RESPONSE_MZ: sp[ii] += -VyI * x + 0.5 * wy * x * x + wy * a * (0.5 * a - x); break;
          default:                  break;
          }
        }
      }

    } else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py     = data(0) * loadFactor;
      double Pz     = data(1) * loadFactor;
      double N      = data(2) * loadFactor;
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      double a = aOverL * L;

      double Vy1 = Py * (1.0 - aOverL);
      double Vy2 = Py * aOverL;

      double Vz1 = Pz * (1.0 - aOverL);
      double Vz2 = Pz * aOverL;

      for (int ii = 0; ii < nsr; ii++) {

        if (x <= a) {
          switch (scheme[ii]) {
          case SECTION_RESPONSE_P:  sp[ii] += N; break;
          case SECTION_RESPONSE_VY: sp[ii] -= Vy1; break;
          case SECTION_RESPONSE_VZ: sp[ii] -= Vz1; break;
          case SECTION_RESPONSE_MY: sp[ii] += x * Vz1; break;
          case SECTION_RESPONSE_MZ: sp[ii] -= x * Vy1; break;
          default:                  break;
          }
        } else {
          switch (scheme[ii]) {
          case SECTION_RESPONSE_MZ: sp[ii] -= (L - x) * Vy2; break;
          case SECTION_RESPONSE_VY: sp[ii] += Vy2; break;
          case SECTION_RESPONSE_MY: sp[ii] += (L - x) * Vz2; break;
          case SECTION_RESPONSE_VZ: sp[ii] += Vz2; break;
          default:                  break;
          }
        }
      }
    } else {
      opserr << "ForceFrame3d::addLoad -- load type unknown for element with tag: "
             << this->getTag() << "\n";
    }
  }
}

void
ForceFrame3d::computeSectionForceSensitivity(Vector& dspdh, int isec, int gradNumber)
{
  int type;

  double L    = theCoordTransf->getInitialLength();
  double dLdh = theCoordTransf->getdLdh();

  double dxidh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dxidh);

  double x    = xi[isec] * L;
  double dxdh = xi[isec] * dLdh + dxidh[isec] * L;

  int order      = sections[isec]->getOrder();
  const ID& code = sections[isec]->getType();

  for (auto[load, loadFactor] : eleLoads) {
    const  Vector& data = load->getData(type, loadFactor);

//  const Vector& data = eleLoads[i]->getData(type, 1.0);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0) * 1.0; // Transverse
      double wz = data(1) * 1.0; // Transverse
      double wa = data(2) * 1.0; // Axial

      const Vector& sens = load->getSensitivityData(gradNumber);
      double dwydh       = sens(0);
      double dwzdh       = sens(1);
      double dwadh       = sens(2);

      for (int ii = 0; ii < order; ii++) {

        switch (scheme[ii]) {
        case SECTION_RESPONSE_P:
          //sp(ii) += wa*(L-x);
          dspdh(ii) += dwadh * (L - x) + wa * (dLdh - dxdh);
          break;
        case SECTION_RESPONSE_MZ:
          //sp(ii) += wy*0.5*x*(x-L);
          //dspdh(ii) += 0.5*(dwydh*x*(x-L) + wy*dxdh*(x-L) + wy*x*(dxdh-dLdh));
          dspdh(ii) += 0.5 * (dwydh * x * (x - L) + wy * (dxdh * (2 * x - L) - x * dLdh));
          break;
        case SECTION_RESPONSE_VY:
          //sp(ii) += wy*(x-0.5*L);
          dspdh(ii) += dwydh * (x - 0.5 * L) + wy * (dxdh - 0.5 * dLdh);
          break;
        case SECTION_RESPONSE_MY:
          //sp(ii) += wz*0.5*x*(L-x);
          //dspdh(ii) += 0.5*(dwzdh*x*(L-x) + wz*dxdh*(L-x) + wz*x*(dLdh-dxdh));
          dspdh(ii) += 0.5 * (dwzdh * x * (L - x) + wz * (dxdh * (L - 2 * x) + x * dLdh));
          break;
        case SECTION_RESPONSE_VZ:
          //sp(ii) += wz*(x-0.5*L);
          dspdh(ii) += dwzdh * (0.5 * L - x) + wz * (0.5 * dLdh - dxdh);
          break;
        default: break;
        }
      }
    } else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py     = data(0) * 1.0;
      double Pz     = data(1) * 1.0;
      double N      = data(2) * 1.0;
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      const Vector& sens = load->getSensitivityData(gradNumber);
      double dPydh       = sens[0];
      double dPzdh       = sens[1];
      double dNdh        = sens[2];
      double daLdh       = sens[3];

      double a = aOverL * L;

      double Vy1    = Py * (1.0 - aOverL);
      double Vy2    = Py * aOverL;
      double dVy1dh = Py * (0.0 - daLdh) + dPydh * (1.0 - aOverL);
      double dVy2dh = Py * daLdh + dPydh * aOverL;

      double Vz1    = Pz * (1.0 - aOverL);
      double Vz2    = Pz * aOverL;
      double dVz1dh = Pz * (0.0 - daLdh) + dPzdh * (1.0 - aOverL);
      double dVz2dh = Pz * daLdh + dPzdh * aOverL;

      for (int ii = 0; ii < order; ii++) {

        if (x <= a) {
          switch (scheme[ii]) {
          case SECTION_RESPONSE_P:
            //sp(ii) += N;
            dspdh(ii) += dNdh;
            break;
          case SECTION_RESPONSE_MZ:
            //sp(ii) -= x*Vy1;
            dspdh(ii) -= (dxdh * Vy1 + x * dVy1dh);
            break;
          case SECTION_RESPONSE_VY:
            //sp(ii) -= Vy1;
            dspdh(ii) -= dVy1dh;
            break;
          case SECTION_RESPONSE_MY:
            //sp(ii) += x*Vz1;
            dspdh(ii) += (dxdh * Vz1 + x * dVz1dh);
            break;
          case SECTION_RESPONSE_VZ:
            //sp(ii) -= Vz1;
            dspdh(ii) -= dVz1dh;
            break;
          default: break;
          }
        } else {
          switch (scheme[ii]) {
          case SECTION_RESPONSE_MZ:
            //sp(ii) -= (L-x)*Vy2;
            dspdh(ii) -= (dLdh - dxdh) * Vy2 + (L - x) * dVy2dh;
            break;
          case SECTION_RESPONSE_VY:
            //sp(ii) += Vy2;
            dspdh(ii) += dVy2dh;
            break;
          case SECTION_RESPONSE_MY:
            //sp(ii) += (L-x)*Vz2;
            dspdh(ii) += (dLdh - dxdh) * Vz2 + (L - x) * dVz2dh;
            break;
          case SECTION_RESPONSE_VZ:
            //sp(ii) += Vz2;
            dspdh(ii) += dVz2dh;
            break;
          default: break;
          }
        }
      }
    } else {
      opserr << "ForceFrame3d::computeSectionForceSensitivity -- load type unknown for element "
                "with tag: "
             << this->getTag() << "\n";
    }
  }
}

int
ForceFrame3d::sendSelf(int commitTag, Channel& theChannel)
{
  // place the integer data into an ID
  int dbTag = this->getDbTag();
  int i, j, k;
  int loc = 0;

  ID idData(11);
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = numSections;
  idData(4) = maxIters;
  idData(5) = initialFlag;

  idData(7)          = theCoordTransf->getClassTag();
  int crdTransfDbTag = theCoordTransf->getDbTag();
  if (crdTransfDbTag == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag != 0)
      theCoordTransf->setDbTag(crdTransfDbTag);
  }
  idData(8) = crdTransfDbTag;


  idData(9)           = beamIntegr->getClassTag();
  int beamIntegrDbTag = beamIntegr->getDbTag();
  if (beamIntegrDbTag == 0) {
    beamIntegrDbTag = theChannel.getDbTag();
    if (beamIntegrDbTag != 0)
      beamIntegr->setDbTag(beamIntegrDbTag);
  }
  idData(10) = beamIntegrDbTag;

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "ForceFrame3d::sendSelf() - failed to send ID data\n";
    return -1;
  }

  // send the coordinate transformation

  if (theCoordTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ForceFrame3d::sendSelf() - failed to send crdTranf\n";
    return -1;
  }

  if (beamIntegr->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ForceFrame3d::sendSelf() - failed to send beamIntegr\n";
    return -1;
  }

  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2 * numSections);
  loc = 0;
  for (i = 0; i < numSections; i++) {
    int sectClassTag = sections[i]->getClassTag();
    int sectDbTag    = sections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      sections[i]->setDbTag(sectDbTag);
    }

    idSections(loc)     = sectClassTag;
    idSections(loc + 1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0) {
    opserr << "ForceFrame3d::sendSelf() - failed to send ID data\n";
    return -1;
  }

  //
  // send the sections
  //

  for (j = 0; j < numSections; j++) {
    if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "ForceFrame3d::sendSelf() - section " << j << "failed to send itself\n";
      return -1;
    }
  }

  // into a vector place distrLoadCommit, rho, UeCommit, q_save and K_save
  int secDefSize = 0;
  for (i = 0; i < numSections; i++) {
    int size = sections[i]->getOrder();
    secDefSize += size;
  }

  Vector dData(1 + 1 + 1 + nq + nq * nq + secDefSize + 4);
  loc = 0;

  // place double variables into Vector
  dData(loc++) = density;
  dData(loc++) = tol;

  // put  distrLoadCommit into the Vector
  //  for (i=0; i<NL; i++)
  //dData(loc++) = distrLoadcommit(i);

  // place K_save into vector
  for (i = 0; i < nq; i++)
    dData(loc++) = q_save[i];

  // place K_save into vector
  for (i = 0; i < nq; i++)
    for (j = 0; j < nq; j++)
      dData(loc++) = K_save(i, j);

  // place es_save into vector
  for (k = 0; k < numSections; k++)
    for (i = 0; i < sections[k]->getOrder(); i++)
      dData(loc++) = es_save[k][i];

  // send damping coefficients
  dData(loc++) = alphaM;
  dData(loc++) = betaK;
  dData(loc++) = betaK0;
  dData(loc++) = betaKc;

  if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
    opserr << "ForceFrame3d::sendSelf() - failed to send Vector data\n";

    return -1;
  }
  return 0;
}

int
ForceFrame3d::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i, j, k;

  ID idData(11); // one bigger than needed

  if (theChannel.recvID(dbTag, commitTag, idData) < 0) {
    opserr << "ForceFrame3d::recvSelf() - failed to recv ID data\n";

    return -1;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  maxIters                  = idData(4);
  initialFlag               = idData(5);

  int crdTransfClassTag = idData(7);
  int crdTransfDbTag    = idData(8);

  int beamIntegrClassTag = idData(9);
  int beamIntegrDbTag    = idData(10);

  // create a new crdTransf object if one needed
  if (theCoordTransf == 0 || theCoordTransf->getClassTag() != crdTransfClassTag) {
    if (theCoordTransf != 0)
      delete theCoordTransf;

    // TODO(cmp) - add FrameTransform to ObjBroker
    theCoordTransf = nullptr; // theBroker.getNewFrameTransform3d(crdTransfClassTag);

    if (theCoordTransf == nullptr) {
      opserr << "ForceFrame3d::recvSelf() - failed to obtain a CrdTrans object with classTag"
             << crdTransfClassTag << "\n";
      return -2;
    }
  }

  theCoordTransf->setDbTag(crdTransfDbTag);
  // invoke recvSelf on the crdTransf obkject
  if (theCoordTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "ForceFrame3d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }

  // create a new beamIntegr object if one needed
  if (beamIntegr == 0 || beamIntegr->getClassTag() != beamIntegrClassTag) {
    if (beamIntegr != 0)
      delete beamIntegr;

    beamIntegr = theBroker.getNewBeamIntegration(beamIntegrClassTag);

    if (beamIntegr == 0) {
      opserr
          << "ForceFrame3d::recvSelf() - failed to obtain the beam integration object with classTag"
          << beamIntegrClassTag << "\n";
      return -1;
    }
  }

  beamIntegr->setDbTag(beamIntegrDbTag);

  // invoke recvSelf on the beamIntegr object
  if (beamIntegr->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "ForceFrame3d::sendSelf() - failed to recv beam integration\n";

    return -3;
  }

  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2 * idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0) {
    opserr << "ForceFrame3d::recvSelf() - failed to recv ID data\n";
    return -1;
  }

  //
  // now receive the sections
  //
  if (numSections != idData(3)) {

    //
    // we do not have correct number of sections, must delete the old and create
    // new ones before can recvSelf on the sections
    //

    // delete the old
    if (numSections != 0) {
      for (int i = 0; i < numSections; i++)
        delete sections[i];
      delete[] sections;
    }

    // create a section and recvSelf on it
    numSections = idData(3);

    // Delete the old
    if (es_save != 0)
      delete[] es_save;

    // Allocate
    es_save = new VectorND<nsr>[numSections]{};

    // Delete the old
    if (fs != 0)
      delete[] fs;

    fs = new MatrixND<nsr,nsr>[numSections]{};

    // Delete the old
    if (es != 0)
      delete[] es;

    // Allocate the right number
    es = new VectorND<nsr>[numSections]{};

    // Delete the old
    if (Ssr != 0)
      delete[] Ssr;

    // Allocate
    Ssr = new VectorND<nsr>[numSections]{};


    // create a new array to hold pointers
    sections = new FrameSection*[idData(3)];

    loc = 0;

    for (i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag    = idSections(loc + 1);
      loc += 2;
      // TODO(cmp) add FrameSection to broker
//    sections[i] = theBroker.getNewSection(sectClassTag);
      if (sections[i] == 0) {
        opserr << "ForceFrame3d::recvSelf() - Broker could not create Section of class type"
               << sectClassTag << "\n";
        return -1;
      }

      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "ForceFrame3d::recvSelf() - section " << i << " failed to recv itself\n";
        return -1;
      }
    }

    this->initializeSectionHistoryVariables();

  } else {

    //
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //

    loc = 0;
    for (i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag    = idSections(loc + 1);
      loc += 2;

      // check of correct type
      if (sections[i]->getClassTag() != sectClassTag) {
        // delete the old section[i] and create a new one
        delete sections[i];
      // TODO(cmp) add FrameSection to broker
//      sections[i] = theBroker.getNewSection(sectClassTag);
        if (sections[i] == 0) {
          opserr << "ForceFrame3d::recvSelf() - Broker could not create Section of class type "
                 << sectClassTag << "\n";
          ;
          return -1;
        }
      }

      // recvvSelf on it
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "ForceFrame3d::recvSelf() - section " << i << " failed to recv itself\n";
        return -1;
      }
    }
  }

  // into a vector place distrLoadCommit, rho, UeCommit, q_save and K_save
  int secDefSize = 0;
  for (int ii = 0; ii < numSections; ii++) {
    int size = sections[ii]->getOrder();
    secDefSize += size;
  }

  Vector dData(1 + 1 + nq + nq * nq + secDefSize + 4);

  if (theChannel.recvVector(dbTag, commitTag, dData) < 0) {
    opserr << "ForceFrame3d::recvSelf() - failed to send Vector data\n";
    return -1;
  }

  loc = 0;

  // place double variables into Vector
  density = dData(loc++);
  tol = dData(loc++);

  // put  distrLoadCommit into the Vector
  //for (i=0; i<NL; i++)
  // distrLoad(i) = dData(loc++);

  // place q_save into vector
  for (i = 0; i < nq; i++)
    q_save[i] = dData(loc++);

  // place K_save into matrix
  for (i = 0; i < nq; i++)
    for (j = 0; j < nq; j++)
      K_save(i, j) = dData(loc++);

  kv = K_save;
  q_pres = q_save;

  for (k = 0; k < numSections; k++) {
    int order = sections[k]->getOrder();

    // place es_save into vector
    es_save[k] = VectorND<nsr>{};
    for (i = 0; i < order; i++)
      es_save[k][i] = dData(loc++);
  }

  // set damping coefficients
  alphaM = dData(loc++);
  betaK  = dData(loc++);
  betaK0 = dData(loc++);
  betaKc = dData(loc++);

  initialFlag = 2;
  return 0;
}

// addBFsB(F, s, i)
int
ForceFrame3d::getInitialFlexibility(MatrixND<nq,nq>& fe)
{
  fe.zero();

  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  // Flexibility from elastic interior
  {
    Matrix wrapper(fe);
    beamIntegr->addElasticFlexibility(L, wrapper);
  }

  for (int i = 0; i < numSections; i++) {


    double xL  = xi[i];
    double xL1 = xL - 1.0;
    double wtL = wt[i] * L;

    const MatrixND<nsr,nsr> fSec = sections[i]->getFlexibility<nsr, scheme>(State::Init);

    thread_local MatrixND<nsr, nq> FB;
    FB.zero();
    double tmp;
    int ii, jj;
    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:
        for (int jj = 0; jj < nsr; jj++)
          FB(jj, 0) += fSec(jj, ii) * wtL;
        break;
      case SECTION_RESPONSE_MZ:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = fSec(jj, ii) * wtL;
          FB(jj, 1) += xL1 * tmp;
          FB(jj, 2) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = oneOverL * fSec(jj, ii) * wtL;
          FB(jj, 1) += tmp;
          FB(jj, 2) += tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = fSec(jj, ii) * wtL;
          FB(jj, 3) += xL1 * tmp;
          FB(jj, 4) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (int jj = 0; jj < nsr; jj++) {
          tmp = oneOverL * fSec(jj, ii) * wtL;
          FB(jj, 3) += tmp;
          FB(jj, 4) += tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (int jj = 0; jj < nsr; jj++)
          FB(jj, 5) += fSec(jj, ii) * wtL;
        break;
      default: break;
      }
    }
    //
    for (int ii = 0; ii < nsr; ii++) {
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:
        for (int jj = 0; jj < nq; jj++)
          fe(0, jj) += FB(ii, jj);
        break;
      case SECTION_RESPONSE_MZ:
        for (int jj = 0; jj < nq; jj++) {
          tmp = FB(ii, jj);
          fe(1, jj) += xL1 * tmp;
          fe(2, jj) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int jj = 0; jj < nq; jj++) {
          tmp = oneOverL * FB(ii, jj);
          fe(1, jj) += tmp;
          fe(2, jj) += tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (int jj = 0; jj < nq; jj++) {
          tmp = FB(ii, jj);
          fe(3, jj) += xL1 * tmp;
          fe(4, jj) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (int jj = 0; jj < nq; jj++) {
          tmp = oneOverL * FB(ii, jj);
          fe(3, jj) += tmp;
          fe(4, jj) += tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (int jj = 0; jj < nq; jj++)
          fe(5, jj) += FB(ii, jj);
        break;
      default: break;
      }
    }
  }
  return 0;
}

int
ForceFrame3d::getInitialDeformations(Vector& v0)
{
  v0.Zero();
  if (eleLoads.size() < 1 || (this->setState(State::Init) != 0))
    return 0;

  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;


  for (int i = 0; i < numSections; i++) {

    double xL  = xi[i];
    double xL1 = xL - 1.0;
    double wtL = wt[i] * L;

    thread_local VectorND<nsr> sp;
    sp.zero();

    this->addLoadAtSection(sp, i);

    MatrixND<nsr,nsr> fse = sections[i]->getFlexibility<nsr,scheme>(State::Init);

    VectorND<nsr> e;

    e.addMatrixVector(0.0, fse, sp, 1.0);

    double dei, tmp;
    for (int ii = 0; ii < nsr; ii++) {
      dei = e[ii] * wtL;
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P: v0(0) += dei; break;
      case SECTION_RESPONSE_MZ:
        v0[1] += xL1 * dei;
        v0[2] += xL * dei;
        break;
      case SECTION_RESPONSE_VY:
        tmp = oneOverL * dei;
        v0[1] += tmp;
        v0[2] += tmp;
        break;
      case SECTION_RESPONSE_MY:
        v0[3] += xL1 * dei;
        v0[4] += xL * dei;
        break;
      case SECTION_RESPONSE_VZ:
        tmp = oneOverL * dei;
        v0[3] += tmp;
        v0[4] += tmp;
        break;
      default: break;
      }
    }
  }

  return 0;
}

void
ForceFrame3d::compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const
{
  // get basic displacements and increments
  static Vector ub(nq);
  ub = theCoordTransf->getBasicTrialDisp();

  double L = theCoordTransf->getInitialLength();

  // setup Vandermode and CBDI influence matrices

  // get CBDI influence matrix
  Matrix ls(numSections, numSections);
  getCBDIinfluenceMatrix(numSections, xi, L, ls);

  // get section curvatures
  Vector kappa_y(numSections);
  Vector kappa_z(numSections);

  for (int i = 0; i < numSections; i++) {
    // get section deformations
    VectorND<nsr> es = sections[i]->getDeformation<nsr,scheme>();

    for (int j = 0; j < nsr; j++) {
      if (scheme[j] == SECTION_RESPONSE_MZ)
        kappa_z(i) = es[j];
      if (scheme[j] == SECTION_RESPONSE_MY)
        kappa_y(i) = es[j];
    }
  }

  Vector v(numSections), w(numSections);
  static VectorND<ndm> xl, uxb;
  static Vector xg(ndm), uxg(ndm);
  // double theta;                             // angle of twist of the sections

  // v = ls * kappa_z;
  v.addMatrixVector(0.0, ls, kappa_z, 1.0);
  // w = ls * kappa_y *  (-1);
  w.addMatrixVector(0.0, ls, kappa_y, -1.0);

  for (int i = 0; i < numSections; i++) {

    xl(0) = xi[i] * L;
    xl(1) = 0;
    xl(2) = 0;

    // get section global coordinates
    sectionCoords[i] = theCoordTransf->getPointGlobalCoordFromLocal(xl);

    // compute section displacements
    //theta  = xi * ub(5); // consider linear variation for angle of twist. CHANGE LATER!!!!!!!!!!
    uxb(0) = xi[i] * ub(0); // consider linear variation for axial displacement. CHANGE LATER!!!!!!!!!!
    uxb(1) = v[i];
    uxb(2) = w[i];

    // get section displacements in global system
    sectionDispls[i] = theCoordTransf->getPointGlobalDisplFromBasic(xi[i], uxb);
  }
  return;
}

void
ForceFrame3d::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"ForceFrame3d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " 
                        << connectedExternalNodes(1) << "], ";

    // Mass
    double mass;
    if (getIntegral(Field::Density, State::Init, mass) == 0)
      s << ", \"mass\": " << mass << ", ";
    else
      s << ", \"massperlength\": " << density << ", ";

    s << "\"sections\": [";
    for (int i = 0; i < numSections - 1; i++)
      s << sections[i]->getTag() << ", ";
    s << sections[numSections - 1]->getTag() << "], ";
    s << "\"crdTransformation\": " << theCoordTransf->getTag()  << ", ";
    s << "\"integration\": ";
    beamIntegr->Print(s, flag);
    s << "}";
  }

  // flags with negative values are used by GSA
  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s << sections[0]->getTag() << "\t" << sections[numSections - 1]->getTag();
    s << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
    s << "\t0\t0.0000000\n";
  }

  // flags with negative values are used by GSA
  else if (flag < -1) {
    int eleTag   = this->getTag();
    int counter  = (flag + 1) * -1;
    double P     = q_save[0];
    double MZ1   = q_save[1];
    double MZ2   = q_save[2];
    double MY1   = q_save[3];
    double MY2   = q_save[4];
    double L     = theCoordTransf->getInitialLength();
    double VY    = (MZ1 + MZ2) / L;
    theVector(1) = VY;
    theVector(4) = -VY;
    double VZ    = (MY1 + MY2) / L;
    double T     = q_save[5];

    double p0[5];
    p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
    if (eleLoads.size() > 0)
      this->computeReactions(p0);

    s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
    s << "\t" << -P + p0[0] << "\t" << VY + p0[1] << "\t" << -VZ + p0[3] << "\n";
    s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
    s << "\t" << P << ' ' << -VY + p0[2] << ' ' << VZ + p0[4] << "\n";
    s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
    s << "\t" << -T << "\t" << MY1 << "\t" << MZ1 << "\n";
    s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
    s << "\t" << T << ' ' << MY2 << ' ' << MZ2 << "\n";
  }

  // flag set to 2 used to print everything .. used for viewing data for UCSD renderer
  else if (flag == 2) {
    static Vector xAxis(3);
    static Vector yAxis(3);
    static Vector zAxis(3);

    theCoordTransf->getLocalAxes(xAxis, yAxis, zAxis);

    s << "#ForceFrame3D\n";
    s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2);
    s << " " << yAxis(0) << " " << yAxis(1) << " " << yAxis(2) << " ";
    s << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << "\n";

    const Vector& node1Crd  = theNodes[0]->getCrds();
    const Vector& node2Crd  = theNodes[1]->getCrds();
    const Vector& node1Disp = theNodes[0]->getDisp();
    const Vector& node2Disp = theNodes[1]->getDisp();

    s << "#NODE " << node1Crd(0) << " " << node1Crd(1) << " " << node1Crd(2) << " " << node1Disp(0)
      << " " << node1Disp(1) << " " << node1Disp(2) << " " << node1Disp(3) << " " << node1Disp(4)
      << " " << node1Disp(5) << "\n";

    s << "#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2) << " " << node2Disp(0)
      << " " << node2Disp(1) << " " << node2Disp(2) << " " << node2Disp(3) << " " << node2Disp(4)
      << " " << node2Disp(5) << "\n";

    double P     = q_save[0];
    double MZ1   = q_save[1];
    double MZ2   = q_save[2];
    double MY1   = q_save[3];
    double MY2   = q_save[4];
    double L     = theCoordTransf->getInitialLength();
    double VY    = (MZ1 + MZ2) / L;
    theVector(1) = VY;
    theVector(4) = -VY;
    double VZ    = (MY1 + MY2) / L;
    double T     = q_save[5];

    double p0[5];
    p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
    if (eleLoads.size() > 0)
      this->computeReactions(p0);

    s << "#END_FORCES " << -P + p0[0] << ' ' << VY + p0[1] << ' ' << -VZ + p0[3] << ' ' << -T << ' '
      << MY1 << ' ' << MZ1 << "\n";
    s << "#END_FORCES " << P << ' ' << -VY + p0[2] << ' ' << VZ + p0[4] << ' ' << T << ' ' << MY2
      << ' ' << MZ2 << "\n";

    // plastic hinge rotation
    static Vector vp(6);
    static MatrixND<nq,nq> fe;
    this->getInitialFlexibility(fe);
    vp = theCoordTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, q_pres, -1.0);
    s << "#PLASTIC_HINGE_ROTATION " << vp[1] << " " << vp[2] << " " << vp[3] << " " << vp[4] << " "
      << 0.1 * L << " " << 0.1 * L << "\n";

    // allocate array of vectors to store section coordinates and displacements
    static int maxNumSections = 0;
    static Vector* coords     = 0;
    static Vector* displs     = 0;
    if (maxNumSections < numSections) {
      if (coords != 0)
        delete[] coords;
      if (displs != 0)
        delete[] displs;

      coords = new Vector[numSections];
      displs = new Vector[numSections];

      for (int i = 0; i < numSections; i++)
        coords[i] = Vector(ndm);

      for (int i = 0; i < numSections; i++)
        displs[i] = Vector(ndm);

      maxNumSections = numSections;
    }

    // compute section location & displacements
    this->compSectionDisplacements(coords, displs);

    // print the section location & invoke print on the section
    for (int i = 0; i < numSections; i++) {
      s << "#SECTION " << (coords[i])(0) << " " << (coords[i])(1) << " " << (coords[i])(2);
      s << " " << (displs[i])(0) << " " << (displs[i])(1) << " " << (displs[i])(2) << "\n";
      sections[i]->Print(s, flag);
    }
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nElement: " << this->getTag() << " Type: ForceFrame3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << density << "\n";
    beamIntegr->Print(s, flag);
    double P     = q_save[0];
    double MZ1   = q_save[1];
    double MZ2   = q_save[2];
    double MY1   = q_save[3];
    double MY2   = q_save[4];
    double L     = theCoordTransf->getInitialLength();
    double VY    = (MZ1 + MZ2) / L;
    theVector[1] = VY;
    theVector[4] = -VY;
    double VZ    = (MY1 + MY2) / L;
    double T     = q_save[5];

    double p0[5];
    p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
    if (eleLoads.size() > 0)
      this->computeReactions(p0);

    s << "\tEnd 1 Forces (P MZ VY MY VZ T): " << -P + p0[0] << " " << MZ1 << " " << VY + p0[1]
      << " " << MY1 << " " << -VZ + p0[3] << " " << T << "\n";
    s << "\tEnd 2 Forces (P MZ VY MY VZ T): " << P << " " << MZ2 << " " << -VY + p0[2] << " " << MY2
      << " " << VZ + p0[4] << " " << -T << "\n";

    for (int i = 0; i < numSections; i++)
      sections[i]->Print(s, flag);
  }
}

OPS_Stream&
operator<<(OPS_Stream& s, ForceFrame3d& E)
{
  E.Print(s);
  return s;
}


Response*
ForceFrame3d::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  Response* theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "ForceFrame3d");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  //
  // we compare argv[0] for known response types
  //

  // global force -
  if (strcmp(argv[0], "forces") == 0 || strcmp(argv[0], "force") == 0 ||
      strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Pz_1");
    output.tag("ResponseType", "Mx_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Pz_2");
    output.tag("ResponseType", "Mx_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");


    theResponse = new ElementResponse(this, 1, theVector);

    // local force -
  } else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "Vy_1");
    output.tag("ResponseType", "Vz_1");
    output.tag("ResponseType", "T_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "Vy_2");
    output.tag("ResponseType", "Vz_2");
    output.tag("ResponseType", "T_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, 2, theVector);

    // basic force -
  } else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {

    output.tag("ResponseType", "N");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Mz_2");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "T");

    theResponse = new ElementResponse(this, 7, Vector(6));

    // basic stiffness -
  } else if (strcmp(argv[0], "basicStiffness") == 0) {

    output.tag("ResponseType", "N");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Mz_2");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "T");

    theResponse = new ElementResponse(this, 19, Matrix(6, 6));

    // chord rotation -
  } else if (strcmp(argv[0], "chordRotation") == 0 || strcmp(argv[0], "chordDeformation") == 0 ||
             strcmp(argv[0], "basicDeformation") == 0) {

    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "thetaZ_1");
    output.tag("ResponseType", "thetaZ_2");
    output.tag("ResponseType", "thetaY_1");
    output.tag("ResponseType", "thetaY_2");
    output.tag("ResponseType", "thetaX");

    theResponse = new ElementResponse(this, 3, Vector(6));

    // plastic rotation -
  } else if (strcmp(argv[0], "plasticRotation") == 0 ||
             strcmp(argv[0], "plasticDeformation") == 0) {

    output.tag("ResponseType", "epsP");
    output.tag("ResponseType", "thetaZP_1");
    output.tag("ResponseType", "thetaZP_2");
    output.tag("ResponseType", "thetaYP_1");
    output.tag("ResponseType", "thetaYP_2");
    output.tag("ResponseType", "thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(6));

    // point of inflection
  } else if (strcmp(argv[0], "inflectionPoint") == 0) {
    theResponse = new ElementResponse(this, 5, Vector(2));

    // tangent drift
  } else if (strcmp(argv[0], "tangentDrift") == 0) {
    theResponse = new ElementResponse(this, 6, Vector(4));

  } else if (strcmp(argv[0], "getRemCriteria1") == 0) {
    theResponse = new ElementResponse(this, 77, Vector(2));

  } else if (strcmp(argv[0], "getRemCriteria2") == 0) {
    theResponse = new ElementResponse(this, 8, Vector(2), ID(6));

  } else if (strcmp(argv[0], "RayleighForces") == 0 || strcmp(argv[0], "rayleighForces") == 0) {

    theResponse = new ElementResponse(this, 12, theVector);

  } else if (strcmp(argv[0], "sections") == 0) {
    CompositeResponse* theCResponse = new CompositeResponse();
    int numResponse                 = 0;
    double xi[maxNumSections];
    double L = theCoordTransf->getInitialLength();
    beamIntegr->getSectionLocations(numSections, L, xi);

    for (int i = 0; i < numSections; i++) {

      output.tag("GaussPointOutput");
      output.attr("number", i + 1);
      output.attr("eta", xi[i] * L);

      Response* theSectionResponse = sections[i]->setResponse(&argv[1], argc - 1, output);

      if (theSectionResponse != 0) {
        numResponse = theCResponse->addResponse(theSectionResponse);
      }
    }

    if (numResponse == 0) // no valid responses found
      delete theCResponse;
    else
      theResponse = theCResponse;
  }

  else if (strcmp(argv[0], "integrationPoints") == 0)
    theResponse = new ElementResponse(this, 10, Vector(numSections));

  else if (strcmp(argv[0], "integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(numSections));

  else if (strcmp(argv[0], "sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(numSections));

  else if (strcmp(argv[0], "sectionDisplacements") == 0) {
    if (argc > 1 && strcmp(argv[1], "local") == 0)
      theResponse = new ElementResponse(this, 1111, Matrix(numSections, 3));
    else
      theResponse = new ElementResponse(this, 111, Matrix(numSections, 3));
  }

  else if (strcmp(argv[0], "cbdiDisplacements") == 0)
    theResponse = new ElementResponse(this, 112, Matrix(20, 3));


  else if (strstr(argv[0], "section") != 0) {

    if (argc > 1) {

      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {
        double xi[maxNumSections];
        double L = theCoordTransf->getInitialLength();
        beamIntegr->getSectionLocations(numSections, L, xi);

        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", 2.0 * xi[sectionNum - 1] - 1.0);

        if (strcmp(argv[2], "dsdh") != 0) {
          theResponse = sections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);
        } else {
          int order         = sections[sectionNum - 1]->getOrder();
          theResponse       = new ElementResponse(this, 76, Vector(order));
          Information& info = theResponse->getInformation();
          info.theInt       = sectionNum;
        }

        output.endTag();

      } else if (sectionNum == 0) { // argv[1] was not an int, we want all sections,

        CompositeResponse* theCResponse = new CompositeResponse();
        int numResponse                 = 0;
        double xi[maxNumSections];
        double L = theCoordTransf->getInitialLength();
        beamIntegr->getSectionLocations(numSections, L, xi);

        for (int i = 0; i < numSections; i++) {

          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", xi[i] * L);

          Response* theSectionResponse = sections[i]->setResponse(&argv[1], argc - 1, output);

          if (theSectionResponse != 0) {
            numResponse = theCResponse->addResponse(theSectionResponse);
          }
        }

        if (numResponse == 0) // no valid responses found
          delete theCResponse;
        else
          theResponse = theCResponse;
      }
    }
  }
  //by SAJalali
  else if (strcmp(argv[0], "energy") == 0) {
    return new ElementResponse(this, 2000, 0.0);
  }

  if (theResponse == nullptr) {
    theResponse = theCoordTransf->setResponse(argv, argc, output);
  }

  output.endTag();

  return theResponse;
}

int
ForceFrame3d::getResponse(int responseID, Information& eleInfo)
{
  static Vector vp(6);
  static MatrixND<nq,nq> fe;

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2) {
    double p0[5];
    p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
    if (eleLoads.size() > 0)
      this->computeReactions(p0);
    // Axial
    double N     = q_pres[0];
    theVector(6) =  N;
    theVector(0) = -N + p0[0];

    // Torsion
    double T     = q_pres[5];
    theVector(9) =  T;
    theVector(3) = -T;

    // Moments about z and shears along y
    double M1     = q_pres[1];
    double M2     = q_pres[2];
    theVector(5)  = M1;
    theVector(11) = M2;
    double L      = theCoordTransf->getInitialLength();
    double V      = (M1 + M2) / L;
    theVector(1)  =  V + p0[1];
    theVector(7)  = -V + p0[2];

    // Moments about y and shears along z
    M1            = q_pres[3];
    M2            = q_pres[4];
    theVector(4)  = M1;
    theVector(10) = M2;
    V             = (M1 + M2) / L;
    theVector(2)  = -V + p0[3];
    theVector(8)  =  V + p0[4];

    return eleInfo.setVector(theVector);

  }

  // Chord rotation
  else if (responseID == 3) {
    vp = theCoordTransf->getBasicTrialDisp();
    return eleInfo.setVector(vp);
  }

  else if (responseID == 7)
    return eleInfo.setVector(q_pres);

  else if (responseID == 19)
    return eleInfo.setMatrix(kv);

  // Plastic rotation
  else if (responseID == 4) {
    this->getInitialFlexibility(fe);
    vp = theCoordTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, q_pres, -1.0);
    Vector v0(6);
    this->getInitialDeformations(v0);
    vp.addVector(1.0, v0, -1.0);
    return eleInfo.setVector(vp);
  }

  else if (responseID == 10) {
    double L = theCoordTransf->getInitialLength();
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i] * L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 11) {
    double L = theCoordTransf->getInitialLength();
    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i] * L;
    return eleInfo.setVector(weights);
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = sections[i]->getTag();
    return eleInfo.setID(tags);
  }

  else if (responseID == 111 || responseID == 1111) {
    double L = theCoordTransf->getInitialLength();
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    // CBDI influence matrix
    Matrix ls(numSections, numSections);
    getCBDIinfluenceMatrix(numSections, pts, L, ls);
    // Curvature vector
    Vector kappaz(numSections); // about section z
    Vector kappay(numSections); // about section y
    for (int i = 0; i < numSections; i++) {
      const ID& code  = sections[i]->getType();
      const Vector& e = sections[i]->getSectionDeformation();
      int order       = sections[i]->getOrder();
      for (int j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_MZ)
          kappaz(i) += e(j);
        if (code(j) == SECTION_RESPONSE_MY)
          kappay(i) += e(j);
      }
    }
    // Displacement vector
    Vector dispsy(numSections); // along local y
    Vector dispsz(numSections); // along local z
    dispsy.addMatrixVector(0.0, ls, kappaz, 1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);
    beamIntegr->getSectionLocations(numSections, L, pts);
    static Vector uxb(3);
    static Vector uxg(3);
    Matrix disps(numSections, 3);
    vp = theCoordTransf->getBasicTrialDisp();
    for (int i = 0; i < numSections; i++) {
      uxb(0) = pts[i] * vp(0); // linear shape function
      uxb(1) = dispsy(i);
      uxb(2) = dispsz(i);
      if (responseID == 111)
        uxg = theCoordTransf->getPointGlobalDisplFromBasic(pts[i], uxb);
      else
        uxg = theCoordTransf->getPointLocalDisplFromBasic(pts[i], uxb);
      disps(i, 0) = uxg(0);
      disps(i, 1) = uxg(1);
      disps(i, 2) = uxg(2);
    }
    return eleInfo.setMatrix(disps);
  }

  else if (responseID == 112) {
    double L = theCoordTransf->getInitialLength();
    double ipts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, ipts);
    // CBDI influence matrix
    double pts[20];
    for (int i = 0; i < 20; i++)
      pts[i] = 1.0 / (20 - 1) * i;
    Matrix ls(20, numSections);
    getCBDIinfluenceMatrix(20, pts, numSections, ipts, L, ls);
    // Curvature vector
    Vector kappaz(numSections); // about section z
    Vector kappay(numSections); // about section y
    for (int i = 0; i < numSections; i++) {
      const ID& code  = sections[i]->getType();
      const Vector& e = sections[i]->getSectionDeformation();
      int order       = sections[i]->getOrder();
      for (int j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_MZ)
          kappaz(i) += e(j);
        if (code(j) == SECTION_RESPONSE_MY)
          kappay(i) += e(j);
      }
    }
    // Displacement vector
    Vector dispsy(20); // along local y
    Vector dispsz(20); // along local z
    dispsy.addMatrixVector(0.0, ls, kappaz, 1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);
    static Vector uxb(3);
    static Vector uxg(3);
    Matrix disps(20, 3);
    vp = theCoordTransf->getBasicTrialDisp();
    for (int i = 0; i < 20; i++) {
      uxb(0)      = pts[i] * vp(0); // linear shape function
      uxb(1)      = dispsy(i);
      uxb(2)      = dispsz(i);
      uxg         = theCoordTransf->getPointGlobalDisplFromBasic(pts[i], uxb);
      disps(i, 0) = uxg(0);
      disps(i, 1) = uxg(1);
      disps(i, 2) = uxg(2);
    }
    return eleInfo.setMatrix(disps);
  }

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());

  // Point of inflection
  else if (responseID == 5) {
    static Vector LI(2);
    LI(0) = 0.0;
    LI(1) = 0.0;

    double L = theCoordTransf->getInitialLength();

    if (fabs(q_pres[1] + q_pres[2]) > DBL_EPSILON)
      LI(0) = q_pres[1] / (q_pres[1] + q_pres[2]) * L;

    if (fabs(q_pres[3] + q_pres[4]) > DBL_EPSILON)
      LI(1) = q_pres[3] / (q_pres[3] + q_pres[4]) * L;

    return eleInfo.setVector(LI);
  }

  // Tangent drift
  else if (responseID == 6) {
    double d2z = 0.0;
    double d2y = 0.0;
    double d3z = 0.0;
    double d3y = 0.0;

    double L = theCoordTransf->getInitialLength();

    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);

    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);

    // Location of inflection point from node I
    double LIz = 0.0;
    if (fabs(q_pres[1] + q_pres[2]) > DBL_EPSILON)
      LIz = q_pres[1] / (q_pres[1] + q_pres[2]) * L;

    double LIy = 0.0;
    if (fabs(q_pres[3] + q_pres[4]) > DBL_EPSILON)
      LIy = q_pres[3] / (q_pres[3] + q_pres[4]) * L;

    for (int i = 0; i < numSections; i++) {
      double x       = pts[i] * L;
      const ID& type = sections[i]->getType();
      int order      = sections[i]->getOrder();
      double kappa   = 0.0;
      if (x < LIz) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MZ)
            kappa += es[i][j];
        double b = -LIz + x;
        d2z += (wts[i] * L) * kappa * b;
      }
      kappa = 0.0;
      if (x < LIy) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MY)
            kappa += es[i][j];
        double b = -LIy + x;
        d2y += (wts[i] * L) * kappa * b;
      }
    }

    d2z += beamIntegr->getTangentDriftI(L, LIz, q_pres[1], q_pres[2]);
    d2y += beamIntegr->getTangentDriftI(L, LIy, q_pres[3], q_pres[4], true);

    for (int i = numSections - 1; i >= 0; i--) {
      double x       = pts[i] * L;
      const ID& type = sections[i]->getType();
      int order      = sections[i]->getOrder();
      double kappa   = 0.0;
      if (x > LIz) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MZ)
            kappa += es[i][j];
        double b = x - LIz;
        d3z += (wts[i] * L) * kappa * b;
      }
      kappa = 0.0;
      if (x > LIy) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MY)
            kappa += es[i][j];
        double b = x - LIy;
        d3y += (wts[i] * L) * kappa * b;
      }
    }

    d3z += beamIntegr->getTangentDriftJ(L, LIz, q_pres[1], q_pres[2]);
    d3y += beamIntegr->getTangentDriftJ(L, LIy, q_pres[3], q_pres[4], true);

    static Vector d(4);
    d(0) = d2z;
    d(1) = d3z;
    d(2) = d2y;
    d(3) = d3y;

    return eleInfo.setVector(d);

  } else if (responseID == 77) { // Why is this here?
    return -1;

  } else if (responseID == 8) {

    ID* eleInfoID = eleInfo.theID;

    int compID     = (*eleInfoID)(0);
    int critID     = (*eleInfoID)(1);
    int nTagbotn11 = (*eleInfoID)(2);
    int nTagmidn11 = (*eleInfoID)(3);
    int nTagtopn11 = (*eleInfoID)(4);
    int globgrav11 = (*eleInfoID)(5);

    const char* filenamewall = eleInfo.theString;

    // int returns
    double value       = 0.0;
    double checkvalue1 = 0.0;

    if (critID == 7) {
      Domain* theDomain = this->getDomain();

      double oofwallresp;
      // determine the in plane horizontal deformation axis
      // and the out of plane horizontal deformation axis
      Node* theNode1a      = theDomain->getNode(nTagbotn11);
      Node* theNode3a      = theDomain->getNode(nTagtopn11);
      const Vector& crdIa1 = theNode1a->getCrds();
      const Vector& crdJa1 = theNode3a->getCrds();
      int indwdir1;
      int indwdir2;
      if (globgrav11 == 1) {
        indwdir1 = 1;
        indwdir2 = 2;
      } else if (globgrav11 == 2) {
        indwdir1 = 0;
        indwdir2 = 2;
      } else if (globgrav11 == 3) {
        indwdir1 = 0;
        indwdir2 = 1;
      }

      double dir1a1    = crdJa1(indwdir1) - crdIa1(indwdir1);
      double dir2a1    = crdJa1(indwdir2) - crdIa1(indwdir2);
      double dirsumdum = sqrt(dir1a1 * dir1a1 + dir2a1 * dir2a1);
      double dir1inp   = dir1a1 / dirsumdum;
      double dir2inp   = dir2a1 / dirsumdum;

      double dir1oop = -dir2inp;
      double dir2oop = dir1inp;

      Node* theNode1                = theDomain->getNode(nTagbotn11);
      const Vector& theResponsewall = theNode1->getTrialDisp();
      double valbotinfn = theResponsewall(indwdir1) * dir1inp + theResponsewall(indwdir2) * dir2inp;
      double valbotoutfn =
          theResponsewall(indwdir1) * dir1oop + theResponsewall(indwdir2) * dir2oop;

      Node* theNode2                 = theDomain->getNode(nTagmidn11);
      const Vector& theResponsewall2 = theNode2->getTrialDisp();
      double valmidinfn =
          theResponsewall2(indwdir1) * dir1inp + theResponsewall2(indwdir2) * dir2inp;
      double valmidoutfn =
          theResponsewall2(indwdir1) * dir1oop + theResponsewall2(indwdir2) * dir2oop;

      Node* theNode3                 = theDomain->getNode(nTagtopn11);
      const Vector& theResponsewall3 = theNode3->getTrialDisp();
      double valtopinfn =
          theResponsewall3(indwdir1) * dir1inp + theResponsewall3(indwdir2) * dir2inp;
      double valtopoutfn =
          theResponsewall3(indwdir1) * dir1oop + theResponsewall3(indwdir2) * dir2oop;

      value             = sqrt(pow((valtopinfn - valbotinfn), 2.0));
      double valoutchck = valmidoutfn - (valtopoutfn + valbotoutfn) / 2.0;
      oofwallresp       = sqrt(pow(valoutchck, 2.0));
      //
      double outplanevaldat;  // variable for input value
      double inplanevaldat;   // variable for input value
      double outplanevaldat1; // variable for input value
      double inplanevaldat1;  // variable for input value
      std::ifstream indata;

      if (filenamewall != nullptr) {
        //
        indata.open(filenamewall); // opens the file
        if (!indata) {             // file couldn't be opened
          opserr << "ForceFrame3d::getResponse"
                 << " file for infill wall (" << filenamewall << " could not be opened" << "\n";
          return -1;
        }
        checkvalue1    = 0.0;
        int counterdum = 0;
        while (!indata.eof()) { // keep reading until end-of-file
          counterdum = counterdum + 1;
          indata >> outplanevaldat >> inplanevaldat; // sets EOF flag if no value found
          if (counterdum != 1) {
            if (oofwallresp >= outplanevaldat1 && oofwallresp <= outplanevaldat) {
              checkvalue1 = inplanevaldat1 + (oofwallresp - outplanevaldat1) /
                                                 (outplanevaldat - outplanevaldat1) *
                                                 (inplanevaldat - inplanevaldat1);
              break;
            }
          }
          indata >> outplanevaldat1 >> inplanevaldat1;
          if (oofwallresp >= outplanevaldat && oofwallresp <= outplanevaldat1) {
            checkvalue1 = inplanevaldat + (oofwallresp - outplanevaldat) /
                                              (outplanevaldat1 - outplanevaldat) *
                                              (inplanevaldat1 - inplanevaldat);
            break;
          }
        }
        indata.close();
      }

      static Vector result8(2);
      result8(0) = value;
      result8(1) = checkvalue1;

      return eleInfo.setVector(result8);
    }

    return -1;
  }
  //by SAJalali
  else if (responseID == 2000) {
    double xi[maxNumSections];
    double L = theCoordTransf->getInitialLength();
    beamIntegr->getSectionWeights(numSections, L, xi);
    double energy = 0;
    for (int i = 0; i < numSections; i++) {
      energy += sections[i]->getEnergy() * xi[i] * L;
    }
    return eleInfo.setDouble(energy);
  }

  return -1;
}

int
ForceFrame3d::getResponseSensitivity(int responseID, int gradNumber, Information& eleInfo)
{
  // Basic deformation sensitivity
  if (responseID == 3) {
    const Vector& dvdh = theCoordTransf->getBasicDisplSensitivity(gradNumber);
    return eleInfo.setVector(dvdh);
  }

  // Basic force sensitivity
  else if (responseID == 7) {
    static Vector dqdh(6);

    const Vector& dvdh = theCoordTransf->getBasicDisplSensitivity(gradNumber);

    dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);

    dqdh.addVector(1.0, this->computedqdh(gradNumber), 1.0);

    return eleInfo.setVector(dqdh);
  }

  // dsdh
  else if (responseID == 76) {

    int sectionNum = eleInfo.theInt;
    int order      = sections[sectionNum - 1]->getOrder();

    Vector dsdh(order);
    dsdh.Zero();

    if (eleLoads.size() > 0) {
      this->computeSectionForceSensitivity(dsdh, sectionNum - 1, gradNumber);
    }
    //opserr << "FBC3d::getRespSens dspdh: " << dsdh;
    static Vector dqdh(6);

    const Vector& dvdh = theCoordTransf->getBasicDisplSensitivity(gradNumber);

    dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);

    dqdh.addVector(1.0, this->computedqdh(gradNumber), 1.0);

    //opserr << "FBC3d::getRespSens dqdh: " << dqdh;

    double L        = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);

    const ID& code = sections[sectionNum - 1]->getType();

    double xL  = pts[sectionNum - 1];
    double xL1 = xL - 1.0;

    for (int ii = 0; ii < order; ii++) {
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:  dsdh(ii) += dqdh(0); break;
      case SECTION_RESPONSE_MZ: dsdh(ii) += xL1 * dqdh(1) + xL * dqdh(2); break;
      case SECTION_RESPONSE_VY: dsdh(ii) += oneOverL * (dqdh(1) + dqdh(2)); break;
      case SECTION_RESPONSE_MY: dsdh(ii) += xL1 * dqdh(3) + xL * dqdh(4); break;
      case SECTION_RESPONSE_VZ: dsdh(ii) += oneOverL * (dqdh(3) + dqdh(4)); break;
      case SECTION_RESPONSE_T:  dsdh(ii) += dqdh(5); break;
      default:                  dsdh(ii) += 0.0; break;
      }
    }

    double dLdh   = theCoordTransf->getdLdh();
    double d1oLdh = theCoordTransf->getd1overLdh();

    double dptsdh[maxNumSections];
    beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);
    double dxLdh = dptsdh[sectionNum - 1]; // - xL/L*dLdh;

    for (int j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
        dsdh(j) += dxLdh * (q_pres[1] + q_pres[2]);
        //dsdh(j) -= dLdh*xL/L*(Se(1)+Se(2));
        break;
      case SECTION_RESPONSE_VY: dsdh(j) += d1oLdh * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_MY: dsdh(j) += dxLdh  * (q_pres[3] + q_pres[4]); break;
      case SECTION_RESPONSE_VZ: dsdh(j) += d1oLdh * (q_pres[3] + q_pres[4]); break;
      default:                  break;
      }
    }

    return eleInfo.setVector(dsdh);
  }

  // Plastic deformation sensitivity
  else if (responseID == 4) {
    static Vector dvpdh(6);

    const Vector& dvdh = theCoordTransf->getBasicDisplSensitivity(gradNumber);

    dvpdh = dvdh;

    static MatrixND<nq,nq> fe;
    this->getInitialFlexibility(fe);

    const Vector& dqdh = this->computedqdh(gradNumber);

    dvpdh.addMatrixVector(1.0, fe, dqdh, -1.0);

    dvpdh.addMatrixVector(1.0, fe*kv, dvdh, -1.0);

    const Matrix& dfedh = this->computedfedh(gradNumber);

    dvpdh.addMatrixVector(1.0, dfedh, q_pres, -1.0);

    return eleInfo.setVector(dvpdh);
  }

  else
    return -1;
}

int
ForceFrame3d::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  // If the parameter belongs to the element itself
  if ((strcmp(argv[0], "rho") == 0) 
     || (strcmp(argv[0], "density") == 0)) {
    param.setValue(density);
    return param.addObject(1, this);
  }

  // section response -
  if (strstr(argv[0], "sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = theCoordTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);

      sectionLoc /= L;

      float minDistance = fabs(xi[0] - sectionLoc);
      int sectionNum    = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i] - sectionLoc) < minDistance) {
          minDistance = fabs(xi[i] - sectionLoc);
          sectionNum  = i;
        }
      }

      return sections[sectionNum]->setParameter(&argv[2], argc - 2, param);
    }
  }

  // If the parameter belongs to a particular section or lower
  if (strstr(argv[0], "section") != 0) {

    if (argc < 3)
      return -1;

    // Get section number
    int sectionNum = atoi(argv[1]);

    if (sectionNum > 0 && sectionNum <= numSections)
      return sections[sectionNum - 1]->setParameter(&argv[2], argc - 2, param);

    else
      return -1;
  }

  // If the parameter belongs to all sections or lower
  if (strstr(argv[0], "allSections") != 0) {

    if (argc < 2)
      return -1;

    int ok;
    for (int i = 0; i < numSections; i++) {
      ok = sections[i]->setParameter(&argv[1], argc - 1, param);
      if (ok != -1)
        result = ok;
    }

    return result;
  }

  if (strstr(argv[0], "integration") != 0) {

    if (argc < 2)
      return -1;

    return beamIntegr->setParameter(&argv[1], argc - 1, param);
  }

  // Default, send to everything
  int ok;

  for (int i = 0; i < numSections; i++) {
    ok = sections[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  ok = beamIntegr->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  return result;
}

int
ForceFrame3d::updateParameter(int parameterID, Information& info)
{
  if (parameterID == 1) {
    this->density = info.theDouble;
    return 0;
  } else
    return -1;
}

int
ForceFrame3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;
}

const Matrix&
ForceFrame3d::getKiSensitivity(int gradNumber)
{
  theMatrix.Zero();
  return theMatrix;
}

const Matrix&
ForceFrame3d::getMassSensitivity(int gradNumber)
{
  theMatrix.Zero();

  double L = theCoordTransf->getInitialLength();
  if (density != 0.0 && parameterID == 1)
    theMatrix(0, 0) = theMatrix(1, 1) = theMatrix(2, 2) = theMatrix(6, 6) = theMatrix(7, 7) =
        theMatrix(8, 8)                                                   = 0.5 * L;

  return theMatrix;
}

const Vector&
ForceFrame3d::getResistingForceSensitivity(int gradNumber)
{
  static Vector dqdh(6);
  dqdh = this->computedqdh(gradNumber);

  // Transform forces
  double dp0dh[6];
  dp0dh[0] = 0.0;
  dp0dh[1] = 0.0;
  dp0dh[2] = 0.0;
  dp0dh[3] = 0.0;
  dp0dh[4] = 0.0;
  dp0dh[5] = 0.0;
  this->computeReactionSensitivity(dp0dh, gradNumber);
  Vector dp0dhVec(dp0dh, 6);

  static Vector P(12);
  P.Zero();

  if (theCoordTransf->isShapeSensitivity()) {
    // dAdh^T q
    P = theCoordTransf->getGlobalResistingForceShapeSensitivity(q_pres, dp0dhVec, gradNumber);
    // k dAdh u
    const Vector& dAdh_u = theCoordTransf->getBasicTrialDispShapeSensitivity();
    dqdh.addMatrixVector(1.0, kv, dAdh_u, 1.0);
  }

  // A^T (dqdh + k dAdh u)
  P += theCoordTransf->getGlobalResistingForce(dqdh, dp0dhVec);

  return P;
}

int
ForceFrame3d::commitSensitivity(int gradNumber, int numGrads)
{
  int err = 0;

  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double pts[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, pts);

  double wts[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wts);

  double dLdh = theCoordTransf->getdLdh();

  double dptsdh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double d1oLdh = theCoordTransf->getd1overLdh();

  static Vector dqdh(6);
  dqdh = this->computedqdh(gradNumber);

  // dvdh = A dudh + dAdh u
  const Vector& dvdh = theCoordTransf->getBasicDisplSensitivity(gradNumber);
  dqdh.addMatrixVector(1.0, kv, dvdh, 1.0); // A dudh

  if (theCoordTransf->isShapeSensitivity()) {
    //const Vector &dAdh_u = theCoordTransf->getBasicTrialDispShapeSensitivity(gradNumber);
    //dqdh.addMatrixVector(1.0, kv, dAdh_u, 1.0);  // dAdh u
  }

  // Loop over integration points
  for (int i = 0; i < numSections; i++) {

    int order      = sections[i]->getOrder();
    const ID& code = sections[i]->getType();

    double xL  = pts[i];
    double xL1 = xL - 1.0;

    double dxLdh = dptsdh[i];

    Vector ds(nsr);
    ds.Zero();

    // Add sensitivity wrt element loads
    if (eleLoads.size() > 0) {
      this->computeSectionForceSensitivity(ds, i, gradNumber);
    }

    int j;
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:  ds(j) += dqdh(0); break;
      case SECTION_RESPONSE_VY: ds(j) += oneOverL * (dqdh(1) + dqdh(2)); break;
      case SECTION_RESPONSE_VZ: ds(j) += oneOverL * (dqdh(3) + dqdh(4)); break;
      case SECTION_RESPONSE_T:  ds(j) += dqdh(5); break;
      case SECTION_RESPONSE_MY: ds(j) += xL1 * dqdh(3) + xL * dqdh(4); break;
      case SECTION_RESPONSE_MZ: ds(j) += xL1 * dqdh(1) + xL * dqdh(2); break;
      default:                  ds(j) += 0.0; break;
      }
    }

    const Vector& dsdh = sections[i]->getStressResultantSensitivity(gradNumber, true);
    ds -= dsdh;

    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ: ds(j) += dxLdh  * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_VY: ds(j) += d1oLdh * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_MY: ds(j) += dxLdh  * (q_pres[3] + q_pres[4]); break;
      case SECTION_RESPONSE_VZ: ds(j) += d1oLdh * (q_pres[3] + q_pres[4]); break;
      default:                  break;
      }
    }

    Vector de(order);
    const Matrix& fs = sections[i]->getSectionFlexibility();
    de.addMatrixVector(0.0, fs, ds, 1.0);

    err += sections[i]->commitSensitivity(de, gradNumber, numGrads);
  }

  return err;
}

const Vector&
ForceFrame3d::computedqdh(int gradNumber)
{

  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double wts[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wts);

  double dLdh = theCoordTransf->getdLdh();

  double dptsdh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double dwtsdh[maxNumSections];
  beamIntegr->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

  double d1oLdh = theCoordTransf->getd1overLdh();

  static Vector dvdh(6);
  dvdh.Zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order      = sections[i]->getOrder();
    const ID& code = sections[i]->getType();

    double xL  = xi[i];
    double xL1 = xL - 1.0;
    double wtL = wt[i] * L;

    double dxLdh  = dptsdh[i]; // - xL/L*dLdh;
    double dwtLdh = wt[i] * dLdh + dwtsdh[i] * L;


    // Get section stress resultant gradient
    Vector dsdh(order);
    dsdh = sections[i]->getStressResultantSensitivity(gradNumber, true);

    Vector dspdh(order);
    dspdh.Zero();
    // Add sensitivity wrt element loads
    if (eleLoads.size() > 0)
      this->computeSectionForceSensitivity(dspdh, i, gradNumber);

    dsdh.addVector(1.0, dspdh, -1.0);

    int j;
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ: dsdh(j) -= dxLdh  * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_VY: dsdh(j) -= d1oLdh * (q_pres[1] + q_pres[2]); break;
      case SECTION_RESPONSE_MY: dsdh(j) -= dxLdh  * (q_pres[3] + q_pres[4]); break;
      case SECTION_RESPONSE_VZ: dsdh(j) -= d1oLdh * (q_pres[3] + q_pres[4]); break;
      default:                  break;
      }
    }

    Vector dedh(order);
    const Matrix& fs = sections[i]->getSectionFlexibility();
    dedh.addMatrixVector(0.0, fs, dsdh, 1.0);

    for (j = 0; j < order; j++) {
      double dei = dedh(j) * wtL;
      switch (code(j)) {
      case SECTION_RESPONSE_P: dvdh(0) += dei; break;
      case SECTION_RESPONSE_MZ:
        dvdh(1) += xL1 * dei;
        dvdh(2) += xL * dei;
        break;
      case SECTION_RESPONSE_VY:
        dei = oneOverL * dei;
        dvdh(1) += dei;
        dvdh(2) += dei;
        break;
      case SECTION_RESPONSE_MY:
        dvdh(3) += xL1 * dei;
        dvdh(4) += xL * dei;
        break;
      case SECTION_RESPONSE_VZ:
        dei = oneOverL * dei;
        dvdh(3) += dei;
        dvdh(4) += dei;
        break;
      case SECTION_RESPONSE_T: dvdh(5) += dei; break;
      default:                 break;
      }
    }

    const Vector& e = es[i];
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P: dvdh(0) -= e(j) * dwtLdh; break;
      case SECTION_RESPONSE_MZ:
        dvdh(1) -= xL1 * e(j) * dwtLdh;
        dvdh(2) -= xL * e(j) * dwtLdh;

        dvdh(1) -= dxLdh * e(j) * wtL;
        dvdh(2) -= dxLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_VY:
        dvdh(1) -= oneOverL * e(j) * dwtLdh;
        dvdh(2) -= oneOverL * e(j) * dwtLdh;

        dvdh(1) -= d1oLdh * e(j) * wtL;
        dvdh(2) -= d1oLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_MY:
        dvdh(3) -= xL1 * e(j) * dwtLdh;
        dvdh(4) -= xL * e(j) * dwtLdh;

        dvdh(3) -= dxLdh * e(j) * wtL;
        dvdh(4) -= dxLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_VZ:
        dvdh(3) -= oneOverL * e(j) * dwtLdh;
        dvdh(4) -= oneOverL * e(j) * dwtLdh;

        dvdh(3) -= d1oLdh * e(j) * wtL;
        dvdh(4) -= d1oLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_T: dvdh(5) -= e(j) * dwtLdh; break;
      default:                 break;
      }
    }
  }

  static Matrix dfedh(6, 6);
  dfedh.Zero();

  if (beamIntegr->addElasticFlexDeriv(L, dfedh, dLdh) < 0)
    dvdh.addMatrixVector(1.0, dfedh, q_pres, -1.0);


  static Vector dqdh(6);
  dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);


  return dqdh;
}

const Matrix&
ForceFrame3d::computedfedh(int gradNumber)
{
  static Matrix dfedh(6, 6);

  dfedh.Zero();

  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double dLdh   = theCoordTransf->getdLdh();
  double d1oLdh = theCoordTransf->getd1overLdh();

  beamIntegr->addElasticFlexDeriv(L, dfedh, dLdh);

  double dptsdh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double dwtsdh[maxNumSections];
  beamIntegr->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

  for (int i = 0; i < numSections; i++) {

    int order      = sections[i]->getOrder();
    const ID& code = sections[i]->getType();

    Matrix fb(order, nq);
    Matrix fb2(order, nq);

    double xL  = xi[i];
    double xL1 = xL - 1.0;
    double wtL = wt[i] * L;

    double dxLdh  = dptsdh[i];
    double dwtLdh = wt[i] * dLdh + dwtsdh[i] * L;

    const Matrix& fs    = sections[i]->getInitialFlexibility();
    const Matrix& dfsdh = sections[i]->getInitialFlexibilitySensitivity(gradNumber);
    fb.Zero();
    fb2.Zero();

    double tmp;
    int ii, jj;
    for (ii = 0; ii < order; ii++) {
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:
        for (jj = 0; jj < order; jj++) {
          fb(jj, 0) += dfsdh(jj, ii) * wtL; // 1

          //fb(jj,0) += fs(jj,ii)*dwtLdh; // 3

          //fb2(jj,0) += fs(jj,ii)*wtL; // 4
        }
        break;
      case SECTION_RESPONSE_MZ:
        for (jj = 0; jj < order; jj++) {
          tmp = dfsdh(jj, ii) * wtL; // 1
          fb(jj, 1) += xL1 * tmp;
          fb(jj, 2) += xL * tmp;

          tmp = fs(jj, ii) * wtL; // 2
          //fb(jj,1) += dxLdh*tmp;
          //fb(jj,2) += dxLdh*tmp;

          tmp = fs(jj, ii) * dwtLdh; // 3
          //fb(jj,1) += xL1*tmp;
          //fb(jj,2) += xL*tmp;

          tmp = fs(jj, ii) * wtL; // 4
          //fb2(jj,1) += xL1*tmp;
          //fb2(jj,2) += xL*tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (jj = 0; jj < order; jj++) {
          tmp = oneOverL * dfsdh(jj, ii) * wtL;
          fb(jj, 1) += tmp;
          fb(jj, 2) += tmp;
          // TODO: Need to complete for dLdh != 0
        }
        break;
      default: break;
      }
    }
    for (ii = 0; ii < order; ii++) {
      switch (scheme[ii]) {
      case SECTION_RESPONSE_P:
        for (jj = 0; jj < nq; jj++)
          dfedh(0, jj) += fb(ii, jj);
        break;
      case SECTION_RESPONSE_MZ:
        for (jj = 0; jj < nq; jj++) {
          tmp = fb(ii, jj); // 1,2,3
          dfedh(1, jj) += xL1 * tmp;
          dfedh(2, jj) += xL * tmp;

          tmp = fb2(ii, jj); // 4
          //dfedh(1,jj) += dxLdh*tmp;
          //dfedh(2,jj) += dxLdh*tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (jj = 0; jj < nq; jj++) {
          tmp = oneOverL * fb(ii, jj);
          dfedh(1, jj) += tmp;
          dfedh(2, jj) += tmp;
          // TODO: Need to complete for dLdh != 0
        }
        break;
      default:
        break;
      }
    }
  }

  return dfedh;
}

void
ForceFrame3d::setSectionPointers(int numSec, FrameSection** secPtrs)
{
  // TODO: change to assert and make sure its enforced by consrtructor
  if (numSec > maxNumSections) {
    opserr << "Error: ForceFrame3d::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (secPtrs == nullptr) {
    opserr << "Error: ForceFrame3d::setSectionPointers -- invalid section pointer";
  }

  sections = new FrameSection*[numSections];

  for (int i = 0; i < numSections; i++) {

    if (secPtrs[i] == 0) {
      opserr << "Error: ForceFrame3d::setSectionPointers -- null section pointer " << i << "\n";
    }

    sections[i] = secPtrs[i]->getFrameCopy(scheme);

    // Check sections
    // TODO: return int, dont exit
    int sectionKey1 = -1;
    int sectionKey2 = -1;
    const ID& code = sections[i]->getType();
    for (int j = 0; j < code.Size(); j++) {
      if (code(j) == SECTION_RESPONSE_MZ)
        sectionKey1 = j;
      if (code(j) == SECTION_RESPONSE_MY)
        sectionKey2 = j;
    }
    if (sectionKey1 == 0) {
      opserr << "FATAL ForceFrame3d::compSectionResponse - section does not provide Mz response\n";
      exit(-1);
    }
    if (sectionKey2 == 0) {
      opserr << "FATAL ForceFrame3d::compSectionResponse - section does not provide My response\n";
      exit(-1);
    }
  }


  // allocate section flexibility matrices and section deformation vectors
  fs       = new MatrixND<nsr,nsr>[numSections]{};
  es       = new VectorND<nsr>[numSections]{};
  Ssr      = new VectorND<nsr>[numSections]{};
  es_save  = new VectorND<nsr>[numSections]{};
}

#if 0
const Vector &
ForceFrame3d::getResistingForce()
{
  // Will remove once we clean up the corotational 3d transformation -- MHS
  // theCoordTransf->update();

  double p0[5];
  Vector p0Vec(p0, 5);
  p0Vec.Zero();
  
  if (eleLoads.size() > 0)
    this->computeReactions(p0);
  
  theVector =  theCoordTransf->getGlobalResistingForce(Se, p0Vec);
  
  if (density != 0)
    theVector.addVector(1.0, load, -1.0);
  
  return theVector;
}


void 
ForceFrame3d::zeroLoad()
{
  // This is a semi-hack -- MHS
  numEleLoads = 0;
  
  return;
}

void
ForceFrame3d::getDistrLoadInterpolatMatrix(double xi, Matrix& bp, const ID& code)
{
  bp.Zero();

  double L = theCoordTransf->getInitialLength();
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ: // Moment, Mz, interpolation
      bp(i, 1) = xi * (xi - 1) * L * L / 2;
      break;
    case SECTION_RESPONSE_P: // Axial, P, interpolation
      bp(i, 0) = (1 - xi) * L;
      break;
    case SECTION_RESPONSE_VY: // Shear, Vy, interpolation
      bp(i, 1) = (xi - 0.5) * L;
      break;
    case SECTION_RESPONSE_MY: // Moment, My, interpolation
      bp(i, 2) = xi * (1 - xi) * L * L / 2;
      break;
    case SECTION_RESPONSE_VZ: // Shear, Vz, interpolation
      bp(i, 2) = (0.5 - xi) * L;
      break;
    case SECTION_RESPONSE_T: // Torsion, T, interpolation
      break;
    default: break;
    }
  }
}

void
ForceFrame3d::getForceInterpolatMatrix(double xi, Matrix& b, const ID& code)
{
  b.Zero();

  double L = theCoordTransf->getInitialLength();
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ: // Moment, Mz, interpolation
      b(i, 1) = xi - 1.0;
      b(i, 2) = xi;
      break;
    case SECTION_RESPONSE_P: // Axial, P, interpolation
      b(i, 0) = 1.0;
      break;
    case SECTION_RESPONSE_VY: // Shear, Vy, interpolation
      b(i, 1) = b(i, 2) = 1.0 / L;
      break;
    case SECTION_RESPONSE_MY: // Moment, My, interpolation
      b(i, 3) = xi - 1.0;
      b(i, 4) = xi;
      break;
    case SECTION_RESPONSE_VZ: // Shear, Vz, interpolation
      b(i, 3) = b(i, 4) = 1.0 / L;
      break;
    case SECTION_RESPONSE_T: // Torque, T, interpolation
      b(i, 5) = 1.0;
      break;
    default: break;
    }
  }
}

#endif

