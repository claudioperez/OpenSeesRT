/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Purpose: This file contains the class implementation of CTestFixedNumIter.
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 09/05
// Revision: A
//
#include <CTestFixedNumIter.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <Logging.h>


CTestFixedNumIter::CTestFixedNumIter()
    : ConvergenceTest(CONVERGENCE_TEST_CTestFixedNumIter),
    theSOE(0), maxNumIter(0), currentIter(0), printFlag(0),
    norms(1), nType(2)
{

}


CTestFixedNumIter::CTestFixedNumIter(int maxIter, int printIt, int normType)
    : ConvergenceTest(CONVERGENCE_TEST_CTestFixedNumIter),
    theSOE(0), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
    norms(maxNumIter), nType(normType)
{

}


CTestFixedNumIter::~CTestFixedNumIter()
{

}


ConvergenceTest* CTestFixedNumIter::getCopy(int iterations)
{
    CTestFixedNumIter *theCopy ;
    theCopy = new CTestFixedNumIter(iterations, this->printFlag, this->nType) ;

    theCopy->theSOE = this->theSOE ;

    return theCopy ;
}


int CTestFixedNumIter::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();

    return 0;
}


int CTestFixedNumIter::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the
    // return from start() is checked
    if (theSOE == 0)  {
        opserr << "WARNING: CTestFixedNumIter::test() - no SOE set.\n";
        return -1;
    }

    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0)  {
        opserr << "WARNING: CTestFixedNumIter::test() - start() was never invoked.\n";
        return -2;
    }

    // determine the energy & save value in norms vector
    const Vector &b = theSOE->getB();
    const Vector &x = theSOE->getX();
    double product = x ^ b;
    if (product < 0.0)
        product *= -0.5;
    else
        product *= 0.5;

    if (currentIter <= maxNumIter)
        norms(currentIter-1) = product;

    // print the data if required
    if (printFlag & ConvergenceTest::PrintTest)  {
        opserr << LOG_ITERATE << "Iter: " << pad(currentIter);
        opserr << ", EnergyIncr: " << product;
        opserr << " (Norm deltaX: " << x.pNorm(nType) << ", Norm deltaR: " << b.pNorm(nType) << ")\n";
    }

    if (printFlag & ConvergenceTest::PrintTest02)  {
        opserr << LOG_ITERATE << "Iter: " << pad(currentIter);
        opserr << ", EnergyIncr: " << product;
        opserr << " (Norm deltaX: " << x.pNorm(nType) << ", Norm deltaR: " << b.pNorm(nType) << ")\n";
        opserr << "\tdeltaX: " << x << "\tdeltaR: " << b;
    }

    //
    // check if the algorithm converged
    //

    // if converged - print & return ok
    if (currentIter == maxNumIter)  {
        if (printFlag & ConvergenceTest::PrintTest || printFlag & ConvergenceTest::PrintTest02)
            opserr << endln;

        if (printFlag & ConvergenceTest::PrintSuccess)  {
            opserr << LOG_SUCCESS << "Iter: " << pad(currentIter);
            opserr << " last EnergyIncr: " << product;
            opserr << " (Norm deltaX: " << x.pNorm(nType) << ", Norm deltaR: " << b.pNorm(nType) << ")\n";
        }

        // return the number of times test has been called
        return currentIter;
    }

    // algorithm not yet converged - increment counter and return -1
    else {
        currentIter++;
        return ConvergenceTest::Continue;
    }
}


int CTestFixedNumIter::start(void)
{
    if (theSOE == 0) {
        opserr << "WARNING: CTestFixedNumIter::test() - no SOE returning true\n";
        return -1;
    }

    // set iteration count = 1
    currentIter = 1;
    norms.Zero();
    return 0;
}


int CTestFixedNumIter::getNumTests()
{
    return currentIter;
}


int CTestFixedNumIter::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestFixedNumIter::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestFixedNumIter::getNorms()
{
    return norms;
}


int CTestFixedNumIter::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    Vector x(3);
    x(0) = maxNumIter;
    x(1) = printFlag;
    x(2) = nType;
    res = theChannel.sendVector(this->getDbTag(), cTag, x);
    if (res < 0)
        opserr << "CTestFixedNumIter::sendSelf() - failed to send data\n";

    return res;
}


int CTestFixedNumIter::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    Vector x(3);
    res = theChannel.recvVector(this->getDbTag(), cTag, x);

    if (res < 0) {
        opserr << "CTestFixedNumIter::sendSelf() - failed to send data\n";
        maxNumIter = 25;
        printFlag = 0;
        nType = 2;
    } else  {
        maxNumIter = (int) x(0);
        printFlag = (int) x(1);
        nType = (int) x(2);
        norms.resize(maxNumIter);
    }
    return res;
}
