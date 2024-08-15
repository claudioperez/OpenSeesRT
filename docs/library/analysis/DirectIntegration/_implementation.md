

<div style="font-family:monospace">

```cpp
int DirectIntegrationAnalysis::analyzeStep(double dT)
{
  int result = 0;
  Domain *the_Domain = this->getDomainPtr();

  result = theAnalysisModel->analysisStep(dT) {
    return myDomain->analysisStep(dT) { /* NOTHING?? */ };
  };

  if (result < 0)
    return the_Domain->revertToLastCommit(), -2;
  
  int stamp = the_Domain->hasDomainChanged();
  if (stamp != this?->domainStamp) {
    domainStamp = stamp;        
    if (this->domainChanged() < 0)
      return -1;
  }

  result = theIncrIntegrator->newStep(dT) {

// <div class="function">
      theModel->setDisp(...);
      theModel->setVel(...);
      theModel->setAccel(...);

      double newTime = dT + theModel->getCurrentTime()
      theModel->updateDomain(newTime, dT) {
        int res = 0;
        myDomain->applyLoad(newTime) {
            currentTime = newTime;
            dT = currentTime - committedTime; // ?????????????

            NodeIter &theNodeIter = this->getNodes();
            while ((nodePtr = theNodeIter()) != 0)
              nodePtr->zeroUnbalancedLoad();

            ElementIter &theElemIter = this->getElements();    
            while ((elePtr = theElemIter()) != 0)
              if (elePtr->isSubdomain() == false)
                elePtr->zeroLoad();    

            // now loop over load patterns, invoking applyLoad on them
            LoadPatternIter &thePatterns = this->getLoadPatterns();
            while ((thePattern = thePatterns()) != 0)
              thePattern->(UniformExcitation::)applyLoad(newTime) {
                  while ((theNode = (theDomain->...)theNodes()) != 0) {
                      theNode->setNumColR(1);
                      const Vector &crds = theNode->getCrds();
                      double yCrd = crds(1), zCrd = crds(2);
                      int ndm = crds.Size();
                      ... if (ndm == 3) {
                          if (theDof < 3) {
                              theNode->setR(theDof, 0, fact);
                          } else if (theDof == 3) {
                              theNode->setR(1, 0, -fact*zCrd);
                              theNode->setR(2, 0,  fact*yCrd);
                              theNode->setR(3, 0,  fact);
                          } else if (theDof == 4) {
                              theNode->setR(0, 0,  fact*zCrd);
                              theNode->setR(2, 0, -fact*xCrd);
                              theNode->setR(4, 0,  fact);
                          } else if (theDof == 5) {
                              theNode->setR(0, 0, -fact*yCrd);
                              theNode->setR(1, 0,  fact*xCrd);
                              theNode->setR(5, 0,  fact);
                          }
                      }
                  }

                  this->EarthquakePattern::applyLoad(time) {
                    currentTime = time;
                    for (int i=0; i<numMotions; i++)
                      (*uDotDotG)(i) = theMotions[i]->(GroundMotionRecord::)getAccel(currentTime) {
                        return theAccelTimeSeries->getFactor(currentTime);
                      };

                    while ((theNode = theDomain->theNodes()) != 0) 
                      theNode->addInertiaLoadToUnbalance(*uDotDotG, 1.0);
          
                    while ((theElement = theDomain->theElements()) != 0) 
                      theElement->addInertiaLoadToUnbalance(*uDotDotG);
                  };
            };

            MP_ConstraintIter &theMPs = this->getMPs();
            while ((theMP = theMPs()) != 0)
              theMP->applyConstraint(newTime);
            
            SP_ConstraintIter &theSPs = this->getSPs();
            while ((theSP = theSPs()) != 0)
              theSP->applyConstraint(newTime);
            ops_Dt = dT;
        };
        if (res == 0)  res = myHandler->applyLoad() {
          // PlainHandler DOES NOTHING
        };
        if (res == 0)  res = myDomain->update() {
          ops_Dt = dT;
          ops_TheActiveDomain = this;
          while ((theEle = theEles()) != 0) {
            ops_TheActiveElement = theEle;
            ok += theEle->update();
          }
        };
        if (res == 0)  res = myHandler->update() {
          // PlainHandler DOES NOTHING
        };
      };
  };

// </div>

  if (result < 0) {
    the_Domain->revertToLastCommit();
    theIntegrator->revertToLastStep();
    return -2;
  }
  
  result = theAlgorithm->solveCurrentStep() {
      AnalysisModel   *theAnaModel = this->getAnalysisModelPtr();
      LinearSOE  *theSOE = this->getLinearSOEptr();

      if (this->theIncrIntegrator->formUnbalance() < 0) {
        return -2;
      }	    

      // set itself as the ConvergenceTest objects EquiSolnAlgo
      theTest->setEquiSolnAlgo(*this);
      if (theTest->start() < 0)
        return -3;

      int result = -1;
      numIterations = 0;

      do {
          if (tangent == INITIAL_THEN_CURRENT_TANGENT) {
            if (numIterations == 0) {
                SOLUTION_ALGORITHM_tangentFlag = INITIAL_TANGENT;
                if (theIncrIntegrator->formTangent(INITIAL_TANGENT) < 0)
                    return -1;

            } else {
              SOLUTION_ALGORITHM_tangentFlag = CURRENT_TANGENT;
              if (theIncrIntegrator->formTangent(CURRENT_TANGENT) < 0)
                return -1;
            }
          }	else { 
            SOLUTION_ALGORITHM_tangentFlag = tangent;
            if (theIncrIntegrator->formTangent(tangent, iFactor, cFactor) < 0)
                return -1;
          } 
          if (theSOE->solve() < 0)
            return -3;
          if (theIncrIntegrator->update(theSOE->getX()) < 0)
            return -4;
          if (theIncrIntegrator->formUnbalance() < 0)
            return -2;

          result = theTest->test();

          numIterations++;
          this->record(numIterations);

      } while (result == -1);

      if (result == -2)
        return -3;
      // note - if postive result we are returning what the convergence test returned
      // which should be the number of iterations 
      return result;
  };
  if (result < 0) {
    the_Domain->revertToLastCommit();            
    theIncrIntegrator->revertToLastStep();
    return -3;
  }     
  result = theIncrIntegrator->(GeneralizedAlpha::)commit() {

// <div class="function">
    theModel->setResponse(*U,*Udot,*Udotdot);
    double time = theModel->getCurrentDomainTime();
    time += (1.0-alphaF)*deltaT;
    theModel->setCurrentDomainTime(time);
    return theModel->commitDomain() {
       myDomain->commit() {
          while ((nodePtr = theNodeIter()) != 0)
            nodePtr->commitState();

          while ((elePtr = theElemIter()) != 0)
            elePtr->commitState();

          committedTime = currentTime;
          this?->dT = 0.0;

          for (int i=0; i<numRecorders; i++)
            if (theRecorders[i] != 0)
              theRecorders[i]->record(commitTag, currentTime);

          this->Domain::commitTag++;
       };
    };
  };
// </div>

  if (result < 0) {
    the_Domain->revertToLastCommit();            
    theIncrIntegrator->revertToLastStep();
    return -4;
  } 
    
  return result;
}


```
