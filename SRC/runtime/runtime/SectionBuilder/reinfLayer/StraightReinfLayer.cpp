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
// File: StraightReinfLayer.C
// Written by Remo M. de Souza
// December 1998
//
#include <Matrix.h>
#include <OPS_Stream.h>
#include <Vector.h>
#include <math.h>
#include <string>

#include <Cell.h>
#include <StraightReinfLayer.h>

StraightReinfLayer::StraightReinfLayer(int materialID, int numReinfBars, double reinfBarArea,
                                       const VectorND<2>& InitialPosition, const VectorND<2>& FinalPosition)
 : ReinfLayer(materialID, reinfBarArea),
   nReinfBars(numReinfBars),
   initPosit(InitialPosition),
   finalPosit(FinalPosition)
{
}

int
StraightReinfLayer::getNumReinfBars() const
{
  return nReinfBars;
}

std::vector<Cell> 
StraightReinfLayer::getReinfBars() const
{
  VectorND<2> barPosit;
  // ReinfBar* reinfBars;
  std::vector<Cell> bars(nReinfBars);

  if (nReinfBars == 1) {
    VectorND<2> location {
      barPosit(0) = (initPosit(0) + finalPosit(0)) / 2,
      barPosit(1) = (initPosit(1) + finalPosit(1)) / 2
    };

    // reinfBars = new ReinfBar[1];
    bars[0] = Cell(material, area, location);

    // bars[0].setPosition(barPosit);
    // bars[0].setArea(this->area);
  }

  else if (nReinfBars > 1) {
    double dy = (finalPosit(0) - initPosit(0)) / (nReinfBars - 1);
    double dz = (finalPosit(1) - initPosit(1)) / (nReinfBars - 1);

    // reinfBars = new ReinfBar[nReinfBars];

    for (int i = 0; i < nReinfBars; i++) {
      VectorND<2> location {
         initPosit(0) + dy * i,
         initPosit(1) + dz * i
      };

      bars[i] = Cell(material, area, location);
    }
  }

  return bars;
}
