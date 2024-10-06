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
// This class supercedes the original MaterialResponse class, so that
// it is only a thin wrapper around the GenericResponse template.
// Eventually this class should be removed altogether, and only
// GenericResponse should be used. The prior implementation of 
// MaterialResponse was the only thing that depended on having a generic
// Material base class, which was otherwise useless. Migrating to 
// GenericResponse eliminates entirely the need for Material.
//
// Written: cmp
// Created: Oct 2024
//
#ifndef MaterialResponse_h
#define MaterialResponse_h

#include <Response.h>
#include <Information.h>
#include <GenericResponse.h>

class Material;

class ID;
class Vector;
class Matrix;

#if 1

template <typename T>
class MaterialResponse : public GenericResponse<T> {
  public:
  using GenericResponse<T>::GenericResponse;
};
template <typename T, typename ...B>
MaterialResponse(T*, B...)->MaterialResponse<T>;


#else
// Original implementation:
class MaterialResponse : public Response
{
 public:
  MaterialResponse(Material *mat, int id);
  MaterialResponse(Material *mat, int id, int val);
  MaterialResponse(Material *mat, int id, double val);
  MaterialResponse(Material *mat, int id, const ID &val);
  MaterialResponse(Material *mat, int id, const Vector &val);
  MaterialResponse(Material *mat, int id, const Matrix &val);
  ~MaterialResponse();
  
  int getResponse(void);
  int getResponseSensitivity(int gradNumber);

private:
  Material *theMaterial;
  int responseID;
};
#endif

#endif // include guard

