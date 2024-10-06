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
                                                                        
// $Revision: 1.7 $
// $Date: 2009-04-30 23:23:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/OPS_Stream.h,v $

#ifndef _OPS_Stream
#define _OPS_Stream

#include <string>
#include <MovableObject.h>
enum class openMode  {OVERWRITE, APPEND};
class Vector;
class ID;

class OPS_Stream:  public MovableObject
{
 public:
  OPS_Stream(int classTag);
  virtual ~OPS_Stream();

  enum class Float { Fixed, Scientific };

  // output format
  virtual int setFile([[maybe_unused]] const char *fileName, 
                      [[maybe_unused]] openMode mode = openMode::OVERWRITE, 
                      [[maybe_unused]] bool echo = false) {return 0;}
  virtual int setPrecision([[maybe_unused]] int precision) {return 0;}
  virtual int setFloatField(Float) {return 0;}
  virtual int precision([[maybe_unused]] int precision) {return 0;}
  virtual int width([[maybe_unused]] int width) {return 0;}

  // xml stuff
  virtual int tag(const char *) =0;
  virtual int tag(const char *, const char *) =0;
  virtual int endTag() =0;
  virtual int attr(const char *name, int value) =0;
  virtual int attr(const char *name, double value) =0;
  virtual int attr(const char *name, const char *value) =0;
  virtual int write(Vector &data) =0; 
  virtual int flush();

  // regular stuff
  virtual OPS_Stream& write(const char *s, int n);
  virtual OPS_Stream& write(const unsigned char *s, int n);
  virtual OPS_Stream& write(const signed char *s, int n);
  virtual OPS_Stream& write(const void *s, int n);
  virtual OPS_Stream& write(const double *s, int n);

  virtual OPS_Stream& operator<<([[maybe_unused]] char c);
  virtual OPS_Stream& operator<<([[maybe_unused]] unsigned char c);
  virtual OPS_Stream& operator<<([[maybe_unused]] signed char c);
  virtual OPS_Stream& operator<<([[maybe_unused]] const char *s);
  virtual OPS_Stream& operator<<([[maybe_unused]] const unsigned char *s);
  virtual OPS_Stream& operator<<([[maybe_unused]] const signed char *s);
  virtual OPS_Stream& operator<<([[maybe_unused]] const void *p);
  virtual OPS_Stream& operator<<([[maybe_unused]] int n);
  virtual OPS_Stream& operator<<([[maybe_unused]] unsigned int n);
  virtual OPS_Stream& operator<<([[maybe_unused]] long n);
  virtual OPS_Stream& operator<<([[maybe_unused]] unsigned long n);
  virtual OPS_Stream& operator<<([[maybe_unused]] short n);
  virtual OPS_Stream& operator<<([[maybe_unused]] unsigned short n);
  virtual OPS_Stream& operator<<([[maybe_unused]] bool b);
  virtual OPS_Stream& operator<<([[maybe_unused]] double n);
  virtual OPS_Stream& operator<<([[maybe_unused]] float n);
  virtual OPS_Stream& operator<<([[maybe_unused]] std::string const&s) { return *this;};

  // parallel stuff
  virtual void setAddCommon(int);
  virtual int setOrder(const ID &order);
  virtual int sendSelf(int commitTag, Channel &theChannel) =0;  
  virtual int recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker) =0;

 protected:
  int addCommonFlag;

 private:
  void indent();
  int numIndent;
};

#endif
