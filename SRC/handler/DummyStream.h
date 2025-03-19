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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-08-03 23:34:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DummyStream.h,v $

#ifndef _DummyStream
#define _DummyStream

#include <OPS_Stream.h>

#include <fstream>
using std::ofstream;

class DummyStream : public OPS_Stream
{
 public:
  DummyStream();
  ~DummyStream();

  // xml stuff
  int tag(const char *) {return 0;}
  int tag(const char *, const char *) {return 0;}
  int endTag() {return 0;}
  int attr([[maybe_unused]] const char *name, [[maybe_unused]] int value) {return 0;}
  int attr([[maybe_unused]] const char *name, [[maybe_unused]] double value) {return 0;}
  int attr([[maybe_unused]] const char *name, [[maybe_unused]] const char *value) {return 0;}
  int write(Vector &data) {return 0;}


  OPS_Stream& write([[maybe_unused]] const char *s, [[maybe_unused]] int n) {return *this;}
  OPS_Stream& write([[maybe_unused]] const unsigned char *s, [[maybe_unused]] int n) {return *this;}
  OPS_Stream& write([[maybe_unused]] const signed char *s, [[maybe_unused]] int n) {return *this;}
  OPS_Stream& write([[maybe_unused]] const void *s, [[maybe_unused]] int n) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] char c) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] unsigned char c) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] signed char c) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] const char *s) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] const unsigned char *s) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] const signed char *s) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] const void *p) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] int n) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] unsigned int n) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] long n) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] unsigned long n) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] short n) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] unsigned short n) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] bool b) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] double n) {return *this;}
  OPS_Stream& operator<<([[maybe_unused]] float n) {return *this;}

  int sendSelf(int commitTag, Channel &theChannel) {return 0;};  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker) {return 0;}

 private:
};

#endif
