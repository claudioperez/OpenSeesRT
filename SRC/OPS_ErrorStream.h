/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** (C) Copyright 2024, All Rights Reserved.                           **
**                                                                    **
** ****************************************************************** */
//
// Author: cmp
//
#ifndef opserr
#  include <OPS_Stream.h>
   extern OPS_Stream *opserrPtr;
#  define opserr (*opserrPtr)
#  define endln "\n"
#endif
