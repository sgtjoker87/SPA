/*
// Subsystem Description:
//
// This is the definition of a list or dynamic array of quad elements.
// It is implemented as a macro of fmuDynArray.
//
// Class defined in this file:
//
// nodClassDynArray
//
// Subsystem History:
//
// None
//
*/
//

#ifndef ELMQARRAY_HXX
#define ELMQARRAY_HXX

#include "cleaner/fmudynarray.hxx"
class elmQClass;

typedef fmuDynArray< elmQClass* > elmQClassDynArray;

#endif
