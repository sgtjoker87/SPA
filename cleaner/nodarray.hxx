/*

*/

/*

*/

#ifndef NODARRAY_HXX
#define NODARRAY_HXX

#include "cleaner/fmudynarray.hxx"
class nodClass;

typedef fmuDynArray< nodClass* > nodClassDynArray;
typedef fmuDynArray< nodClassDynArray* > nodClassDynArrayDynArray;

namespace nodeDynArray
{
    void dump( FILE* ps, const nodClassDynArray& dynNodes );
    void dump( int context, const nodClassDynArray& dynNodes );
} // end namespace nodeDynArray

#endif
