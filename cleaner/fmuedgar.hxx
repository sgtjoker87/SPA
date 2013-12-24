/*

*/


// Subsystem Description:
//
// This is the definition of a list or dynamic array of edges.
// It is implemented as a macro of fmuDynrray.
//
// Class defined in this file:
//
// nodClassDynArray
//
// Subsystem History:
//
// None
//


#ifndef FMUEDGAR_HXX
#define FMUEDGAR_HXX

#include "cleaner/fmudynarray.hxx"
class fmuEdge;

typedef fmuDynArray< fmuEdge* > fmuEdgeDynArray;

namespace edgeDynArray
{
    void dump( FILE* ps, const fmuEdgeDynArray& dynEdges );
    void dump( int context, const fmuEdgeDynArray& dynEdges );
} // end namespace edgeDynArray

#endif
