#ifndef FMUFHS_HXX
#define FMUFHS_HXX


#include "samcommon/portable.h"
#include "cleaner/nodclass.hxx"

#include "samstlmemory/samSTLAllocator.hxx"

#include <map>
#include <functional>


// Structure to map node label to pointer to node
struct fmu_TzLink
{
   int       iLabel;
   nodClass *pCNode;
};

typedef std::map
<
    int, 
    fmu_TzLink, 
    std::less<int>, 
    SAM_STL_MAP_ALLOC(int, fmu_TzLink) 
> 
NodeTable;

typedef NodeTable::value_type valType;

inline
int 
fhsNodePut( NodeTable* pzNodeLink, fmu_TzLink *pzNode )
{
    if ( pzNode->iLabel == 0 )
        return 1;
    pzNodeLink->insert( valType( pzNode->iLabel, (*pzNode) ) );
    return 0;
}

inline
int 
fhsNodeGet( NodeTable *pzNodeLink, fmu_TzLink *pzNode )
{
    NodeTable::iterator p = pzNodeLink->find( pzNode->iLabel );
    if ( p == pzNodeLink->end() )
        return 1;
    pzNode->pCNode = (*p).second.pCNode;
    return 0;
}

inline
int 
fhsFlushNodeTable( NodeTable *pzNodeLink )
{
    pzNodeLink->clear();
    return 0;
}

#endif
