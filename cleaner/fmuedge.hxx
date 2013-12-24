/*
*/

#ifndef FMUEDGE_HXX
#define FMUEDGE_HXX

#include <vector>
#include "samstlmemory/samSTLAllocator.hxx"

#include <cstdio>

class nodClass;
class elmBase;
class elmQClass;
class fmuMaster;

class fmuEdge
{
public:
    virtual ~fmuEdge();
    fmuEdge( nodClass*, nodClass*, fmuMaster* master=NULL );
    fmuEdge( nodClass*, nodClass*, elmBase*, fmuMaster* master=NULL );

    nodClass* start_node() const { return pCStartNode; }
    nodClass* end_node() const { return pCEndNode; }
    nodClass* other_node( nodClass *the_node ) const;
    nodClass* mid_node() const { return pCMidNode; }
    void      mid_node( nodClass* pCNode ) { pCMidNode = pCNode; }
      
    void frozen( bool new_value ) { mIsFrozen = new_value;}
    bool frozen() const { return mIsFrozen; }
      
    int number_elems() const { return iNumberElems; }
    void add_elem( elmBase* );
    void delete_elem( elmBase* );
    void elems( elmBase**, elmBase** ) const;
    void quads( elmQClass**, elmQClass** ) const;
      
    elmBase* other_elem( elmBase* ) const;
    elmQClass* other_quad( elmQClass* ) const;
    elmBase* shared_elem( fmuEdge* ) const;
    elmQClass* shared_quad( fmuEdge* ) const;

    void draw_me( int ) const;
    void highlight_me() const;

    void id( int iNewId ) { mId = iNewId; }
    int id() const { return mId; }

    int marked() const { return iMarked; }
    void marked( int mark_flag ){ iMarked = mark_flag; }

    void dump( FILE* ps ) const;

private:

    nodClass *pCStartNode;
    nodClass *pCEndNode;
    nodClass *pCMidNode;
    elmBase *apElems[2];

    int iNumberElems;
    int mId;
    int iMarked;
    bool mIsFrozen;

private:
    // disabled
    fmuEdge();
    fmuEdge( const fmuEdge& );
    fmuEdge& operator=( const fmuEdge& );
};

 typedef std::vector< fmuEdge *, SAM_STL_ALLOC( fmuEdge * ) > VEC_fmuEdge;
 typedef VEC_fmuEdge::iterator  VECITER_fmuEdge;

#endif
