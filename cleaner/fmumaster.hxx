/*
*/

#ifndef FMUMASTER_HXX
#define FMUMASTER_HXX

#include "cleaner/nodarray.hxx"
#include "cleaner/fmuedgar.hxx"
#include "cleaner/elmqarray.hxx"
#include "cleaner/elmtarray.hxx"

class fmuMaster
{
public:
    ~fmuMaster(){}
    fmuMaster();

    // may modify contents of arrays
    nodClassDynArray*  global_interior(){ return &CInteriorList; }
    nodClassDynArrayDynArray* global_boundary()     { return &CBoundaryList; }
    nodClassDynArrayDynArray* original_boundary()   { return &COrigBndyList; }

    fmuEdgeDynArray*   global_edges()   { return &CEdgeList; }
    elmQClassDynArray* global_quads()   { return &CQuadList; }
    elmTClassDynArray* global_trias()   { return &CTriaList; }   

    // array contents may only be queried
    const nodClassDynArrayDynArray& global_boundary() const   { return CBoundaryList; }
    const nodClassDynArrayDynArray& original_boundary() const { return COrigBndyList; }
    const nodClassDynArray&  global_interior() const    { return CInteriorList; }
    const fmuEdgeDynArray&   global_edges() const       { return CEdgeList; }
    const elmQClassDynArray& global_quads() const       { return CQuadList; }
    const elmTClassDynArray& global_trias() const       { return CTriaList; }   

    // additional info
    void highest_node_id( int iNewValue )   { iHighestNodeId = iNewValue; }
    int  highest_node_id() const            { return iHighestNodeId; }

    void abort( bool lNewValue ){ lAbortFlag = lNewValue; }
    bool abort() const          { return lAbortFlag; }

    void triqua_mesh_flag( bool fNewFlag )  { fAllowTrias = fNewFlag; }
    bool triqua_mesh_flag() const           { return fAllowTrias; }

    void FaceId( int idFace )   { m_idFace = idFace; }
    int  FaceId() const         { return m_idFace; }

    // debugging
    void dump( int ) const;

private:
    fmuEdgeDynArray     CEdgeList;

    elmQClassDynArray   CQuadList;
    elmTClassDynArray   CTriaList;

    nodClassDynArrayDynArray    CBoundaryList;
    nodClassDynArrayDynArray    COrigBndyList;
    nodClassDynArray            CInteriorList;

    int     iHighestNodeId;
    int     m_idFace;

    bool    fAllowTrias;
    bool    lAbortFlag;

private:
    // disable
    fmuMaster( const fmuMaster& );
    fmuMaster& operator= ( const fmuMaster& );
};

#endif
