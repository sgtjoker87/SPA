#ifndef CSPAMESHVERTEXSUPPORT_HPP
#define CSPAMESHVERTEXSUPPORT_HPP

/*
 * CSPAMeshVertexSupport.h - class defintion for Mesh Vertex Support
 *
 *
 */

/* Generated by Together */

#include "interface_api/IMeshVertexSupport.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CSPAMeshVertexSupport : public IMeshVertexSupport
{
public:

    SPA_DECLARE_NEW_DELETE( CSPAMeshVertexSupport )

    static IMeshVertexSupport* createInstance(SPA_TzMeshVertexSupport& zMeshVertexSupport );

    static IMeshVertexSupport* getInstance();

    static void deleteInstance();
      
    
    /** mesh definitions (return -1 if not defined) */
    virtual int queryElementSize( int iVertexId, double* pdElemSize );

    /** frozen status */
    virtual bool IsFrozen( int iVertexId );

    /** frozen node */
    virtual int queryFrozenNode( int iVertexId, int * piNodeId );

protected:

    CSPAMeshVertexSupport( SPA_TzMeshVertexSupport& zMeshVertexSupport );

    virtual ~CSPAMeshVertexSupport(){}

private:

    static CSPAMeshVertexSupport*   m_pzInstance;

    SPA_TzMeshVertexSupport& m_zMeshVertexSupport;

    // CSPAMeshVertexSupport( const CSPAMeshVertexSupport &X);

    CSPAMeshVertexSupport& operator= ( const CSPAMeshVertexSupport &X);

};

#endif //CSPAMESHVERTEXSUPPORT_HPP
