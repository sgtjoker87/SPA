/*===========================================================================
 * CSPAPatchSupport.h - class definition for Patch mesh surface/edge 
 *                                 Support class
 *
*/
#ifndef CSPAPATCHSUPPORT_HPP
#define CSPAPATCHSUPPORT_HPP
#include "interface_api/IPatchSupport.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CSPAPatchSupport : public IPatchSupport
{
public:
	
    SPA_DECLARE_NEW_DELETE( CSPAPatchSupport )

    static IPatchSupport* createInstance(SPA_TzPatchSupport& zPatchSupport);

    static IPatchSupport* getInstance();

    static void deleteInstance();

    /** query number of frozen element edges */
    virtual int queryPatchNumFrozenElementEdges(int iEdgeId);
        
	/** query frozen element edges */
    virtual int queryPatchFrozenElementEdges(int iEdgeId, int aiElemEdges[]);
	virtual int  queryPatchNumElements( int iFaceId );
    virtual int  queryPatchElements   ( int iFaceId, int aiElems[]);
    virtual int  queryPatchElementSize  ( double* pdPatchElemSize, double *pdPatchElementSizeFactor);
    virtual int  queryPatchTransitionFactor   ( double *pdTransitionFactor );

protected:

    CSPAPatchSupport(SPA_TzPatchSupport& zPatchSupport);

    virtual ~CSPAPatchSupport(){}

private:

    static CSPAPatchSupport* m_pzInstance;

    SPA_TzPatchSupport&  m_zPatchSupport;

private:

    CSPAPatchSupport& operator=(const CSPAPatchSupport& X ) { return *this;}

};
#endif //CSPAPATCHSUPPORT_HPP
