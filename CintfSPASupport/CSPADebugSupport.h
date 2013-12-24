#ifndef CSPADEBUGSUPPORT_HPP
#define CSPADEBUGSUPPORT_HPP

/*
 * CSPADebugSupport.h - implementation for Debug Support class
 * 
 *
 */



#include "interface_api/IDebugSupport.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CSPADebugSupport : public IDebugSupport 
{
public:
    
    SPA_DECLARE_NEW_DELETE(CSPADebugSupport)

    static IDebugSupport* createInstance( SPA_TzDebugSupport& zDebugSupport );

    static IDebugSupport* getInstance();

    static void deleteInstance();

    virtual SPA_TeDebugLevel queryDebugLevel();

    virtual void setDebugLevel(SPA_TeDebugLevel eDebugLevel);

    virtual bool debugSurfaceMesher(void);

    virtual bool debug2DMesher(void);

    virtual bool debugCleaner(void);

    virtual bool debugSmoother(void);

    virtual SPA_TeTimerOption timerOption(void);

    //display meshes
    virtual int display3dNodes( int     iEntType,
                                int     iEntId,
                                int     iNumNodes, 
                                int     aiNodeId[], 
                                double  adXYZ[] );

    virtual int display3dMesh( int      iEntType,
                               int      iEntId,
                               int      iNumNodes, 
                               int      aiNodeId[], 
                               double   adXYZ[],
                               int      iNumElem, 
                               int      aiElemId[], 
                               int      aiElemType[],
                               int      aiElemConn[] );

    virtual int display2dNodes( int     iEntType,
                                int     iEntId,
                                int     iNumNodes, 
                                int     aiNodeId[], 
                                double  adXY[] );

    virtual int display2dMesh( int      iEntType,
                               int      iEntId,
                               int      iNumNodes, 
                               int      aiNodeId[], 
                               double   adXY[],
                               int      iNumElem, 
                               int      aiElemId[], 
                               int      aiElemType[],
                               int      aiElemConn[] );

    virtual int display2dContour( int      iEntType,
                                  int      iEntId,
                                  int      iNumNodes, 
                                  int       aiNodeId[], 
                                  double    adXY[] );

    virtual int remove2dDisplay();

    virtual int remove3dDisplay();

    /** Meshing status . */
    virtual int status( const char *pcMessage );
    virtual int ProcessStatus( int *piProcessType, 
		                       int *piEntityId, 
						       int *piEntityType,
						       int *piMessageId,
						       double *pdProcessComplete 
						     );
    virtual int message( const char *pcMessage );

    virtual int prompt( const char *pcMessage , int iOption);

    virtual const char* meshFileNameRoot( int iOpt=0 ); // root file name string (iOpt: =0, file name; !=0, path name)
    virtual const char* meshSummaryFileName();          // root file name + "_mesh_summary.txt"
	virtual const char* mappedmeshSummaryFileName();    // root file name + "_mappedmesh_summary.txt"

protected:

    CSPADebugSupport( SPA_TzDebugSupport& zDebugSupport );

    virtual ~CSPADebugSupport(){};
    
private:

    static CSPADebugSupport* m_pzInstance;

    SPA_TzDebugSupport& m_zDebugSupport;
  
    // CSPADebugSupport( const CSPADebugSupport &X){}
    CSPADebugSupport& operator= ( const CSPADebugSupport& ){return *this;};
};

#endif //CSPAEDGESUPPORT_HPP
