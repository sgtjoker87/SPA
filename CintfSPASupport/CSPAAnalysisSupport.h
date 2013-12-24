/*===========================================================================
 * CSAMAnalysisSupport.h - class definition for Analysis (ie Adaptivity) 
 *                                 Support class
 *
 * 
*/
#ifndef CSAMANALYSISSUPPORT_HPP
#define CSAMANALYSISSUPPORT_HPP

#include "interface_api/IAnalysisSupport.h"
#include "interface_api/SamInterface.h"
#include "samstlmemory/samNew.hxx"

class CSAMAnalysisSupport : public IAnalysisSupport
{
public:
	
    SAM_DECLARE_NEW_DELETE( CSAMAnalysisSupport )

    static IAnalysisSupport* createInstance(sam_TzAnalysisSupport& zAnalysisSupport);

    static IAnalysisSupport* getInstance();

    static void deleteInstance();

    // query # iterations
    virtual int queryAnalysisNumIterations( );

    /** query number of data points */
    virtual int queryAnalysisNumDataPoints();

    // (Node ID, TargetSize) assumes that the size map is based on an existing mesh stored in the workspace
    virtual int queryAnalysisDataSizeMap(int aiNodeID[], double adTargetSize[]);

    // min elt size - default 0.0 when not in use 
    virtual double queryAnalysisMinElementSize();

    // Max. number of nodes - default 0 when not in use
    virtual int queryAnalysisMaxNodes();

    // User specified adaptivity - local or global remesh desired.
    virtual bool queryAnalysisGlobalReanalyize();



protected:

    CSAMAnalysisSupport(sam_TzAnalysisSupport& zAnalysisSupport);

    virtual ~CSAMAnalysisSupport(){}

private:

    static CSAMAnalysisSupport* m_pzInstance;

    sam_TzAnalysisSupport&  m_zAnalysisSupport;

private:

    CSAMAnalysisSupport& operator=(const CSAMAnalysisSupport& X ) { return *this;}

};
#endif //CSAMAnalysisSupport_HPP
