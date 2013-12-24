/*
 * CSPAEdgeSupport.h - class defintion for Edge Support
 *
 */

#ifndef CSPAEDGESUPPORT_HPP
#define CSPAEDGESUPPORT_HPP

#include "interface_api/IEdgeSupport.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CSPAEdgeSupport : public IEdgeSupport 
{
public:
    
    SPA_DECLARE_NEW_DELETE( CSPAEdgeSupport )

    static IEdgeSupport* createInstance( SPA_TzEdgeSupport&  zEdgeSupport );

    static IEdgeSupport* getInstance();

    static void deleteInstance();
    
    /**
     * Query the start vertex of the input edge
     * @return int Id of the start vertex */
    virtual int queryStartVertex(int iEdgeId);

    /**
     * Query the end vertex of the input edge
     * @return int Id of the end vertex
     */
    virtual int queryEndVertex(int iEdgeId);

    virtual int queryVertices(int iEdgeId, int aiVertexId[]);


    /** geometry query -
     *  Compute the total arc length of the edge;
     *  Return the xyz location given the parametric location;
     *  Return the paramter location (SPAe as arc length from starting point)
     *  given the xyz location
     */
    virtual double queryArcLength(int iEdgeId );

    virtual int queryTangent(int iEdgeId, double dX, double adTang[]);

    virtual int queryParamLocation(int iEdgeId, double adXYZ[],
                             double* pdParamLoc);

    virtual int queryXYZLocation(int iEdgeId, double dParamLoc, double adXYZ[]);
                                           /** for discretization */

    virtual int queryCurvature(int iEdgeId, double dParamLoc, double* pdCurvature);

    // Edge representation as Control points/Polyline 
    virtual int GetNumControlPoints(int iEdgeId);

    virtual int GetXYZControlPoints(int iEdgeId, int aiPts[], double adXYZ[]);

    virtual int GetControlPointsParam(int iEdgeId,double adParam[]);

    virtual int GetLineRadiusCurvatureAtControlPoints(int iEdgeId,double adCurvature[]);

	virtual int IsSeam(int iEdgeId);

    virtual double queryFacetEdgeTolerance( int iEdgeId );

	virtual int queryAppId (int iEdgeId);

protected:

    CSPAEdgeSupport(SPA_TzEdgeSupport&  zEdgeSupport);

    virtual ~CSPAEdgeSupport(){}

private:

    static  CSPAEdgeSupport*    m_pzInstance;

    SPA_TzEdgeSupport& m_zEdgeSupport;

    // CSPAEdgeSupport( const CSPAEdgeSupport &X){}

    CSPAEdgeSupport& operator= ( 
        const CSPAEdgeSupport &X){ return *this; }

};
#endif //CSPAEDGESUPPORT_HPP




