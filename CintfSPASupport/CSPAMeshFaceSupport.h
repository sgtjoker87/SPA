/*
 * CSPAMeshSurfaceSupport.h - class definition for mesh surface Support class
 *
 */

#ifndef CSPAMESHFACESUPPORT_HPP
#define CSPAMESHFACESUPPORT_HPP
#include "interface_api/IMeshFaceSupport.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CSPAMeshFaceSupport : public IMeshFaceSupport
{
public:
	
    SPA_DECLARE_NEW_DELETE( CSPAMeshFaceSupport )

    static IMeshFaceSupport* createInstance(SPA_TzMeshFaceSupport& zMeshFaceSupport);

    static IMeshFaceSupport* getInstance();

    static void deleteInstance();

    /** mesh definitions (return -1 if not defined) */

    /** frozen status */
    virtual bool IsFrozen( int iFaceId);

    virtual int  queryNumFrozenNodes   ( int iFaceId);
    virtual int  queryFrozenNodes      ( int iFaceId, int aiNodes[], double adXYZ[]);
    virtual int  queryNumFrozenElements( int iFaceId);
    virtual int  queryFrozenElements   ( int iFaceId, int aiElems[], int aiElemTypes[], 
                                         int aiElemConn[], int *piElemOrder);

    /** mesh methods */
    virtual int  queryMeshMethod   ( int iFaceId ); // free, free-mapped, mapped
    virtual int  querySpaceStrategy( int iFaceId ); // 1: PS, 2: MAP, 3: Flattening
    virtual void querySmoothCleanOptions( int iFaceId,
		                                   int* piSmoothOption, int* piCleanOption ); // 0: off, 1: on

    virtual int  queryElementSize  ( int iFaceId, double* pdElemSize );
    virtual int  queryElementType  ( int iFaceId ); // =3, tria; =4, quad
    virtual int  queryElementOrder ( int iFaceId ); // =0, linear; =2, parabolic

    virtual int  queryStraightEdgesOption(int iFaceId);
    virtual void queryMidNodeQualityChecks( int iFaceId,
		                                    int*    piCheckJacobian,
                                            double* pdMinJacobian,
                                            int*    piCheckDeviationRatio,
                                            double* pdMaxDeviationRatio,
                                            int*    piCheckNastranEPLR );

    /** element size */
    virtual int queryNumLocalSizePts( int iFaceId );
    virtual int queryLocalSizePts   ( int iFaceId, double adXYZ[],
                                      double adLocalSize[]);

    /** hard points and hard lines */
    virtual int queryNumHardPoints(int iFaceId);
    virtual int queryHardPoints(int iFaceId, int aiHardPtIds[], double adXYZ[]);
    virtual int queryNumHardLines(int iFaceId);
    virtual int queryNumPointsOnHardLines(int iFaceId, int aiHardLineId[], int aiNumPts[]);
    virtual int queryHardLinePoints(int iFaceId, int iHardLineId, double adXYZ[]);

    // virtual int queryCylinder(int iFaceId);
    virtual int queryNumMinElmOnEdg(int iFaceId, int *piCylCircumType);
    virtual int queryNumCorners	   ( int iFaceId );
    virtual int queryCorners	   ( int iFaceId, double a3DCorners[]);
	virtual int queryElemSizeAlongCylHeight( int iFaceId, double *pdSize );
	virtual double queryElemSizeAlongFilletAxis( int iFaceId );
    virtual double queryMinimumElementSize (int iFaceId);
	virtual double queryAspectRatioLimit (int iFaceId);
    virtual double queryElemSizeAlongCircumference(int iFaceId);
	virtual int    queryElemCountAlongCircumference(int iFaceId);
    virtual int    queryFilletSelectionParams(int iFaceId, double adFilletSelectionParam[]);
	virtual bool   isFeatureFace(int iFaceId);

protected:

    CSPAMeshFaceSupport(SPA_TzMeshFaceSupport& zMeshFaceSupport);

    virtual ~CSPAMeshFaceSupport(){}

private:

    static CSPAMeshFaceSupport* m_pzInstance;

    SPA_TzMeshFaceSupport&  m_zMeshFaceSupport;

private:

    // CSPAMeshFaceSupport( const CSPAMeshFaceSupport& X ){}

    CSPAMeshFaceSupport& operator=(const CSPAMeshFaceSupport& X ) { return *this;}

};
#endif //CSPAMESHFACESUPPORT_HPP
