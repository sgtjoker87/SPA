#ifndef CSPAMESHOPTIONSSUPPORT_HPP
#define CSPAMESHOPTIONSSUPPORT_HPP

/*
 * CSPAMeshOptionsSupport.h - implementation for Mesh Option Support class
 *
*/

#include "interface_api/IMeshOptionsSupport.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CMeshOptions;

class CSPAMeshOptionsSupport : public IMeshOptionsSupport 
{
public:

    SPA_DECLARE_NEW_DELETE( CSPAMeshOptionsSupport )

    static IMeshOptionsSupport* createInstance(SPA_TzMeshOptionsSupport&  zMeshOptionsSupport );

    static IMeshOptionsSupport* getInstance();

    static void deleteInstance();
    
    virtual double            getElementSize(void);
    virtual SPA_TeElementShape  getElementShape(void);
    virtual SPA_TeElementOrder  getElementOrder(void);
    virtual SPA_TeElementFamily getElementFamily(void);
    virtual int               getCleanFlagOption(void);
    virtual int               getSmoothFlagOption(void);
    virtual int               getStraightEdgesFlagOption(void);
    virtual int               getDebugFlagOption();
    virtual int               getCheckMinJacobianOption();
    virtual int               getCheckMaxDeviationOption();
    virtual double            getMaxDeviationRatioOption();
    virtual double            getMinJacobianRatioOption();

    virtual int               get2DMesherOption();
    virtual int               getCleanerOption();
	virtual int               getMeshingStrategy();
    virtual int               getSmootherOption();
    virtual int               getSpaceDevelopmentOption();

	virtual double            getMinElementSizeReduction(void);
	virtual int               getSizeControlOption(void);

    virtual int               getLocalRemeshOption(void);

    virtual bool              getQualityCheckOption();
    virtual void              getQualityMeasureOptions( 
                                                        double *dMinSkew,
                                                        double *dMinWarp,
                                                        double *dMinTaper, 
                                                        double *dMinAspect
                                                      );
    virtual void initializeMapMeshingOptions();
    virtual void setMapMeshingOptions( SPA_TzMapMeshingOptions zMMO);
    virtual SPA_TzMapMeshingOptions queryMapMeshingOptions();
    virtual void mmTurnOnSmoothing();
    virtual void mmTurnOffSmoothing();
    virtual void mmTurnOnClamping();
    virtual void mmTurnOffClamping();
    virtual void mmTurnOnSubmapping();
    virtual void mmTurnOffSubmapping();
    virtual void setProcessNeed(SPA_TeMeshingProcessStage eProcessNeed);
    virtual SPA_TeMeshingProcessStage getProcessNeed();

    virtual int getHighNodeLabel();
    virtual int getHighElementLabel();

    virtual void getSplitQuadOption(int *pfSpliQuad, double *pdMaxWarp, double *pdMaxJacobian);
    virtual int  getAttemptMapMeshFlagOption(void);
    virtual void setAttemptMapMeshFlagOption(int fAMFlag);
	virtual int  getVertexMappingOption(void);
    virtual void setVertexMappingOption(int fVMFlag);
	virtual int  getFlipDiagonalsOption(void);
    virtual void setFlipDiagonalsOption(int fFDFlag);
    virtual int previewBoundary(); // 1 - yes, 0 - no
	virtual void setSmallFeatureFactor(double dSFF);
	virtual double querySmallFeatureFactor();
    virtual int    queryHigherOrderInterpolation();
    virtual double queryAbstractionTolerance(void);
    virtual void   setAbstractionTolerance(double dAbstractionTolerance);
    virtual double getAbstractionTolerance(void);
	virtual int    queryVolumeRemeshingOption();
	virtual int    queryVolumePyramidTransitionOn(void);
	virtual int    queryVolumeTwoElementsThruThicknessOn(void);
	virtual void   setMinElementSizeReduction(double dMinElementSizeReduction);
	virtual void setNumLayers(int nLayers);
	virtual int  queryNumLayers();
	virtual void initAppsCachedMeshOptions(void);
	virtual void  setVolumePyramidTransitionOn(void);
	virtual void  setVolumeTwoElementsThruThicknessOn(void);
    virtual int getAllQuadOption(void);
	virtual void setAllQuadOption(int fQuadOnly);
	virtual void  setVolumeRemeshingOption(int fVolRemeshing);
    virtual int   getVolumeRemeshingOption();
    virtual int getSmoothTargetOption(void);
	virtual void setSmoothTargetOption(int fSmoothTarget);
    virtual int  getMidNodeUpdateOption(void);
    virtual void setMidNodeUpdateOption(int fMidNodeUpdateOption);
    virtual int  queryFluidDomainSurfaceMeshOption(void);
    virtual void   setFluidDomainSurfaceMeshOption(int iFluidDomainSurfaceMeshOption);
	virtual int  queryTrgMinElmEdgLenOption(void);
	virtual void setTrgMinElmEdgLenOption(int fFlag);
	virtual int  queryFeatureAbstractionOption(void);
	virtual void setFeatureAbstractionOption(int fFlag);
	virtual double queryGlobalMinimumElementLength(void);
	virtual void setGlobalMinimumElementLength(double dMEL);
	virtual void setFeatureElemSizeFactor(double dFESF);
	virtual double queryFeatureElemSizeFactor();

protected:

    CSPAMeshOptionsSupport( SPA_TzMeshOptionsSupport&  zMeshOptionsSupport );

    virtual ~CSPAMeshOptionsSupport() {};
             

private:

    static CSPAMeshOptionsSupport*  m_pzInstance;

    SPA_TzMeshOptionsSupport& m_zMeshOptionsSupport;

    // CSPAMeshOptionsSupport( const CSPAMeshOptionsSupport& )

    CSPAMeshOptionsSupport& operator=( const CSPAMeshOptionsSupport& ){return *this;}

};
#endif //CSPAMESHOPTIONSSUPPORT_HPP

